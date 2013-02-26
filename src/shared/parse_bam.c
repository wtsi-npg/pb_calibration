#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sam.h>
#include <version.h>
#include <smalloc.h>
#include <io_lib/hash_table.h>
#include "die.h"
#include "shared.h"
#include "parse_bam.h"

/*
 * cts - parse bam file line
 *
 * returns 0 on success, 1 on expected failure.
 */
int parse_bam_readinfo( samfile_t *fp,
                    bam1_t *bam,
                    int *bam_lane,
                    int *bam_tile,
                    int *bam_x,
                    int *bam_y,
                    int *bam_read,
                    size_t *bam_offset) 
{

    char *name;
    const char *sep = ":#/";
    const char *cp;
	uint8_t *ci_ptr;
    int lane, tile, x, y, read;
	int offset;

    if (0 > samread(fp, bam)) return 1;

    lane = -1;
    tile = -1;
    x = -1;
    y = -1;

    name = bam1_qname(bam);
    cp = strchr(name,':');
    if (NULL == cp) die("ERROR: Invalid run in name: \"%s\"\n",name);
    cp++;
    cp = parse_next_int(cp,&lane,sep);
    cp = parse_next_int(cp,&tile,sep);
    cp = parse_next_int(cp,&x,sep);
    cp = parse_next_int(cp,&y,sep);

	if (bam_offset) {	/* look for offset, if we want it */
		/* look for ci tag */
		ci_ptr = bam_aux_get(bam, "ci");
		if (NULL == ci_ptr) {
			/* no ci tag get offset from name */
			cp = parse_next_int(cp,&offset,sep);
			if (NULL == cp) die("ERROR: No ci tag and no offset in name: \"%s\"\n",name);
		} else {
			offset = bam_aux2i(ci_ptr);
			/* name offset is 0 based but ci is 1 based */
			offset--;
		}
	}

	if (lane < 1) die("ERROR: Invalid lane value in name: \"%s\"\n",name);

	if (tile <= 0) die("ERROR: Invalid tile value in name: \"%s\"\n",name);

    read = 0;
    if(BAM_FPAIRED & bam->core.flag){
        if(BAM_FREAD1 & bam->core.flag)
            read = 1;
        if(BAM_FREAD2 & bam->core.flag)
            read = 2;
        if(read == 0){
            die("ERROR: Unable to determine read from flag %d for read: \"%s\"\n",bam->core.flag,name);
        }
    }

    *bam_lane = lane;
    *bam_tile = tile;
    *bam_x = x;
    *bam_y = y;
    *bam_read = read;
	if (bam_offset) *bam_offset = offset;
    return 0;
}


/*
 * cts - parse bam file line
 *
 * returns 0 on success, 1 on expected failure
 */
int
parse_bam_alignments(
             samfile_t *fp,
             bam1_t *bam,
             char *read_seq,
             int *read_qual,
             char *read_ref,
             int *read_mismatch, const int read_buff_size,
	HashTable *snp_hash)
{

    char *name;
    int32_t pos;
    uint32_t *cigar;
    uint8_t *seq, *qual, *m_ptr;
    char *mismatch = NULL;
    const char *sep = "^ACGTKMRYSWBVHDNacgtkmryswbvhdn";
    const char *cp;
    int i, j, skip;
    HashItem *hi;

    if (0 == read_buff_size)
        die("ERROR: Invalid read_buff_size");

    name = bam1_qname(bam);
    pos = bam->core.pos;
    cigar = bam1_cigar(bam);
    seq = bam1_seq(bam);
    qual = bam1_qual(bam);
    m_ptr = bam_aux_get(bam, "MD");

    if (NULL == m_ptr) {
        die("ERROR: No mismatch for read: \"%s\"\n", name);
    } else {
        mismatch = bam_aux2Z(m_ptr);
        if (NULL == mismatch) {
            die("ERROR: Invalid mismatch %s for read: \"%s\"\n",
                mismatch, name);
        }
    }

    memset(read_mismatch, 0, bam->core.l_qseq * sizeof(int));

    for (i = 0; i < bam->core.l_qseq; i++) {
        read_seq[i] = bam_nt16_rev_table[bam1_seqi(seq, i)];
        read_qual[i] = qual[i];
    }
    read_seq[i] = 0;
	if (read_ref) read_ref[i] = 0;

    j = 0;
    for (i = 0; i < bam->core.n_cigar; i++) {
        int l = cigar[i] >> 4, op = cigar[i] & 0xf, k;
        switch (op) {
        case BAM_CMATCH:
            // CIGAR: alignment match;
            for (k = 0; k < l; j++, k++) {
                read_mismatch[j] |= BASE_ALIGN;
				if (read_ref) read_ref[j] = read_seq[j];
			}
            break;
        case BAM_CDEL:
            // CIGAR: deletion from the reference
            if (j == bam->core.l_qseq) {
                die("ERROR: Trailing deletion for read: %s\n",
                    name);
            }
            read_mismatch[j] |= BASE_DELETION;
            break;
        case BAM_CINS:
            // CIGAR: insertion to the reference 
            for (k = 0; k < l; j++, k++) {
                read_mismatch[j] |= BASE_INSERTION;
				if (read_ref) read_ref[j] = 'I';
			}
            break;
        case BAM_CSOFT_CLIP:
            // CIGAR: clip on the read with clipped sequence present in qseq
            for (k = 0; k < l; j++, k++) {
                read_mismatch[j] |= BASE_SOFT_CLIP;
				if (read_ref) read_ref[j] = 'S';
			}
            break;
        default:
            die("ERROR: Unexpected CIGAR operation: %c\n", op);
        }
    }
    if (j != bam->core.l_qseq) {
        die("ERROR: Inconsistent cigar string %d > %d for read: \"%s\"\n", j, bam->core.l_qseq, name);
    }

    /* clipped sequence or insertions are missing from MD */
    for (i = 0, skip = 0; i < bam->core.l_qseq; i++, skip++)
      if (0 == (read_mismatch[i] & (BASE_SOFT_CLIP | BASE_INSERTION)))
            break;

    cp = mismatch;
    while (NULL != (cp = parse_next_int(cp, &j, sep))) {
        /* skip matching bases, exclude insertions which are missing from MD */
        for (; j > 0; i++)
            if (read_mismatch[i] & BASE_INSERTION)
                skip++;
            else
                j--;

        if (0 == strlen(cp))
            /* reached end of MD string */
            break;

        /* skip insertions which are missing from MD */
        for (; i < bam->core.l_qseq; i++, skip++)
            if (0 == (read_mismatch[i] & BASE_INSERTION))
                break;
        if (i == bam->core.l_qseq) {
            die("ERROR: Invalid MD string %s for read: \"%s\"\n",
                mismatch, name);
        }

        switch (*cp) {
        case '^':
            /* deletions missing from read_seq */
            while(NULL != strchr(sep, *(++cp)))
                skip--;
            break;
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            /* mismatch */
            if (0 == (read_mismatch[i] & BASE_ALIGN)) {
                die("ERROR: Inconsistent cigar string expect alignment at mismatch for read: \"%s\"\n", name);
            }
			if (read_ref) read_ref[i] = *cp;
            read_mismatch[i] |= BASE_MISMATCH;

            pos = bam->core.pos + (i - skip);

            if (NULL != snp_hash) {
                char *chrom =
                    fp->header->target_name[bam->core.tid];
                char key[100];
                /* N.B bam->core.pos is 0 based */
                snprintf(key, sizeof(key), "%s:%d", chrom, pos);
                if (NULL !=
                    (hi =
                     HashTableSearch(snp_hash, key,
                             strlen(key)))) {
                    hi->data.i++;
                    read_mismatch[i] |= BASE_KNOWN_SNP;
                }
            }
            i++;
            break;
        default:
            /* treat all other reference bases as a known SNP */
            while(NULL != strchr(sep, *cp)) 
            {
                read_mismatch[i++] |= BASE_KNOWN_SNP;
                cp++;
            }
            break;
        }
    }

    /* clipped sequence or insertions are missing from MD */
    for (; i < bam->core.l_qseq; i++)
        if (0 == (read_mismatch[i] & (BASE_SOFT_CLIP | BASE_INSERTION)))
            break;
    if (i != bam->core.l_qseq) {
        die("ERROR: Inconsistent MD string %d != %d for read: \"%s\"\n",
            i, bam->core.l_qseq, name);
    }

    if (BAM_FREVERSE & bam->core.flag) {
        read_seq = reverse_seq(read_seq);
        read_seq = complement_seq(read_seq);
        read_qual = reverse_int(read_qual, bam->core.l_qseq);
		if (read_ref) {
			read_ref = reverse_seq(read_ref);
			read_ref = complement_seq(read_ref);
		}
        read_mismatch = reverse_int(read_mismatch, bam->core.l_qseq);
    }

    return 0;
}

/*
 * bam_header_add_pg - add PG record to BAM header
 */

static void append_header_text(bam_header_t *header, char* text, int len)
{
	int x = header->l_text + 1;
	int y = header->l_text + len + 1; // 1 byte null
	if (text == 0) return;
#if 0
	/* this is only sensible if l_text was rounded up when text was first allocated
           We call append_header_text() with a header generated by bam_header_dup()
           which doesn't round up l_text so we should not round up */
        kroundup32(x); 
	kroundup32(y);
#endif
	if (x < y) header->text = (char*)realloc(header->text, y);
	strncpy(header->text + header->l_text, text, len); // we cannot use strcpy() here.
	header->l_text += len;
	header->text[header->l_text] = 0;
}

void bam_header_add_pg(char *id, char *pn, char *ds, char *cl, bam_header_t *bam_header)
{
	char *pp = NULL;
	char *id2, *pg;
	int id2size, id2num, pgsize, pglen;

	if (NULL == bam_header) {
		die("ERROR: No bam header\n");
	}

	if (NULL == bam_header->text) {
		die("ERROR: No text in bam header\n");
	}

        id2size = 10 + strlen(id);
        id2 = smalloc(id2size);
        id2num = 0;
        while(id2num < 10) {
                char *text = strdup(bam_header->text);
                char *hl, *endl, *endt;
                sprintf(id2, (id2num > 0 ? "%s%d" : "%s"), id, id2num);
                hl = text;
                while (0 < strlen(hl)) {
                        if (NULL == (endl = strchr(hl, '\n'))) {
                                die("ERROR: Corrupt bam header \"%s\"\n", hl);
                        }
                        *endl = 0;
                        if (0 == memcmp(hl, "@PG", 3)) {
                                if (NULL == (pp = strstr(hl, "ID:"))) {
                                        die("ERROR: No ID in PG line \"%s\"\n", hl);
                                }
                                pp += 3;
                                if (NULL != (endt = strchr(pp, '\t'))) {
                                        *endt = 0;
                                }
                                if (0 == strcmp(pp, id2)) {
                                        /* an existing ID matches the new ID */
                                        break;
                                }
                        }
                        hl = endl + 1;
                }
                if( 0 == strlen(hl)) {
                        /* the new ID doesn't match an existing ID */
                        free(text);
                        break;
                }
                /* the new ID matches an existing ID, increment the number and try again */
                id2num++;
                free(text);
        }
        if (id2num == 10) {
            die("ERROR: Header already contains PG lines with ID=%s .. ID=%s\n", id, id2);
        }

	pgsize = 128 + strlen(pn) + strlen(pn) + strlen(ds) + strlen(version) + strlen(cl);
	if (NULL != pp) {
		pgsize += strlen(pp);
		pg = smalloc(pgsize);
		pglen =
		    snprintf(pg, pgsize,
			     "@PG\tID:%s\tPN:%s\tPP:%s\tDS:%s\tVN:%s\tCL:%s\n", id2, pn, pp, ds, version, cl);
	} else {
		pg = smalloc(pgsize);
		pglen =
		    snprintf(pg, pgsize,
			     "@PG\tID:%s\tPN:%s\tDS:%s\tVN:%s\tCL:%s\n", id2, pn, ds, version, cl);
	}
	assert(pglen < pgsize);

	append_header_text(bam_header, pg, pglen);

	free(pg);
}
