/*  -*- mode: c; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 8; -*- */

/* TO-DO:
 *
 * Argument parsing needs tidying up. It should use optarg maybe.
 */

/*
 * Copyright (c) 2006-2009, Genome Research Ltd (GRL).
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *
 *     * Neither the name of the Genome Research Limited nor the
 *       names of its contributors may be used to endorse or promote
 *       products derived from this software without specific prior
 *       written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GRL ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL GRL BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Parts of this code were inherited from sequenceread, see
 * http://sequenceread.svn.sourceforge.net/viewvc/sequenceread, which
 * contained code edited by Illumina. As no illumina code is used here
 * the Illumina copyright notice has been removed
 */

/*
 * Author: Steven Leonard, Jan 2009
 *
 * This code applies a purity based calibration table
 *
 */

#ifdef HAVE_CONFIG_H
#include "pb_config.h"
#endif

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <fcntl.h>
#include <math.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <glob.h>
#include <regex.h>
#include <stdarg.h>

/* To turn off assert for a small speed gain, uncomment this line */
/* #define NDEBUG */

/* Hack to stop io_lib from trying to include its own config.h */
#ifdef HAVE_CONFIG_H
#undef HAVE_CONFIG_H
#endif
#include <io_lib/misc.h>
#include <io_lib/hash_table.h>

#include <sam.h>

#include <smalloc.h>
#include <aprintf.h>
#include <die.h>
#include <rts.h>

#include <version.h>

#define MAX_CIF_CHUNK_BYTES 4194304

/* if we split data by state, using a filter file rather than tile, use tile as a place holder for state */
#define N_STATES 2

#define PHRED_QUAL_OFFSET 33  // phred quality values offset
#define ILL_QUAL_OFFSET   64  // illumina quality values offset

typedef struct {
    int         tile;
    int         read;
    int         cycle;
    float       offset;
    float       delta;
    float       scale;
    int         nbins;
    float       *purity;
    long        *num_bases;
    long        *num_errors;
    float       *frac_bases;
    float       *error_rate;
    float       *quality;
    int         *ibin;
} CalTable;

typedef struct {
    long        lane;
    long        cycle;
    long        read;
    const char *dir;
} CifDir;

typedef struct {
    char *cmdline;
    char *intensity_dir;
    char *output;
    int read_length[3];
    int cstart[3];
    int quiet;
    int spatial_filter;
    int region_size;
    int nregions_x;
    int nregions_y;
    char *working_dir;
    CifDir *cif_dirs;
    size_t  n_cif_dirs;
    size_t *cif_lane_index;
    size_t  n_cif_lanes;
} Settings;

typedef union {
    int8_t  *i8;
    int16_t *i16;
    int32_t *i32;
} ChanData;

typedef struct {
    char *filename;
    int fd;
    int data_type;
    int num_channels;
    int cycle_number;
    off_t  cycle_start_pos;
    size_t num_entries;
    size_t chunk_num_entries;
    size_t chunk_start;
    ChanData *chan_data;
} CifCycleData;

typedef struct {
    CifCycleData *cycles;
    size_t        ncycles;
    size_t        size;
    size_t        num_spots;
} CifData;


static
void
checked_chdir(const char* dir){
    if(chdir(dir)){
        fprintf(stderr,"ERROR: failed to change directory to: %s\n",dir);
        exit(EXIT_FAILURE);
    }
}

static void initialiseCalTable(CalTable *ct, int tile, int read, int cycle)
{
    int nbins = 76;
    int nvals = 76;
    int ival;

    ct->tile = tile;
    ct->read = read;
    ct->cycle = cycle;

    ct->offset = 0.25;
    ct->delta  = 0.01;
    ct->scale  = 1.0 / ct->delta;

    ct->purity = smalloc(nbins * sizeof(float));
    ct->num_bases = smalloc(nbins * sizeof(long));
    ct->num_errors = smalloc(nbins * sizeof(long));
    ct->frac_bases = smalloc(nbins * sizeof(float));
    ct->error_rate = smalloc(nbins * sizeof(float));
    ct->quality = smalloc(nbins * sizeof(float));

    ct->ibin = smalloc(nvals * sizeof(int));
    for(ival=0;ival<nvals;ival++)
        ct->ibin[ival] = -1;

    ct->nbins = 0;

    return;
}

static void freeCalTable(CalTable *ct)
{
    if(ct->nbins) {
        free(ct->purity);
        free(ct->num_bases);
        free(ct->num_errors);
        free(ct->frac_bases);
        free(ct->error_rate);
        free(ct->quality);
        free(ct->ibin);
    }
    free(ct);
}

static int restoreCalTable(Settings *s, const char* calibrationFile, HashTable *ct_hash)
{
    int nct = 0;
    static const int line_size = 1028;
    char line[line_size];
    FILE *fp;
    CalTable *current_ct = NULL;

    fp = fopen(calibrationFile, "r");
    if( NULL == fp ){
        fprintf(stderr, "ERROR: can't open CT file %s: %s\n",
                calibrationFile, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if( !s->quiet )
        fprintf(stderr, "reading calibration table %s\n", calibrationFile);

    while( fgets(line, line_size, fp) ){
        int tile, read, cycle, num_bases, num_errors;
        float purity, frac_bases, error_rate, quality;
        char key[100];
        HashItem *hi;
        HashData hd;
        int ibin;
        
        if( 9 != sscanf(line, "%f %d %d %d %d %d %f %f %f", &purity, &read, &cycle, &tile,
                        &num_bases, &num_errors, &frac_bases, &error_rate, &quality) ){
            perror("restoreCalTable");
            exit(EXIT_FAILURE);
        }

        if( read < 0 || read >= N_READS ){
            fprintf(stderr,"ERROR: Invalid read in CT file (%s) %d < 0 or > %d.\n",
                    calibrationFile, read, N_READS);
            exit(EXIT_FAILURE);
        }

        if( cycle < 0 ){
            fprintf(stderr,"ERROR: Invalid cycle in CT file (%s) %d < 0.\n",
                    calibrationFile, cycle);
            exit(EXIT_FAILURE);
        }

        if( s->spatial_filter ){
            if( tile < 0 || tile >= N_STATES ){
                fprintf(stderr,"ERROR: Invalid state in CT file (%s) %d < 0 or > %d.\n",
                        calibrationFile, tile, N_STATES);
                exit(EXIT_FAILURE);
            }
        }
    
        snprintf(key, sizeof(key), "%d:%d:%d", tile, read, cycle);
        if( NULL == (hi = HashTableSearch(ct_hash, key, strlen(key))) ){
            hd.p = smalloc(sizeof(CalTable));
            if( NULL == HashTableAdd(ct_hash, key, strlen(key), hd, NULL) ) {
                fprintf(stderr, "ERROR: building tile cycle hash table\n");
                exit(EXIT_FAILURE);
            }

            if( s->spatial_filter ){
                if( tile == 1 ) fprintf(stderr, "bad region read=%d cycle=%3d\n", read, cycle);
            }else{
                if( tile < 0 ) fprintf(stderr, "bad tile=%4d read=%d cycle=%3d\n", tile, read, cycle);
            }
            
            /* current ct is the new ct */
            current_ct = (CalTable *)hd.p;

            initialiseCalTable(current_ct, tile, read, cycle);

            nct++;
        }else{
            current_ct = (CalTable *)hi->data.p;
        }

        if( current_ct->nbins == 76 ){
            fprintf(stderr, "ERROR: number of lines in CT file (%s) for tile/state=%d read=%d cycle=%d exceed maximum %d > %d\n",
                    calibrationFile, current_ct->tile, current_ct->read, current_ct->cycle, current_ct->nbins, 76);
            exit(EXIT_FAILURE);
        }

        ibin = current_ct->nbins;

        current_ct->purity[ibin] = purity;

        current_ct->num_bases[ibin] = num_bases;
        current_ct->num_errors[ibin] = num_errors;
        current_ct->frac_bases[ibin] = frac_bases;
        current_ct->error_rate[ibin] = error_rate;
        current_ct->quality[ibin] = quality;

        current_ct->nbins++;
    }

    fclose(fp);

    if( 0 == nct ){
        fprintf(stderr, "ERROR: no rows in calibration table %s\n",
                calibrationFile);
        exit(EXIT_FAILURE);
    }

    return nct;
}

////////////////////////////////////////////////////
// Calculates the purity
// data - array of intensities
// return purity
////////////////////////////////////////////////////
float GetPu (int n, int data[])
{
    int Isum, Imin, Imax;
    int i;

    Isum = Imin = Imax = data[0];
    for (i=1; i<n; i++)
    {
        if(Imin > data[i])
            Imin=data[i];
        if( Imax < data[i])
            Imax = data[i];
        Isum += data[i];
    }
    if( Imin <= 0 ) {
        Imin--;
        Imax -= Imin;
        Isum -= 4 * Imin;
    }

    return (float)Imax/(float)Isum;
}

int Get_bin_purity(CalTable *ct, int value)
{
    int ibin = 0;

    if (ct->ibin[value] < 0) {
        float purity = ct->offset + ct->delta * value;
        for(ibin=0;ibin<ct->nbins;ibin++)
            if (purity<=ct->purity[ibin])
                break;
        if (ibin == ct->nbins)
            ibin--;
        ct->ibin[value] = ibin;
    }

    ibin = ct->ibin[value];
    return ibin;
}

/**
 * Reverses the direction of the array of ints
 *
 * @param int is the input array. <b>NB.</b> this is destructively modified.
 *
 * @returns a pointer to the original storage but with the contents
 * modified
 */
int *
reverse_int(int *num, int n)
{
    int *t = num, *s = num + n - 1;
    int c;

    while (t < s) {
        c = *s;
        *s = *t;
        *t = c;
        t++;
        s--;
    }

    return num;
}


// lifted from /nfs/users/nfs_j/jkb/work/tracetools/sequtil/seq_ops.c

/**
 * Reverses the direction of the string.
 *
 * @param seq is the input sequence. <b>NB.</b> this is destructively modified.
 *
 * @returns a pointer to the original string storage but with the contents
 * modified
 */
char *
reverse_seq(char *seq)
{
    char *t = seq, *s = seq + strlen(seq) - 1;
    char c;

    while (t < s) {
        c = *s;
        *s = *t;
        *t = c;
        t++;
        s--;
    }

    return seq;
}

/* cts -simplification of parse_4_int code for single int parse
 */
inline static
const char *
parse_next_int(const char *str, int *val, const char *sep) {
    const static char spaces[256] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    const char* const start = str;
    int minus = 0;
    int ival = 0;
    char c;

    if (NULL == sep) {
        while (*str && spaces[(unsigned) *str]) ++str;
    } else {
        while (*str && NULL != strchr(sep, *str)) ++str;
    }

    c = *str;

    if (!c) {
      /*
        fprintf(stderr, "Error: expected to parse int from string \"%s\"\n", start);
        exit(EXIT_FAILURE);
	*/
      return NULL;
    }

    if (c == '-' || c == '+') {
        minus = (c == '-');
        c = *++str;
    }
    
    while (c >= '0' && c <= '9') {
        ival = ival * 10 + (c-'0');
        c = *++str;
    }

    if (NULL == sep) {
        switch(c) {
        case '\n': case '\r': case '\t': case ' ': case '\0':
            if (minus) ival = -ival;
            *val = ival;
            return str;
        }
    } else {
        if (NULL != strchr(sep, *str)) {
            if (minus) ival = -ival;
            *val = ival;
            return str;
        }
    }
    fprintf(stderr, "Error: expected to parse int from string \"%s\"\n", start);
    exit(EXIT_FAILURE);
}


/*
 * cts - parse bam file line
 *
 * returns 0 on success, 1 on expected failure.
 */
static
int
parse_bam_file_line(Settings *s,
                    samfile_t *fp,
                    bam1_t *bam,
                    int *bam_lane,
                    int *bam_tile,
                    int *bam_x,
                    int *bam_y,
                    size_t *bam_offset,
                    int *bam_read) {

    char *name;
    uint8_t *ci_ptr;
    const char *sep = ":#/";
    const char *cp;
    int lane, tile, x, y, offset, read;

    if( 0 > samread(fp, bam)) return 1;
    
    lane = -1;
    tile = -1;
    x = -1;
    y = -1;
    offset = -1;

    name = bam1_qname(bam);
    cp = strchr(name,':');
    if(NULL == cp){
        fprintf(stderr,"ERROR: Invalid run in name: \"%s\"\n",name);
        exit(EXIT_FAILURE);
    }
    cp++;
    cp = parse_next_int(cp,&lane,sep);
    cp = parse_next_int(cp,&tile,sep);
    cp = parse_next_int(cp,&x,sep);
    cp = parse_next_int(cp,&y,sep);

    /* look for ci tag */
    ci_ptr = bam_aux_get(bam, "ci");
    if (NULL == ci_ptr){
        /* no ci tag get offset from name */
        cp = parse_next_int(cp,&offset,sep);
        if(NULL == cp){
            fprintf(stderr,"ERROR: No ci tag and no offset in name: \"%s\"\n",name);
            exit(EXIT_FAILURE);
        }
    }else{
        offset = bam_aux2i(ci_ptr);
        /* name offset is 0 based but ci is 1 based */
        offset--;
    }
       
    if(lane > N_LANES+1 || lane < 1){
        fprintf(stderr,"ERROR: Invalid lane value in name: \"%s\"\n",name);
        exit(EXIT_FAILURE);
    }

    if(tile <= 0){
        fprintf(stderr,"ERROR: Invalid tile value in name: \"%s\"\n",name);
        exit(EXIT_FAILURE);
    }

    if(offset < 0){
        fprintf(stderr,"ERROR: Invalid offset value in name: \"%s\"\n",name);
        exit(EXIT_FAILURE);
    }

    read = 0;
    if(BAM_FPAIRED & bam->core.flag){
        if(BAM_FREAD1 & bam->core.flag)
            read = 1;
        if(BAM_FREAD2 & bam->core.flag)
            read = 2;
        if(read == 0){
            fprintf(stderr,"ERROR: Unable to determine read from flag %d for read: \"%s\"\n",bam->core.flag,name);
            exit(EXIT_FAILURE);
        }
    }

    *bam_lane = lane;
    *bam_tile = tile;
    *bam_x = x;
    *bam_y = y;
    *bam_offset = offset;
    *bam_read = read;

    return 0;
}


/* copied from sam.c */

static bam_header_t *bam_header_dup(const bam_header_t *h0)
{
	bam_header_t *h;
	int i;
	h = bam_header_init();
	*h = *h0;
	h->hash = h->dict = h->rg2lib = 0;
	h->text = (char*)calloc(h->l_text + 1, 1);
	memcpy(h->text, h0->text, h->l_text);
	h->target_len = (uint32_t*)calloc(h->n_targets, 4);
	h->target_name = (char**)calloc(h->n_targets, sizeof(void*));
	for (i = 0; i < h->n_targets; ++i) {
		h->target_len[i] = h0->target_len[i];
		h->target_name[i] = strdup(h0->target_name[i]);
	}
	return h;
}
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

static void bam_header_add_pg(Settings *s, bam_header_t *bam_header) {
    char *text;
    char *hl, *endl, *endt;
    char *id = "pb_cal";
    char *pn = "predictor_pu";
    char *pp = NULL;
    char *ds = "A program to apply a calibration table";
    char *pg;
    int pgsize, pglen;

    if (NULL == bam_header) {
        fprintf(stderr,"ERROR: No bam header\n");
        exit(EXIT_FAILURE);
    }

    if (NULL == bam_header->text) {
        fprintf(stderr,"ERROR: No text in bam header\n");
        exit(EXIT_FAILURE);
    }

    text = strdup(bam_header->text);

    hl = text;
    while (0 < strlen(hl)) {
        if (NULL == (endl = strchr(hl, '\n'))) {
            fprintf(stderr,"ERROR: Corrupt bam header \"%s\"\n", hl);
            exit(EXIT_FAILURE);
        }
        *endl = 0;
        if (0 == memcmp(hl, "@PG", 3)) {
            if (NULL == (pp = strstr(hl, "PN:"))) {
                fprintf(stderr,"ERROR: No ID in PG line \"%s\"\n", hl);
                exit(EXIT_FAILURE);
            }
            pp += 3;
            if (NULL != (endt = strchr(pp, '\t'))) {
                *endt = 0;
            }
        }
        hl = endl + 1;
    }
    
    pgsize = 128 + strlen(pn) + strlen(pn) + strlen(ds) + strlen(version) + strlen(s->cmdline);
    if( NULL != pp ){
        pgsize += strlen(pp);
        pg = smalloc(pgsize);
        pglen = snprintf(pg, pgsize, "@PG\tID:%s\tPN:%s\tPP:%s\tDS:%s\tVN:%s\tCL:%s\n", id, pn, pp, ds, version, s->cmdline);
    }else{
        pg = smalloc(pgsize);
        pglen = snprintf(pg, pgsize, "@PG\tID:%s\tPN:%s\tDS:%s\tVN:%s\tCL:%s\n", pn, pn, ds, version, s->cmdline);
    }
    assert(pglen < pgsize);
    
    append_header_text(bam_header, pg, pglen);

    free(text);
    free(pg);
}

int cif_dir_compare(const void *va, const void *vb) {
    const CifDir *a = (const CifDir *) va;
    const CifDir *b = (const CifDir *) vb;
    
    if (a->lane  < b->lane)  return -1;
    if (a->lane  > b->lane)  return  1;
    if (a->read  < b->read)  return -1;
    if (a->read  > b->read)  return  1;
    if (a->cycle < b->cycle) return -1;
    if (a->cycle > b->cycle) return  1;
    return strcmp(a->dir, b->dir);
}

static int get_cif_dirs(Settings *s) {
    char   *pattern = NULL;
    size_t  pattern_sz = 0;
    glob_t  glob_buf;
    CifDir *cif_dirs;
    size_t *lane_index;
    regex_t    re;
    regmatch_t matches[4];
    size_t  ndirs;
    size_t  i, l;
    long    max_lane = 0;
    int     res;

    pattern_sz = strlen(s->intensity_dir) + 30;
    pattern    = smalloc(pattern_sz);

    snprintf(pattern, pattern_sz,
             "%s/L*/C*", s->intensity_dir);

    memset(&glob_buf, 0, sizeof(glob_buf));
    memset(&re,       0, sizeof(re));

    res = glob(pattern,
               GLOB_ERR|GLOB_NOSORT|GLOB_NOESCAPE, NULL, &glob_buf);
    if (0 != res) {
        if (GLOB_NOMATCH != res) {
            perror("Looking for CIF directories");
        }
        globfree(&glob_buf);
        free(pattern);
        return -1;
    }
    
    cif_dirs = smalloc(glob_buf.gl_pathc * sizeof(CifDir));

    res = regcomp(&re, "/L([[:digit:]]+)/C([[:digit:]]+)\\.([[:digit:]]+)/*$",
                  REG_EXTENDED);
    if (0 != res) goto regerr;
    
    for (ndirs = 0, i = 0; i < glob_buf.gl_pathc; i++) {
        const char *p = glob_buf.gl_pathv[i];

        res = regexec(&re, p, sizeof(matches)/sizeof(matches[0]), matches, 0);
        if (REG_NOMATCH == res) continue;
        if (0 != res) goto regerr;
        
        cif_dirs[ndirs].lane  = strtol(p + matches[1].rm_so, NULL, 10);
        cif_dirs[ndirs].cycle = strtol(p + matches[2].rm_so, NULL, 10);
        cif_dirs[ndirs].read  = strtol(p + matches[3].rm_so, NULL, 10);
        cif_dirs[ndirs].dir   = strdup(p);
        if (NULL == cif_dirs[ndirs].dir) goto nomem;
        if (cif_dirs[ndirs].lane > max_lane) max_lane = cif_dirs[ndirs].lane;
        ndirs++;
    }
    regfree(&re);
    globfree(&glob_buf);

    if (0 == ndirs) {
        free(cif_dirs);
        return -1;
    }

    qsort(cif_dirs, ndirs, sizeof(*cif_dirs), cif_dir_compare);

    lane_index = scalloc(max_lane + 2, sizeof(*lane_index));

    for (i = 0, l = 0; i < ndirs; i++) {
        if (cif_dirs[i].lane > l) {
            long j;

            for (j = l + 1; j <= cif_dirs[i].lane; j++) {
                lane_index[j] = i;
            }
            l = cif_dirs[i].lane;
        }
    }
    lane_index[l + 1] = ndirs;

    s->cif_dirs       = cif_dirs;
    s->n_cif_dirs     = ndirs;
    s->cif_lane_index = lane_index;
    s->n_cif_lanes    = max_lane;

    return 0;

 nomem:
    perror("get_cif_dirs");
    exit(EXIT_FAILURE);

 regerr:
    {
        char msg[1024];
        regerror(res, &re, msg, sizeof(msg));
        fprintf(stderr,
                "Regular expression search failed in get_cif_dirs: %s\n",
                msg);
        exit(EXIT_FAILURE);
    }
}

/* Read count bytes at position offset in file fd.  It actually emulates
   pread(2) as the real thing doesn't seem to be any faster, and it
   may not be present. */
static ssize_t pread_bytes(int fd, void *buf, size_t count, off_t offset) {
    char *b = (char *) buf;
    ssize_t res;
    ssize_t total = 0;
    
    if (lseek(fd, offset, SEEK_SET) < 0) return -1;

    do {
        do {
            res = read(fd, b + total, count - total);
        } while (res < 0 && EINTR == errno);
        if (res > 0) total += res;
    } while (res > 0 && total < count);
    return res < 0 ? res : total;
}

static void read_cif_chunk(Settings *s, CifCycleData *cycle, size_t spot_num) {
    size_t num_entries;
    size_t num_bytes;
    size_t chan_size_bytes;
    ssize_t got;
    off_t  file_pos;
    int chan;
    
    assert(spot_num < cycle->num_entries);
    if (spot_num >= cycle->chunk_start
        && spot_num < cycle->chunk_start + cycle->chunk_num_entries) {
        return;
    }

    chan_size_bytes  = cycle->num_entries * cycle->data_type;
    num_entries  = (cycle->num_entries - spot_num < cycle->chunk_num_entries
                    ? cycle->num_entries - spot_num
                    : cycle->chunk_num_entries);
    num_bytes = num_entries * cycle->data_type;

    if (NULL == cycle->chan_data) {
        size_t chunk_size_bytes = cycle->chunk_num_entries * cycle->data_type;
        cycle->chan_data = smalloc(cycle->num_channels * sizeof(cycle->chan_data[0]));
        for (chan = 0; chan < cycle->num_channels; chan++) {
            cycle->chan_data[chan].i8 = smalloc(chunk_size_bytes);
        }
    }

    for (chan = 0; chan < cycle->num_channels; chan++) {
        file_pos = (cycle->cycle_start_pos
                    + chan_size_bytes * chan
                    + spot_num * cycle->data_type);
        got = pread_bytes(cycle->fd, cycle->chan_data[chan].i8, num_bytes, file_pos);
        if (got < 0) {
            die("Error reading %s: %s\n", cycle->filename, strerror(errno));
        }
        if (got < num_bytes) {
            die("Error: Did not get the expected amount of data from %s\n",
                cycle->filename);
        }
    }
    cycle->chunk_start = spot_num;

#ifdef WORDS_BIGENDIAN
    /* CIF files are little-endian, so we need to switch everything
       around on big-endian machines */
    switch (cycle->data_type) {
    case 1:
        break;
    case 2:
        for (chan = 0; chan < cycle->num_channels; chan++) {
            size_t i;
            int16_t *d = cycle->chan_data[chan].i16;
            for (i = 0; i < num_entries; i++) {
                d[i] = bswap_16(d[i]);
            }
        }
        break;
    case 4:
        for (chan = 0; chan < cycle->num_channels; chan++) {
            size_t i;
            int32_t *d = cycle->chan_data[chan].i32;
            for (i = 0; i < num_entries; i++) {
                d[i] = bswap_32(d[i]);
            }
        }
        break;
    default:
        die("Unexpected data type in CIF file %s\n", cycle->filename);
    }
#endif
}

static int read_cif_file(char *name, int fd, CifData *cif_data) {
    uint8_t cif_header[13];
    int    num_cycles;
    int    first_cycle;
    int    num_channels = 4;
    size_t num_entries;
    size_t chunk_num_entries;
    size_t cycle_size_bytes;
    off_t  cycle_start_pos = sizeof(cif_header);
    int    data_type;
    CifCycleData *cycle;
    size_t i;
    
    if (sizeof(cif_header)
        != pread_bytes(fd, cif_header, sizeof(cif_header), 0)) return -1;
    if (0 != memcmp(cif_header, "CIF", 3))                 return -1;
    if (cif_header[3] != '\1')                             return -1;

    data_type   = cif_header[4];
    first_cycle = cif_header[5] | (cif_header[6] << 8);
    num_cycles  = cif_header[7] | (cif_header[8] << 8);
    num_entries = (cif_header[9]
                   | (cif_header[10] << 8)
                   | (cif_header[11] << 16)
                   | (cif_header[12] << 24));

    if (0 == cif_data->num_spots) {
        cif_data->num_spots = num_entries;
    } else if (cif_data->num_spots != num_entries) {
        fprintf(stderr, "Got unexpected number of entries in CIF file %s\n"
                "Expected: %zd; Got %zd\n",
                name, cif_data->num_spots, num_entries);
    }
    if (data_type < 1 || data_type > 4) {
        die("Unexpected data_type in CIF file %s\n", name);
    }

    chunk_num_entries = MAX_CIF_CHUNK_BYTES / (data_type * num_channels);
    if (chunk_num_entries > num_entries) chunk_num_entries = num_entries;
    cycle_size_bytes = num_entries * data_type * num_channels;

    for (i = 0; i < num_cycles; i++) {
        if (cif_data->ncycles == cif_data->size) {
            cif_data->size = cif_data->size > 0 ? cif_data->size * 2 : 128;
            cif_data->cycles = srealloc(cif_data->cycles,
                                        cif_data->size * sizeof(CifCycleData));
        }
        
        cycle = cif_data->cycles + cif_data->ncycles;
        
        cycle->filename     = sstrdup(name);
        cycle->fd           = fd;
        cycle->data_type    = data_type;
        cycle->num_channels = num_channels;
        cycle->cycle_number = i + first_cycle;
        cycle->cycle_start_pos = cycle_start_pos + i * cycle_size_bytes;
        cycle->num_entries  = num_entries;
        cycle->chunk_num_entries = chunk_num_entries;
        cycle->chunk_start  = num_entries; /* Will force first load */
        cycle->chan_data    = NULL; /* Allocate memory on first load */

        cif_data->ncycles++;
    }

    return 0;
}

static CifData *load_cif_data(Settings *s, int lane, int tile, char *suffix) {
    char *cif_file = NULL;
    size_t cif_file_sz = strlen(s->intensity_dir) + strlen(suffix) + 100;
    size_t first_dir;
    size_t end_dir;
    size_t i;
    int cif = -1;
    CifData *cif_data = NULL;
    
    cif_file = smalloc(cif_file_sz);

    cif_data = scalloc(1, sizeof(CifData));

    first_dir = (lane <= s->n_cif_lanes
                 ? s->cif_lane_index[lane]
                 : s->cif_lane_index[s->n_cif_lanes + 1]);
    end_dir = (lane <= s->n_cif_lanes
               ? s->cif_lane_index[lane + 1]
               : s->cif_lane_index[s->n_cif_lanes + 1]);

    for (i = first_dir; i < end_dir; i++) {
        snprintf(cif_file, cif_file_sz,
                 "%s/s_%d_%d.%s", s->cif_dirs[i].dir, lane, tile, suffix);
        cif = open(cif_file, O_RDONLY);
        if (cif < 0) {
            if (ENOENT == errno) continue;
            die("Couldn't open %s : %s\n", cif_file, strerror(errno));
        }

        if (0 != read_cif_file(cif_file, cif, cif_data)) {
            die("Error reading %s\n", cif_file);
        }
    }

    free(cif_file);

    /* Check we have a full set of cycles */
    for (i = 0; i < cif_data->ncycles; i++) {
        if (cif_data->cycles[i].cycle_number > i + 1) {
            die("Error: Missing cycle %zd for lane %d tile %d"
                " from CIF files.\n", i + 1, lane, tile);
        } else if (cif_data->cycles[i].cycle_number != i + 1) {
            die("Error: Unexpected cycle number %d for lane %d tile %d"
                " from CIF files.\n", cif_data->cycles[i].cycle_number,
                lane, tile);
        }
    }

    return cif_data;
}

static void free_cif_data(CifData *cif_data) {
    size_t i;
    int c;

    if (NULL != cif_data->cycles) {
        for (i = 0; i < cif_data->ncycles; i++) {
            if (NULL != cif_data->cycles[i].chan_data) {
                for (c = 0; c < cif_data->cycles[i].num_channels; c++) {
                    if (NULL != cif_data->cycles[i].chan_data[c].i8)
                        free(cif_data->cycles[i].chan_data[c].i8);
                }
                free(cif_data->cycles[i].chan_data);
            }
            if (0 < cif_data->cycles[i].fd) {
                if (0 != close(cif_data->cycles[i].fd)) {
                    if (NULL != cif_data->cycles[i].filename)
                        die("Error when closing %s: %s\n", cif_data->cycles[i].filename, strerror(errno));
                    else
                        die("Error when closing cif file: %s\n", strerror(errno));
                }
            }
            if (NULL != cif_data->cycles[i].filename)
                free(cif_data->cycles[i].filename);
        }
        free(cif_data->cycles);
    }
    
    free(cif_data);
}

/*
  Calculate purity, update quality value and write bam file

  Returns 0 on success
          1 on EOF
  Quits the program if the data in the files is inconsistent.
*/

static int update_bam_qualities(Settings *s, HashTable *ct_hash, CifData *cif_data,
                                size_t spot_num, int tile, int x, int y, int read,
                                samfile_t *fp_bam, bam1_t *bam) {
    int cstart = s->cstart[read];
    int read_length = s->read_length[read];
    int iregion = -1;
    int c, b;
    uint8_t *qual;
    static const int read_buff_size = 1024;
    uint8_t oq[read_buff_size];
    int read_qual[read_buff_size];

    assert((cstart + read_length) <= cif_data->ncycles);

    if (s->spatial_filter) iregion = xy2region(x, y, s->region_size, s->nregions_x, s->nregions_y);

    /* copy original qualities to OQ:Z tag, N.B. have to null terminate aux strings */
    qual = bam1_qual(bam);
    for (b = 0; b < read_length; b++)
        oq[b] = qual[b] + PHRED_QUAL_OFFSET;
    oq[b] = 0;
    bam_aux_append(bam, "OQ", 'Z', read_length+1, oq);

    /* calculate purity and update qualities */
    for (c = cstart, b = 0; b < read_length; c++, b++) {
        CifCycleData *cycle = cif_data->cycles + c;
        char key[100];
        HashItem *hi;
        CalTable *ct;
        int channel;
        int bin[cycle->num_channels];
        float purity = -1;
        int value, ibin, quality;

        /* set cycle ct */
        if (s->spatial_filter) {
            /* state/read/cycle ct */
            int state = (getFilterData(tile, read, b, iregion) & REGION_STATE_MISMATCH) ? 1 : 0;
            snprintf(key, sizeof(key), "%d:%d:%d", state, read, b);
            if (NULL == (hi = HashTableSearch(ct_hash, key, strlen(key)))) {
                fprintf(stderr,"ERROR: no calibration table for state=%d read=%d cycle=%d.\n", state, read, b);
                exit(EXIT_FAILURE);
            }
        }else{
            /* tile/read/cycle ct */
            snprintf(key, sizeof(key), "%d:%d:%d", tile, read, b);
            if (NULL == (hi = HashTableSearch(ct_hash, key, strlen(key)))) {
                /* read/cycle ct */
                snprintf(key, sizeof(key), "%d:%d:%d", -1, read, b);
                if (NULL == (hi = HashTableSearch(ct_hash, key, strlen(key)))) {
                    fprintf(stderr,"ERROR: no calibration table for tile=%d read=%d cycle=%d.\n", tile, read, b);
                    exit(EXIT_FAILURE);
                }
            }
        }
        ct = (CalTable *)hi->data.p;

        read_cif_chunk(s, cycle, spot_num);
        for (channel = 0; channel < cycle->num_channels; channel++) {
            size_t base_pos = spot_num - cycle->chunk_start;
            switch (cycle->data_type) {
            case 1:
                bin[channel] = cycle->chan_data[channel].i8[base_pos];
                break;
            case 2:
                bin[channel] = cycle->chan_data[channel].i16[base_pos];
                break;
            case 4:
                bin[channel] = cycle->chan_data[channel].i32[base_pos];
                if (bin[channel] > 65535) {
                    bin[channel] = 65535;
                } else if (bin[channel] < -65535) {
                    bin[channel] = -65535;
                }
                break;
            default:
                abort();
            }
        }

        purity = GetPu(cycle->num_channels, bin);

        value = ct->scale * (purity - ct->offset) + 0.5;

        ibin = Get_bin_purity(ct, value);
        quality = (int)(ct->quality[ibin] + 0.5);
        read_qual[b] = quality;
    }

    if (BAM_FREVERSE & bam->core.flag)
        reverse_int(read_qual, read_length);

    qual = bam1_qual(bam);
    for (b = 0; b < read_length; b++)
        qual[b] = read_qual[b];

    /* write bam file */
    if( 0 > samwrite(fp_bam, bam)) {
        fprintf(stderr, "Error: writing bam file\n");
        exit(EXIT_FAILURE);
    }

    return 0;
}

/*
 * Takes a bam file as input and outputs a re-calibrated bam file
 * It uses the associated CIF files too to do this.
 *
 * Assumption: within a single input file, all reads are the same length and
 * we're using unclipped data. We then utilise this to construct a common
 * header to reduce overhead of data.
 *
 * Returns: 0 written for success
 *	   -1 for failure
 */
int recalibrate_bam(Settings *s, HashTable *ct_hash, samfile_t *fp_in_bam, samfile_t *fp_out_bam, size_t *nreads) {
    CifData *cif_data = NULL;

    HashTable *tile_hash;

    int ncycles_firecrest = -1;

    int lane = -1;
    int tile = -1;

    int ntiles_bam = 0;
    size_t nreads_bam = 0;

    tile_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE | HASH_FUNC_JENKINS3);
    if (!tile_hash) die("Failed to create tile_hash\n");

    bam1_t *bam = bam_init1();

    checked_chdir(s->intensity_dir);

    /* loop over reads in the input bam file */
    while (1){
        int bam_lane = -1, bam_tile = -1, bam_x = -1, bam_y = -1, bam_read = -1, read_length;
        size_t bam_offset = 0;

        if (0 != parse_bam_file_line(s, fp_in_bam, bam,
                                     &bam_lane, &bam_tile, &bam_x, &bam_y, &bam_offset, &bam_read)) {
            break;
        }

        read_length = bam->core.l_qseq;
        if (0 == read_length) continue;

        if (0 == s->read_length[bam_read]) {
            s->read_length[bam_read] = read_length;
        }
        if (s->read_length[bam_read] != read_length) {
            fprintf(stderr,
                    "Error: inconsistent read lengths "
                    "within bam file for read %d.\n"
                    "have length %ld, previously it was %d.\n",
                    bam_read, (long) read_length, s->read_length[bam_read]);
            exit(EXIT_FAILURE);
        }

        if (lane == -1) {
            lane = bam_lane;
        }
        if (bam_lane != lane){
            fprintf(stderr,"Error: Inconsistent lane.\n");
            fprintf(stderr,"Bam lane %d qseq lane %d.\n",bam_lane, lane);
            exit(EXIT_FAILURE);
        }

        if (bam_tile != tile) {
            tile = bam_tile;
            
            if (!s->quiet) fprintf(stderr, "Processing tile %i (%lu)\n", tile, nreads_bam);

            /* Look for processed trace data */
            if (NULL != cif_data) free_cif_data(cif_data);
            cif_data = load_cif_data(s, lane, tile, "dif");

            /* Check that we actually got some trace data */
            if (NULL == cif_data) {
                fprintf(stderr, "Error: no intensity files found for lane %i tile %i.\n", lane, tile);
                exit(EXIT_FAILURE);
            }

            if(ncycles_firecrest == -1) {
                ncycles_firecrest = cif_data->ncycles;
            } else if(cif_data->ncycles != ncycles_firecrest){
                fprintf(stderr,
                        "ERROR: %lu intensity cycles for tile %i"
                        "with %i cycles expected.\n",
                        cif_data->ncycles, tile, ncycles_firecrest);
                exit(EXIT_FAILURE);
            }

            // lookup itile from tile in tile hash
            HashItem *tileItem = HashTableSearch(tile_hash, (char *)&tile, sizeof(tile));
            if (tileItem) {
                fprintf(stderr,"ERROR: alignments are not sorted by tile.\n");
                exit(EXIT_FAILURE);
            } else {
  	         HashData hd;
                 hd.i = ntiles_bam++;
                 if (HashTableAdd(tile_hash, (char *)&tile, sizeof(tile), hd, NULL) == NULL)
                     die("Failed to add tile %d to tile_hash\n", bam_tile);
            }
        }

        if (0 != update_bam_qualities(s, ct_hash, cif_data,
                                      bam_offset, bam_tile, bam_x, bam_y, bam_read,
                                      fp_out_bam, bam)) {
            fprintf(stderr,"ERROR: updating quality values.\n");
            exit(EXIT_FAILURE);
        }

        nreads_bam++;
    }
    
    *nreads = nreads_bam;

    bam_destroy1(bam);
    
    if (NULL != cif_data) free_cif_data(cif_data);

    return 0;
}

static char * alloc_getcwd(void) {
    size_t sz = 1024;
    char *out = smalloc(sz);
    
    while (NULL == getcwd(out, sz)) {
        if (ERANGE != errno) {
            free(out);
            return NULL;
        }

        sz *= 2;
        out = srealloc(out, sz);
    }

    return out;
}

/*
 * Get the absolute file path and (depending on the libc) the real
 * name for sym-link dirs in path.
 * Safe for in_path == out_path case.
 */
static char * get_real_path_name(const char* in_path) {
    char   *oldwd;
    char   *out_path;

    oldwd = alloc_getcwd();
    if (NULL == oldwd) return NULL;

    checked_chdir(in_path);

    out_path = alloc_getcwd();
    
    checked_chdir(oldwd);
    free(oldwd);

    return out_path;
}



static
void usage(int code) {
    FILE* usagefp = stderr;

    fprintf(usagefp, "pb_predictor v%s\n\n", version);
    fprintf(usagefp, 
            "Usage: pb_predictor [options] bam file\n"
            ""
            "  outputs a bam file with updated quality values\n"
            "");
    fprintf(usagefp, "  options:\n");
    fprintf(usagefp, "    -ct file\n");
    fprintf(usagefp, "             calibration table to apply\n");
    fprintf(usagefp, "    -o output\n");
    fprintf(usagefp, "             Output bam file name\n");
    fprintf(usagefp, "               no default will write output to stdout\n");
    fprintf(usagefp, "    -u\n");
    fprintf(usagefp, "             do not compress the output bam file\n");
    fprintf(usagefp, "               default compress\n");
    fprintf(usagefp, "    -intensity-dir dir\n");
    fprintf(usagefp, "             Intensity directory\n");
    fprintf(usagefp, "               no default\n");
    fprintf(usagefp, "    -filter_file file\n");
    fprintf(usagefp, "               spatial filter file\n");
    fprintf(usagefp, "    -cstart int\n");
    fprintf(usagefp, "             intensity cycle number of first base of read in single-end bam file\n");
    fprintf(usagefp, "               no default\n");
    fprintf(usagefp, "    -cstart1 int\n");
    fprintf(usagefp, "             intensity cycle number of first base of read 1 in paired-end bam file\n");
    fprintf(usagefp, "               no default\n");
    fprintf(usagefp, "    -cstart2 int\n");
    fprintf(usagefp, "             intensity cycle number of first base of read 2 in paired-end bam file\n");
    fprintf(usagefp, "               no default\n");
    fprintf(usagefp, 
            "    -q       Quiet \n");
    fprintf(usagefp, "\n");

    exit(code);
}



/* arg parsing helper function
 */
static
void
check_arg(const int i,
          const int argc,
          const char* option){
    if ((i+1) >= argc) {
        fprintf(stderr, "no argument for %s option.\n",option);
        usage(1);
    }
}


char *get_command_line(int argc, char **argv) {
    char *cmdline = NULL;
    size_t sz = argc; /* All the spaces & the terminating \0 */
    size_t pos = 0;
    int i;

    for (i = 0; i < argc; i++) {
        sz += strlen(argv[i]);
    }

    cmdline = smalloc(sz);

    for (i = 0; i < argc && pos < sz; i++) {
        int len = snprintf(cmdline + pos, sz - pos,
                           i > 0 ? " %s" : "%s", argv[i]);
        if (len < 0) {
            perror("snprintf");
            exit(EXIT_FAILURE);
        }
        pos += len;
    }
    assert(pos < sz);
    
    return cmdline;
}

int main(int argc, char **argv) {

    Settings settings;
    int i;
    int bam_compress = 1;
    char *in_bam_file;
    samfile_t *fp_input_bam;
    char out_mode[5] = "wb";
    bam_header_t *out_bam_header = NULL;
    char *out_bam_file = NULL;
    samfile_t *fp_output_bam;
    char *ct_filename = NULL;
    const char *override_intensity_dir = NULL;
    const char *filter_file = NULL;
    size_t nreads = 0;
    HashTable *ct_hash = NULL;
    int nct = 0;

    settings.output = NULL;
    settings.quiet = 0;
    settings.spatial_filter = 0;
    settings.region_size = 0;
    settings.nregions_x = 0;
    settings.nregions_y = 0;
    settings.cstart[0] = 0;
    settings.cstart[1] = 0;
    settings.cstart[2] = 0;
    settings.intensity_dir = NULL;
    settings.read_length[0] = 0;
    settings.read_length[1] = 0;
    settings.read_length[2] = 0;
    settings.working_dir = NULL;
    settings.cif_dirs       = NULL;
    settings.n_cif_dirs     = 0;
    settings.cif_lane_index = NULL;
    settings.n_cif_lanes    = 0;

    settings.cmdline = get_command_line(argc, argv);

    /* Parse args */
    for (i = 1; i < argc && argv[i][0] == '-'; i++) {
	if (!strcmp(argv[i], "-")) {
	    break;
	} else if (!strcmp(argv[i], "-intensity-dir")) {
            if(override_intensity_dir != NULL) {
		fprintf(stderr, "ERROR: -intensity-dir option specified multiple times\n");
                usage(1);
            }
            check_arg(i,argc,"-intensity-dir");
            override_intensity_dir = argv[++i];

	} else if (!strcmp(argv[i], "-filter_file")) {
            if(filter_file != NULL) {
		fprintf(stderr, "ERROR: -filter_file specified multiple times\n");
                usage(1);
            }
            check_arg(i,argc,"-filter_file");
            filter_file = argv[++i];

	} else if (!strcmp(argv[i], "-q")) {
	    settings.quiet = 1;

	} else if (!strcmp(argv[i], "-u")) {
            bam_compress = 0;
	} else if (!strcmp(argv[i], "-o")) {
            if(settings.output != NULL) {
		fprintf(stderr, "ERROR: -o option specified multiple times\n");
                usage(1);
            }
            check_arg(i,argc,"-o");
            settings.output = argv[++i];

	} else if (!strcmp(argv[i], "-ct")){
            check_arg(i,argc,"-ct");
            ct_filename = argv[++i];
	} else if (!strcmp(argv[i], "-cstart1")){
            check_arg(i,argc,"-cstart1");
            settings.cstart[1] = atoi(argv[++i]);
            if(settings.cstart[1] < 1){
                fprintf(stderr,"ERROR: invalid argument to -cstart1\n");
                usage(1);
            }
            /* cycles are indexed from 0 not 1 */
            --settings.cstart[1];
	} else if (!strcmp(argv[i], "-cstart2")){
            check_arg(i,argc,"-cstart2");
            settings.cstart[2] = atoi(argv[++i]);
            if(settings.cstart[2] < 1){
                fprintf(stderr,"ERROR: invalid argument to -cstart2\n");
                usage(1);
            }
            /* cycles are indexed from 0 not 1 */
            --settings.cstart[2];
	} else if (!strcmp(argv[i], "-cstart")){
            check_arg(i,argc,"-cstart");
            settings.cstart[0] = atoi(argv[++i]);
            if(settings.cstart[0] < 1){
                fprintf(stderr,"ERROR: invalid argument to -cstart\n");
                usage(1);
            }
            /* cycles are indexed from 0 not 1 */
            --settings.cstart[0];

	} else if (!strcmp(argv[i], "-h")) {
	    usage(0);
	} else {
            fprintf(stderr,"ERROR: Unknown option %s\n", argv[i]);
	    usage(1);
	}
    }

    if ((argc-i) < 1)
	usage(0);

    in_bam_file = argv[i++];

    /* preserve starting directory b/c recalibrate is going to chdir all over the place */
    settings.working_dir = alloc_getcwd();
    if (NULL == settings.working_dir) {
        fprintf(stderr, "ERROR: can't obtain working directory: %s\n",
                strerror(errno));
        exit(EXIT_FAILURE);
    }

    /* get absolute intensity dir*/
    if(override_intensity_dir){
        settings.intensity_dir = get_real_path_name(override_intensity_dir);
        if (NULL == settings.intensity_dir) {
            fprintf(stderr, "ERROR: can't process intensity dir: %s\n",
                    override_intensity_dir);
            exit(EXIT_FAILURE);
        }
    } else {
        fprintf(stderr,"ERROR: you must specify an intensity dir\n");
        exit(EXIT_FAILURE);
    }

    /* read filter file */
    if (NULL != filter_file) {
        FILE *fp = fopen(filter_file, "rb");
        if (!fp) die("Can't open filter file %s\n", filter_file);
        Header filter_header;
        readHeader(fp, &filter_header);
        readFilterData(fp, &filter_header);
        settings.spatial_filter = 1;
        settings.region_size = filter_header.region_size;
        settings.nregions_x = filter_header.nregions_x;
        settings.nregions_y = filter_header.nregions_y;
    }

    if (NULL == (ct_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE|HASH_FUNC_JENKINS3))) {
        fprintf(stderr, "ERROR: creating tile caltable hash\n");
        exit(EXIT_FAILURE);
    }

    // read the callibration table
    nct = restoreCalTable(&settings, ct_filename, ct_hash);

    /* Look for CIF directories */
    get_cif_dirs(&settings);

    /* open the input bam file */
    fp_input_bam = samopen(in_bam_file, "rb", 0);
    if (NULL == fp_input_bam) {
        fprintf(stderr, "ERROR: can't open bam file file %s: %s\n",
                in_bam_file, strerror(errno));
        exit(EXIT_FAILURE);
    }

    /* open the output bam file */
    if (0 == bam_compress)
        strcat(out_mode, "u");

    out_bam_header = bam_header_dup(fp_input_bam->header);
    bam_header_add_pg(&settings, out_bam_header);

    out_bam_file = (NULL == settings.output ? aprintf("/dev/stdout") : aprintf("%s/%s", settings.working_dir, settings.output));
    fp_output_bam = samopen(out_bam_file, out_mode, out_bam_header);
    if (NULL == fp_output_bam) {
        fprintf(stderr, "ERROR: can't open bam file file %s: %s\n",
                out_bam_file, strerror(errno));
        exit(EXIT_FAILURE);
    }
    free(out_bam_file);

    bam_header_destroy(out_bam_header);

    if (-1 == recalibrate_bam(&settings, ct_hash, fp_input_bam, fp_output_bam, &nreads)) {
        fprintf(stderr,"ERROR: failed to process bam file %s\n", in_bam_file);
        exit(EXIT_FAILURE);
    }

    samclose(fp_input_bam);
    samclose(fp_output_bam);

    if (!settings.quiet) {
        fprintf(stderr, "Wrote %lu traces\n", nreads);
    }

    if (NULL != ct_hash) {
	HashIter *tileIter = HashTableIterCreate();
	HashItem *hashItem;
	while ((hashItem = HashTableIterNext(ct_hash, tileIter)))
            freeCalTable((CalTable *)hashItem->data.p);
    }

    if (NULL != settings.cif_dirs) free(settings.cif_dirs);

    if (NULL != settings.working_dir) free(settings.working_dir);

    return EXIT_SUCCESS;

}
