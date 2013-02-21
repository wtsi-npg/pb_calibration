/*  -*- mode: c; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 8; -*- */

/* TO-DO:
 *
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
 * Author: Steven Leonard, Jan 2009
 *
 * This code looks for spatial features given an aligned bam file
 *
 */

/*
 * various compile time options, either uncomment option
 * or compile with CFLAGS = '-D$OPTION ..'
 *
 * QC_FAIL
 *   ignore reads that fail QC when building calibration table
 *
 * PROPERLY_PAIRED
 *   only use properly paired reads when building calibration table
 *
 * BINARY_FILTER_FILE
 *   write filter file in binary format
 */

#define QC_FAIL
#define PROPERLY_PAIRED
#define BINARY_FILTER_FILE

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
#include <stdarg.h>
#include <getopt.h>

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
#include <shared.h>
#include <snp.h>
#include <parse_bam.h>

#include <version.h>

#define PHRED_QUAL_OFFSET  33  // phred quality values offset

#define REGION_MIN_COUNT            122    // minimum coverage when setting region state

// assuming a region size of 700 and ~125 reads in a region, a threshold of 0.016 equates to 2 reads in a region
#define REGION_MISMATCH_THRESHOLD   0.016  // threshold for setting region mismatch state
#define REGION_INSERTION_THRESHOLD  0.016  // threshold for setting region insertion state
#define REGION_DELETION_THRESHOLD   0.016  // threshold for setting region deletion state

#define TILE_REGION_THRESHOLD  0.75  // threshold for setting region state at tile level

#define REGION_STATE_MASK  (REGION_STATE_INSERTION | REGION_STATE_DELETION)  // region mask used to filter reads

typedef struct {
	char *cmdline;
	char *filter;
	char *intensity_dir;
	char *snp_file;
	char *in_bam_file;
	HashTable *snp_hash;
	HashTable *tile_hash;
	char *working_dir;
	char *output;
	int lane;
	int read_length[3];
	int dump;
	int calculate;
	int apply;
	int qcfail;
	int tileViz;
	int quiet;
	int region_size;
	int nregions_x;
	int nregions_y;
	int nregions;
	int compress;
	int width;
	int height;
	int region_min_count;
	float region_mismatch_threshold;
	float region_insertion_threshold;
	float region_deletion_threshold;
} Settings;

static void setRegions(Settings * s)
{
	// if an intensity directory is given read the ImageSize.dat file to get the tile size
	if (s->intensity_dir) {
                char *file = NULL;
                size_t file_sz = strlen(s->intensity_dir) + 30;
                int fd = -1;
                ssize_t res;
                uint32_t width, height;
		file = smalloc(file_sz);
		snprintf(file, file_sz, "%s/../ImageSize.dat", s->intensity_dir);
		fd = open(file, O_RDONLY);
		if (fd < 0) die("Couldn't open %s : %s\n", file, strerror(errno));
		res = read(fd, &width, 4);
		if (4 != res) die("Error reading %s\n", file);
		res = read(fd, &height, 4);
		if (4 != res) die("Error reading %s\n", file);
		s->width = width;
		s->height = height;
   	        if (!s->quiet) display("Read tile width=%u height=%u from file %s\n", s->width, s->height, file);
		close(fd);
	        free(file);
	}

	s->nregions_x = 1 + (int)(s->width / s->region_size);
	s->nregions_y = 1 + (int)(s->height / s->region_size);

	s->nregions = s->nregions_x * s->nregions_y;

	if (!s->quiet) display("nregions_x=%d nregions_y=%d nregions=%d\n", s->nregions_x, s->nregions_y, s->nregions);

	return;
}

//
// generate tileviz like summary plots for read number tileViz
//
static void displayTileViz(Settings *s, int ntiles, int tile)
{
	HashItem *hashItem;
	int itile;
	int read = s->tileViz;
	int *coverage  = smalloc(s->nregions * sizeof(int));
	int *mismatch  = smalloc(s->nregions * sizeof(int));
	int *insertion = smalloc(s->nregions * sizeof(int));
	int *deletion  = smalloc(s->nregions * sizeof(int));
        int cycle, iregion, i, j;
	hashItem = HashTableSearch(s->tile_hash, (char *)&tile, sizeof(tile));
	if (!hashItem) die("Can't find tile %d in tile_hash\n", tile);
	itile = hashItem->data.i;
 
        for (iregion = 0; iregion < s->nregions; iregion++) {
            coverage[iregion] = 0;
            mismatch[iregion]  = 0;
            insertion[iregion] = 0;
            deletion[iregion]  = 0;
        }

        for (cycle = 0; cycle < s->read_length[read]; cycle++) {
                for (iregion = 0; iregion < s->nregions; iregion++) {
                        RegionTable *rt = getRTS(itile,iregion,read,cycle);
                        int n = rt->align + rt->insertion + rt->deletion + rt->soft_clip + rt->known_snp;
                        // correct for sparse bins by assuming ALL bins have atleast REGION_MIN_COUNT clusters
 	                n = max(n, s->region_min_count);
                        // coverage should be the same for all cycles
                        coverage[iregion] = n;
                        // for all other values take the maximum over all cycles
                        mismatch[iregion]  = max(mismatch[iregion],  rt->mismatch);
                        insertion[iregion] = max(insertion[iregion], rt->insertion);
                        deletion[iregion]  = max(deletion[iregion],  rt->deletion);
    	        }
	}

	display ("tile=%-4d coverage       mismatch       insertion      deletion\n", tile);
	for (i = 0; i < s->nregions_y; i++) {
		display("          ");
                // truncate coverage to 15 so it displays as a single hex value
		for (j = 0; j < s->nregions_x; j++) {
			iregion = j * s->nregions_y + i;
			int v = min(coverage[iregion], 15);
			display("%1x", v);
		}
		display("    ");
		// convert all other values to a percentage and bin 0(0%), 1(10%), .. 10(100%)
		// correct for sparse bins by assuming ALL bins have atleast REGION_MIN_COUNT clusters
		for (j = 0; j < s->nregions_x; j++) {
			iregion = j * s->nregions_y + i;
			int v = (int)(10.0 * mismatch[iregion] / coverage[iregion]);
			display("%1x", v);
		}
		display("    ");
		for (j = 0; j < s->nregions_x; j++) {
			iregion = j * s->nregions_y + i;
			int v = (int)(10.0 * insertion[iregion] / coverage[iregion]);
			display("%1x", v);
		}
		display("    ");
		for (j = 0; j < s->nregions_x; j++) {
			iregion = j * s->nregions_y + i;
			int v = (int)(10.0 * deletion[iregion] / coverage[iregion]);
			display("%1x", v);
		}
		display("\n");
	}
	display("\n");

        free(coverage);
        free(mismatch);
        free(insertion);
        free(deletion);
}

/*
 * identify bad tiles quality < mean quality - filter * stdev quality
 * filter is typically 2.
*/

static void findBadRegions(Settings *s, int ntiles)
{
	HashItem *hashItem;

	if (0 >= ntiles)
		return;

	HashIter *iter = HashTableIterCreate();
	while ((hashItem = HashTableIterNext(s->tile_hash, iter))) {
		int tile = *(int*)(hashItem->key);
		int itile = hashItem->data.i;

		if (s->tileViz) displayTileViz(s,ntiles,tile);

		int read;
		for (read = 0; read < N_READS; read++) {
			int cycle;
			for (cycle = 0; cycle < s->read_length[read]; cycle++) {
				// set the state for each region
				int iregion;
				for (iregion = 0; iregion < s->nregions; iregion++) {
					RegionTable *rt = getRTS(itile,iregion,read,cycle);
					rt->state = 0;

   					// coverage
					int n = rt->align + rt->insertion + rt->deletion + rt->soft_clip + rt->known_snp;

					// coverage - mark spare bins
					if (n < s->region_min_count) rt->state |= REGION_STATE_COVERAGE;

					// correct for sparse bins by assuming ALL bins have atleast REGION_MIN_COUNT clusters
					n = max(n, s->region_min_count);

					// mismatch - mark bins with maximum mismatch rate > threshold
					if (((float)rt->mismatch  / (float)n) >= s->region_mismatch_threshold)  rt->state |= REGION_STATE_MISMATCH;
					// insertion - mark bins with maximum insertion rate > threshold
					if (((float)rt->insertion / (float)n) >= s->region_insertion_threshold) rt->state |= REGION_STATE_INSERTION;
					// deletion - mark bins with maximum deletion rate > threshold
					if (((float)rt->deletion  / (float)n) >= s->region_deletion_threshold)  rt->state |= REGION_STATE_DELETION;
					//if (rt->state) display("%d:%d:%d:%d:%d:%d\t%02x\n", s->lane, tileArray[itile], iregion, read, cycle, n, rt->state);
				}
				// ignoring low coverage, if all regions with a non-zero state have the same state and the
                                // fraction of regions with this state exceeds a theshold set the state for the whole tile
				int tile_state = -1, nregions = 0;
				for (iregion = 0; iregion < s->nregions; iregion++) {
					RegionTable *rt = getRTS(itile,iregion,read,cycle);
					int state = rt->state & ~REGION_STATE_COVERAGE;
					if (!state) continue;
					if (tile_state == -1) tile_state = state;
					if (state != tile_state) break;
					nregions++;
				}
				if (iregion == s->nregions && (((float)nregions/(float)s->nregions) >= TILE_REGION_THRESHOLD))
                                        for (iregion = 0; iregion < s->nregions; iregion++) {
                                                RegionTable *rt = getRTS(itile,iregion,read,cycle);
                                                rt->state = tile_state | (rt->state & REGION_STATE_COVERAGE);
                                        }
			        if (!s->quiet) {
                                        int mismatch = 0, insertion = 0, deletion = 0, soft_clip = 0;
                                        long quality_bases = 0, quality_errors = 0;
                                        for (iregion = 0; iregion < s->nregions; iregion++) {
                                                RegionTable *rt = getRTS(itile,iregion,read,cycle);
                                                if (rt->state & REGION_STATE_MISMATCH)  mismatch++;
                                                if (rt->state & REGION_STATE_INSERTION) insertion++;
                                                if (rt->state & REGION_STATE_DELETION)  deletion++;
                                                if (rt->state & REGION_STATE_SOFT_CLIP) soft_clip++;
                                                quality_bases  += rt->align;
                                                quality_errors += rt->mismatch;
                                        }
                                        float ssc = 1.0;
                                        float quality = -10.0 * log10((quality_errors + ssc)/(quality_bases + ssc));
                                        display("tile=%-4d read=%1d cycle=%-3d quality=%.2f mismatch=%-4d insertion=%-4d deletion=%-4d soft_clip=%-4d\n",
                                                tile, read, cycle, quality, mismatch, insertion, deletion, soft_clip);
                                }
			}
		}
	}

	return;
}

void printFilter(Settings *s, int ntiles, FILE *fp) 
{
	int read, cycle, iregion;
	Header hdr;
	HashItem *hashItem;
	int *tileArray;

	// create the tile array from the tile hash
	tileArray = smalloc(ntiles * sizeof(int));
	HashIter *tileIter = HashTableIterCreate();
	int i = 0;
	while ((hashItem = HashTableIterNext(s->tile_hash, tileIter))) {
		if (i>=ntiles) die("Too many tiles!\n");
		tileArray[i++] = *(int*)hashItem->key;
	}

#ifdef BINARY_FILTER_FILE
	hdr.region_magic = strdup(REGION_MAGIC);
	hdr.coord_shift = COORD_SHIFT; 
	hdr.coord_factor = COORD_FACTOR;
	hdr.ntiles = ntiles;
	hdr.tileArray = tileArray;
	hdr.nreads = N_READS;
	hdr.region_size = s->region_size;
	hdr.nregions = s->nregions;
	hdr.nregions_x = s->nregions_x;
	hdr.nregions_y = s->nregions_y;
	for (read=0; read < hdr.nreads; read++)
		hdr.readLength[read] = s->read_length[read];
	hdr.cmdLine = strdup(s->cmdline);
	hdr.ncomments = 0;
	writeHeader(fp,&hdr);
#endif
        
	HashIter *iter = HashTableIterCreate();
	while ((hashItem = HashTableIterNext(s->tile_hash, iter))) {
#ifndef BINARY_FILTER_FILE
		int tile = *(int*)(hashItem->key);
#endif
		int itile = hashItem->data.i;
		for (read = 0; read < N_READS; read++) {
			for (cycle = 0; cycle < s->read_length[read]; cycle++) {
				for (iregion = 0; iregion < s->nregions; iregion++) {
					RegionTable *rt = getRTS(itile,iregion,read,cycle);
                                        int state = rt->state;
#ifdef BINARY_FILTER_FILE
					fputc(state, fp);
#else
					if (state) fprintf(fp, "%d:%d:%d:%d:%d:%02x\n", s->lane, tile, iregion, read, cycle, state);
#endif
				}
			}
		}
	}
}

/*
 * cts - parse bam file line
 *
 * returns 0 on success, 1 on expected failure.
 */
int dump_bam_file(Settings *s, samfile_t *fp_bam, size_t *nreads)
{

	size_t nreads_bam = 0;

	static const int bam_read_buff_size = 1024;
	char bam_read_seq[bam_read_buff_size];
	int bam_read_qual[bam_read_buff_size];
	int bam_read_mismatch[bam_read_buff_size];

	bam1_t *bam = bam_init1();

	/* loop over reads in the bam file */
	while (1) {
		int bam_lane = -1, bam_tile = -1, bam_read = -1, bam_x = -1, bam_y = -1;

		if (parse_bam_readinfo(fp_bam, bam, &bam_lane, &bam_tile, &bam_x, &bam_y, &bam_read, NULL)) {
			/* break on end of BAM file */
			break;
		}

		if (BAM_FUNMAP & bam->core.flag) continue;

		parse_bam_alignments(fp_bam, bam, bam_read_seq, bam_read_qual, NULL, bam_read_mismatch,
                                                  bam_read_buff_size, s->snp_hash);

		char *name = bam1_qname(bam);
		uint8_t *oq_ptr;
		char *oq = NULL;
		int i;

		oq_ptr = bam_aux_get(bam, "OQ");
		if (NULL != oq_ptr) {
			oq = bam_aux2Z(oq_ptr);
			if (NULL == oq) {
				die("ERROR: Invalid original qualities %s for read: \"%s\"\n", oq, name);
			}
		}

		printf("%s\t%d\t%s", name, bam->core.flag, bam_read_seq);
		/* stringify quality values and mismatch bitmap into read_seq */
		for (i = 0; i < bam->core.l_qseq; i++)
			bam_read_seq[i] = bam_read_qual[i] + PHRED_QUAL_OFFSET;
		printf("\t%s", bam_read_seq);
		for (i = 0; i < bam->core.l_qseq; i++)
			bam_read_seq[i] =
			    bam_read_mismatch[i] + PHRED_QUAL_OFFSET;
		printf("\t%s", bam_read_seq);
		if (NULL != oq) {
			if (BAM_FREVERSE & bam->core.flag)
				oq = reverse_seq(oq);
			printf("\t%s", oq);
		}
		printf("\n");

		nreads_bam++;
	}

	bam_destroy1(bam);

	*nreads = nreads_bam;

	return 0;
}

/* copied from sam.c */

static bam_header_t *bam_header_dup(const bam_header_t * h0)
{
	bam_header_t *h;
	int i;
	h = bam_header_init();
	*h = *h0;
	h->hash = h->dict = h->rg2lib = 0;
	h->text = (char *)calloc(h->l_text + 1, 1);
	memcpy(h->text, h0->text, h->l_text);
	h->target_len = (uint32_t *) calloc(h->n_targets, 4);
	h->target_name = (char **)calloc(h->n_targets, sizeof(void *));
	for (i = 0; i < h->n_targets; ++i) {
		h->target_len[i] = h0->target_len[i];
		h->target_name[i] = strdup(h0->target_name[i]);
	}
	return h;
}


static int updateRegionTable(Settings *s, int itile, int read, int x, int y, int *read_mismatch)
{

        int iregion = xy2region(x, y, s->region_size, s->nregions_x, s->nregions_y);
	int cycle;

        /* update region table */
	for (cycle = 0; cycle < s->read_length[read]; cycle++) {
                RegionTable *rt = getRTS(itile,iregion,read,cycle);

		if (read_mismatch[cycle] & BASE_INSERTION) rt->insertion++;
		if (read_mismatch[cycle] & BASE_DELETION) rt->deletion++;
		if (read_mismatch[cycle] & BASE_SOFT_CLIP) rt->soft_clip++;
		if (read_mismatch[cycle] & BASE_KNOWN_SNP) { 
			rt->known_snp++;
		} else {
			if (read_mismatch[cycle] & BASE_ALIGN) rt->align++;
			if (read_mismatch[cycle] & BASE_MISMATCH) rt->mismatch++;
		}
	}

	return 0;
}

/*
 * Takes the bam file as input and updates the region table
 *
 * Assumption: within a single input file, all reads are the same length and
 * we're using unclipped data.
 *
 * Returns: 0 written for success
 *	   -1 for failure
 */
void makeRegionTable(Settings *s, samfile_t *fp_bam, int *ntiles, size_t * nreads)
{
	HashData hd;
	int lane = -1;
	int tile = -1;

	int ntiles_bam = 0;
	size_t nreads_bam = 0;

	static const int bam_read_buff_size = 1024;
	char bam_read_seq[bam_read_buff_size];
	int bam_read_qual[bam_read_buff_size];
	int bam_read_mismatch[bam_read_buff_size];

	s->tile_hash = HashTableCreate(N_TILES, HASH_DYNAMIC_SIZE | HASH_FUNC_JENKINS3);
	if (!s->tile_hash) die("Failed to create tile_hash\n");

	bam1_t *bam = bam_init1();

	int itile = -1;

	/* loop over reads in the bam file */
	while (1) {
		int bam_lane = -1, bam_tile = -1, bam_read = -1, bam_x = -1, bam_y = -1, read_length;

		if (parse_bam_readinfo(fp_bam, bam, &bam_lane, &bam_tile, &bam_x, &bam_y, &bam_read, NULL)) {
			break;	/* break on end of BAM file */
		}

		if (BAM_FUNMAP & bam->core.flag) continue;
		if (BAM_FQCFAIL & bam->core.flag) continue;
		if (BAM_FPAIRED & bam->core.flag) {
			if (0 == (BAM_FPROPER_PAIR & bam->core.flag)) {
				continue;
			}
		}

		parse_bam_alignments(fp_bam, bam, bam_read_seq, bam_read_qual, NULL, bam_read_mismatch,
                                                  bam_read_buff_size, s->snp_hash);

		read_length = strlen(bam_read_seq);
		if (0 == s->read_length[bam_read]) {
			s->read_length[bam_read] = read_length;
   	                initRTS(N_TILES, s->nregions, bam_read, read_length);
		}

		if (s->read_length[bam_read] != read_length) {
			die("Error: inconsistent read lengths within bam file for read %d.\nHave length %ld, previously it was %d.\n",
			    bam_read, (long)read_length,
			    s->read_length[bam_read]);
		}

		if (lane == -1) {
			lane = bam_lane;
		}
		if (bam_lane != lane) {
			die("Error: Inconsistent lane within bam file.\nHave %d, previously it was %d.\n", bam_lane, lane);
		}

		// lookup itile from tile in tile hash
		HashItem *tileItem = HashTableSearch(s->tile_hash, (char *)&bam_tile, sizeof(bam_tile));
		if (tileItem) {
			// if tile found, extract itile
			itile = tileItem->data.i;
		} else {
			// if tile not found in hash, add it
			hd.i = ntiles_bam++;
			if (HashTableAdd(s->tile_hash, (char *)&bam_tile, sizeof(bam_tile), hd, NULL) == NULL)
				die("Failed to add tile %d to tile_hash\n", bam_tile);
			if (!s->quiet) display("Processing tile %i (%lu)\n", bam_tile, nreads_bam);
			itile = hd.i;
		}

                if (0 != updateRegionTable(s, itile, bam_read, bam_x, bam_y, bam_read_mismatch)) {
			die("ERROR: updating quality values for tile %i.\n", tile);
		}
		nreads_bam++;
	}

	bam_destroy1(bam);

	*ntiles = ntiles_bam;
	*nreads = nreads_bam;
	s->lane = lane;
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
int filter_bam(Settings * s, samfile_t * fp_in_bam, samfile_t * fp_out_bam,
	       size_t * nreads, size_t * nfiltered)
{
	int lane = -1;
	size_t nreads_bam = 0;
	int nfiltered_bam = 0;


	bam1_t *bam = bam_init1();

	/* loop over reads in the input bam file */
	while (1) {
		int bam_lane = -1, bam_tile = -1, bam_read = -1, read_length, bam_x=-1, bam_y=-1;

		if (parse_bam_readinfo(fp_in_bam, bam, &bam_lane, &bam_tile, &bam_x, &bam_y, &bam_read, NULL)) {
			break;	/* break on end of BAM file */
		}

		read_length = bam->core.l_qseq;

		if (s->read_length[bam_read] != read_length) {
			die("Error: inconsistent read lengths within bam file for read %d.\n"
			    "have length %ld, previously it was %d.\n",
			    bam_read, (long)read_length, s->read_length[bam_read]);
		}

		if (lane == -1) {
			lane = bam_lane;
		}
		if (bam_lane != lane) {
			die("Error: Inconsistent lane: Bam lane %d qseq lane %d.\n", bam_lane, lane);
		}

		nreads_bam++;

		char state = 0, bad_cycle_count = 0;
		int iregion = xy2region(bam_x, bam_y, s->region_size, s->nregions_x, s->nregions_y);
		int read;
		for (read = 0; read < N_READS; read++) {
			int cycle;
			for (cycle = 0; cycle < s->read_length[read]; cycle++) {
				state = getFilterData(bam_tile, read, cycle, iregion);
				if (state & REGION_STATE_MASK)
                                        bad_cycle_count++;
  		        }
		}

                if (bad_cycle_count) {
                        nfiltered_bam++;
 		        if (s->qcfail) 
                                bam->core.flag |= BAM_FQCFAIL;
                        else
                                continue;
                }

		if (0 > samwrite(fp_out_bam, bam)) die("Error: writing bam file\n");
	}

	*nreads = nreads_bam;
	*nfiltered = nfiltered_bam;

	bam_destroy1(bam);

	return 0;
}

/*
 * Get the absolute file path and (depending on the libc) the real
 * name for sym-link dirs in path.
 * Safe for in_path == out_path case.
 */
static char *get_real_path_name(const char *in_path)
{
	char *oldwd;
	char *out_path;

	oldwd = getcwd(NULL,0);
	if (NULL == oldwd)
		return NULL;

	checked_chdir(in_path);

	out_path = getcwd(NULL,0);

	checked_chdir(oldwd);
	free(oldwd);

	return out_path;
}

static void usage(int code)
{
	FILE *usagefp = stderr;

	fprintf(usagefp, "spatial_filter %s\n\n", version);
	fprintf(usagefp,
		"Usage: spatial_filter [command] [options] bam_file\n"
		"" "  calculate or apply a spatial filter\n" "");
	fprintf(usagefp, "  command must be one of:\n");
	fprintf(usagefp, "    -d         just dump bam file in 'mismatch' format\n");
	fprintf(usagefp, "    -D         just dump filter file in ascii text format to stdout\n");
	fprintf(usagefp, "    -c         calculate filter files\n");
	fprintf(usagefp, "    -a         apply filter files\n");
	fprintf(usagefp, "    -v         display version and exit\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "  options:\n");
	fprintf(usagefp, "    -q         Quiet\n");
	fprintf(usagefp, "               Inhibit all progress messages\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "  comand specific options:\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "    dump bam file:\n");
	fprintf(usagefp, "      none:\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "    all other commands require:\n");
	fprintf(usagefp, "      -F --filter file\n");
	fprintf(usagefp, "                  Filter filename e.g. 8088.filter\n");
	fprintf(usagefp, "                  no default: must be supplied\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "    calculate filter:\n");
	fprintf(usagefp, "      -width width\n");
	fprintf(usagefp, "                 Tile width\n");
	fprintf(usagefp, "      -height height\n");
	fprintf(usagefp, "                 Tile height\n");
        fprintf(usagefp, "      -i --intensity-dir dir\n");
        fprintf(usagefp, "                 (Raw) Intensity directory (to find ../ImageSize.dat)\n");
        fprintf(usagefp, "                 required unless tile width and height are supplied\n");
	fprintf(usagefp, "      -s --snp_file file\n");
	fprintf(usagefp, "                 set of snps to be removed\n");
	fprintf(usagefp, "                 file in Reference Ordered Data (ROD) format\n");
	fprintf(usagefp, "      -t --tileviz read\n");
	fprintf(usagefp, "                 generate tileviz like summary plots for read number\n");
	fprintf(usagefp, "      --region_size\n");
	fprintf(usagefp, "                 default %d\n", REGION_SIZE);
	fprintf(usagefp, "      --region_min_count\n");
	fprintf(usagefp, "                 minimum coverage when setting region state\n");
	fprintf(usagefp, "                 default %d\n", REGION_MIN_COUNT);
	fprintf(usagefp, "      --region_mismatch_threshold\n");
	fprintf(usagefp, "                 threshold for setting region mismatch state\n");
	fprintf(usagefp, "                 default %-6.4f\n", REGION_MISMATCH_THRESHOLD);
	fprintf(usagefp, "      --region_insertion_threshold\n");
	fprintf(usagefp, "                 threshold for setting region insertion state\n");
	fprintf(usagefp, "                 default %-6.4f\n", REGION_INSERTION_THRESHOLD);
	fprintf(usagefp, "      --region_deletion_threshold\n");
	fprintf(usagefp, "                 threshold for setting region deletion state\n");
	fprintf(usagefp, "                 default %-6.4f\n", REGION_DELETION_THRESHOLD);
	fprintf(usagefp, "\n");
	fprintf(usagefp, "    apply filter:\n");
	fprintf(usagefp, "      -o         output\n");
	fprintf(usagefp, "                 Output bam file name\n");
	fprintf(usagefp, "                 default: stdout\n");
	fprintf(usagefp, "      -u         do not compress the output bam file\n");
	fprintf(usagefp, "                 default: compress\n");
	fprintf(usagefp, "      -f         mark filtered reads as QCFAIL\n");
	fprintf(usagefp, "                 default: do not output filtered read\n");
	fprintf(usagefp, "\n");

	exit(code);
}

char *get_command_line(int argc, char **argv)
{
	char *cmdline = NULL;
	size_t sz = argc;	/* All the spaces & the terminating \0 */
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

void calculateFilter(Settings *s)
{
	samfile_t *fp_input_bam;
	int ntiles = 0;
	size_t nreads = 0;
	FILE *fp_filter;

	fp_input_bam = samopen(s->in_bam_file, "rb", 0);
	if (NULL == fp_input_bam) {
		die("ERROR: can't open bam file file %s: %s\n", s->in_bam_file, strerror(errno));
	}

        if( NULL == s->filter) {
                display("Writing filter to stdout\n");
                s->filter = "/dev/stdout";
        }

        fp_filter = fopen(s->filter, "w+");
	if (!fp_filter) die("Can't open filter file %s: %s\n", s->filter, strerror(errno));

	makeRegionTable(s, fp_input_bam, &ntiles, &nreads);

	/* close the bam file */
	samclose(fp_input_bam);

	if (!s->quiet) {
		display("Processed %8lu traces\n", nreads);
		if (NULL != s->snp_hash) {
			size_t nsnps = 0;
			int ibucket;
			for (ibucket = 0; ibucket < s->snp_hash->nbuckets; ibucket++) {
				HashItem *hi;
				for (hi = s->snp_hash->bucket[ibucket]; hi; hi = hi->next)
					if (hi->data.i) nsnps += hi->data.i;
			}
			display("Ignored %lu snps\n", nsnps);
		}
	}

	/* back to where we belong */
	checked_chdir(s->working_dir);

	findBadRegions(s, ntiles);

	printFilter(s, ntiles, fp_filter);

	freeRTS();
	fclose(fp_filter);
}

void applyFilter(Settings *s)
{
	samfile_t *fp_input_bam;
	samfile_t *fp_output_bam;
	char out_mode[5] = "wb";
	char *out_bam_file = NULL;
	bam_header_t *out_bam_header = NULL;
	size_t nreads = 0;
	size_t nfiltered = 0;
	Header hdr;
	FILE *fp;
        int read;

	fp = fopen(s->filter, "rb");
	if (!fp) die("Can't open filter file %s\n", s->filter);
	readHeader(fp, &hdr);
	readFilterData(fp, &hdr);
	fclose(fp);

	s->region_size = hdr.region_size;
	s->nregions_x = hdr.nregions_x;
	s->nregions_y = hdr.nregions_y;

        for (read=0;read<hdr.nreads;read++)
                s->read_length[read] = hdr.readLength[read];

	if (0 == s->compress)
		strcat(out_mode, "u");

	fp_input_bam = samopen(s->in_bam_file, "rb", 0);
	if (NULL == fp_input_bam) {
		die("ERROR: can't open bam file file %s: %s\n", s->in_bam_file, strerror(errno));
	}

	out_bam_header = bam_header_dup(fp_input_bam->header);
	char concat_cmd[2048];
	strcat(concat_cmd, hdr.cmdLine);
	strcat(concat_cmd, " ; ");
	strcat(concat_cmd, s->cmdline);
	bam_header_add_pg("spf", "spatial_filter", "A program to apply a spatial filter", concat_cmd, out_bam_header);

	out_bam_file = (NULL == s->output ? aprintf("/dev/stdout") : aprintf("%s/%s", s->working_dir, s->output));
	fp_output_bam = samopen(out_bam_file, out_mode, out_bam_header);
	if (NULL == fp_output_bam) {
		die("ERROR: can't open bam file file %s: %s\n", out_bam_file, strerror(errno));
	}
	free(out_bam_file);

	bam_header_destroy(out_bam_header);

	if (-1 == filter_bam(s, fp_input_bam, fp_output_bam, &nreads, &nfiltered)) {
		die("ERROR: failed to filter bam file %s\n", s->in_bam_file);
	}

	samclose(fp_input_bam);
	samclose(fp_output_bam);

	if (!s->quiet) display("Processed %8lu traces\n", nreads);
	if (!s->quiet) display("%s %8lu traces\n", (s->qcfail ? "QC failed" : "Removed"), nfiltered);

	/* back to where we belong */
	checked_chdir(s->working_dir);
}

void dumpBAM(Settings *s)
{
	samfile_t *fp_input_bam;
		size_t nreads = 0;

		fp_input_bam = samopen(s->in_bam_file, "rb", 0);
		if (NULL == fp_input_bam) {
			die("ERROR: can't open bam file file %s: %s\n", s->in_bam_file, strerror(errno));
		}

		if (0 != dump_bam_file(s, fp_input_bam, &nreads)) {
			die("ERROR: failed to dump bam file %s\n", s->in_bam_file);
		}

		/* close the bam file */
		samclose(fp_input_bam);

		if (!s->quiet) {
			display("Dumped %8lu traces\n", nreads);
		}
}

void dumpFilterFile(char *filename)
{
	FILE *fp;
	Header hdr;
	int i;

	if (!filename) die("dumpFilterFile: no filter filename given\n");
	fp = fopen(filename, "r");
	if (!fp) die("dumpFilterFile: Can't open file %s\n", filename);
	readHeader(fp,&hdr);
	printf("Magic:          %s\n", hdr.region_magic);
	printf("Coord Shift:    %-5d\n", hdr.coord_shift);
	printf("Coord Factor:   %-5d\n", hdr.coord_shift);
	printf("Region Size:    %-5d\n", hdr.region_size);
	printf("Num Regions:    %-5d\n", hdr.nregions);
	printf("Num Tiles:      %-5d\n", hdr.ntiles);
	for (i=0; i < hdr.ntiles; i++) printf("%-5d ", hdr.tileArray[i]);
	printf("\n");
	printf("Read Length:    ");
	for (i=0; i < hdr.nreads; i++) printf("%-5d ", hdr.readLength[i]);
	printf("\n");
	printf("Command Line:   %s\n", hdr.cmdLine);
	for (i=0; i < hdr.ncomments; i++) {
		if (i) printf("                %s\n", hdr.comments[i]);
		else   printf("Comments:       %s\n", hdr.comments[i]);
	}
	fclose(fp);
}


int main(int argc, char **argv)
{
	Settings settings;
	const char *override_intensity_dir = NULL;
	int dumpFilter = 0;
	char *cmd = NULL;

	settings.filter = NULL;
	settings.tileViz = 0;
	settings.quiet = 0;
	settings.dump = 0;
	settings.calculate = 0;
	settings.apply = 0;
	settings.qcfail = 0;
	settings.intensity_dir = NULL;
	settings.snp_file = NULL;
	settings.snp_hash = NULL;
	settings.tile_hash = NULL;
	settings.output = NULL;
	settings.in_bam_file = NULL;
	settings.read_length[0] = 0;
	settings.read_length[1] = 0;
	settings.read_length[2] = 0;
	settings.working_dir = NULL;
	settings.region_size = REGION_SIZE;
	settings.width = 0;
	settings.height = 0;
	settings.nregions_x = 0;
	settings.nregions_y = 0;
	settings.nregions = 0;
	settings.compress = 1;
	settings.region_min_count = REGION_MIN_COUNT;
	settings.region_mismatch_threshold = REGION_MISMATCH_THRESHOLD;
	settings.region_insertion_threshold = REGION_INSERTION_THRESHOLD;
	settings.region_deletion_threshold = REGION_DELETION_THRESHOLD;

	static struct option long_options[] = {
                   {"intensity_dir", 1, 0, 'i'},
                   {"intensity-dir", 1, 0, 'i'},
                   {"snp_file", 1, 0, 's'},
                   {"snp-file", 1, 0, 's'},
                   {"tileviz", 1, 0, 't'},
                   {"help", 0, 0, 'h'},
                   {"filter", 1, 0, 'F'},
                   {"width", 1, 0, 'x'},
                   {"height", 1, 0, 'y'},
                   {"version", 0, 0, 'v'},
                   {"region_min_count", 1, 0, 'm'},
                   {"region-size", 1, 0, 'r'},
                   {"region_size", 1, 0, 'r'},
                   {"region_mismatch_threshold", 1, 0, 'z'},
                   {"region_insertion_threshold", 1, 0, 'b'},
                   {"region_deletion_threshold", 1, 0, 'e'},
                   {0, 0, 0, 0}
               };

	int ncmd = 0;
        char c;
	while ( (c = getopt_long(argc, argv, "vdcafuDF:b:e:o:i:m:p:s:r:x:y:t:z:qh?", long_options, 0)) != -1) {
		switch (c) {
			case 'v':	display("spatial_filter: Version %s\n", version); 
						exit(0);
                        case 'd':	settings.dump = 1; ncmd++; break;
			case 'D':	dumpFilter = 1; ncmd++; break;
			case 'c':	settings.calculate = 1; ncmd++; break;
			case 'a':	settings.apply = 1; ncmd++; break;
			case 't':	parse_next_int(optarg,&settings.tileViz,NULL);	break;
			case 'f':	settings.qcfail = 1;		break;
			case 'u':	settings.compress = 0;		break;
			case 'o':	settings.output = optarg;	break;
			case 'i':	override_intensity_dir = optarg;	break;
			case 's':	settings.snp_file = optarg;	break;
			case 'F':	settings.filter = optarg;	break;
			case 'x':	settings.width = atoi(optarg);	break;
			case 'y':	settings.height = atoi(optarg);	break;
			case 'r':	settings.region_size = atoi(optarg); break;
			case 'm':	settings.region_min_count = atoi(optarg); break;
			case 'z':	settings.region_mismatch_threshold = atof(optarg); break;
			case 'b':	settings.region_insertion_threshold = atof(optarg); break;
			case 'e':	settings.region_deletion_threshold = atof(optarg); break;
			case 'q':	settings.quiet = 1;			break;
			case 'h':
			case '?':	usage(0);					break;
			default:	display("ERROR: Unknown option %c\n", c);
						usage(1);
						break;
		}
	}

	if (ncmd > 1) {
	        display("ERROR: More than one command specified\n", c);
		usage(1);
        }

	if (optind < argc) settings.in_bam_file = argv[optind];

	if (!settings.in_bam_file && !dumpFilter) die("Error: no BAM file specified\n");

        if (!settings.filter && (dumpFilter || settings.apply)) die("Error: no filter file specified\n");

	if (settings.calculate) {
            if (settings.width < 0) die("Error: invalid tile width");
	    if (settings.height < 0) die("Error: invalid tile height");

            if (!override_intensity_dir && (!settings.width || !settings.height))
                die("ERROR: you must specify an intensity dir or tile width and height\n");

   	    if (settings.region_size < 1) die("Error: invalid region size");

	    if (settings.region_min_count < 1) die("Error: invalid region_min_count");

#if 0 // this code will reset region_min_count so that at least 2 reads are required to pass any threshold
            int region_min_count = settings.region_min_count;

	    if ((region_min_count * settings.region_mismatch_threshold) < 1.0 ) {
                region_min_count = ceil(1.0 / settings.region_mismatch_threshold);
                display("Resetting region_min_count to %d\n", region_min_count);
            }
	    if ((region_min_count * settings.region_insertion_threshold) < 1.0 ) {
                region_min_count = ceil(1.0 / settings.region_insertion_threshold);
                display("Resetting region_min_count to %d\n", region_min_count);
            }
	    if ((region_min_count * settings.region_deletion_threshold) < 1.0 ) {
                region_min_count = ceil(1.0 / settings.region_deletion_threshold);
                display("Resetting region_min_count to %d\n", region_min_count);
            }
            settings.region_min_count = region_min_count;
#endif            
        }

	// create pseudo command line
	if (settings.calculate) {
		char arg[64];
		cmd = smalloc(2048);
		strcat(cmd, argv[0]);
		strcat(cmd, " -c ");
		if (override_intensity_dir)              { strcat(cmd, " -i "); strcat(cmd, override_intensity_dir); }
		if (settings.snp_file)                   { strcat(cmd, " -s "); strcat(cmd, settings.snp_file); }
		if (settings.filter)                     { strcat(cmd, " -F "); strcat(cmd, settings.filter); }
		if (settings.width)                      { snprintf(arg, 64, " --width %d", settings.width);
                                                           strcat(cmd, arg); }
		if (settings.height)                     { snprintf(arg, 64, " --height %d", settings.height);
                                                           strcat(cmd, arg); }
		if (settings.region_size)                { snprintf(arg, 64, " --region_size %d", settings.region_size);
                                                           strcat(cmd, arg); }
		if (settings.region_min_count)           { snprintf(arg, 64, " --region_min_count %d", settings.region_min_count);
                                                           strcat(cmd, arg); }
		if (settings.region_mismatch_threshold)  { snprintf(arg, 64, " --region_mismatch_threshold %-6.4f", settings.region_mismatch_threshold);
                                                           strcat(cmd, arg); }
		if (settings.region_insertion_threshold) { snprintf(arg, 64, " --region_insertion_threshold %-6.4f", settings.region_insertion_threshold);
                                                           strcat(cmd, arg); }
		if (settings.region_deletion_threshold)  { snprintf(arg, 64, " --region_deletion_threshold %-6.4f", settings.region_deletion_threshold);
                                                           strcat(cmd, arg); }
		strcat(cmd, " ");
		strcat(cmd, settings.in_bam_file);
		if (strlen(cmd) > 2047) die("Command line too big");
	}
	if (settings.apply) {
		cmd = smalloc(2048);
		strcat(cmd, argv[0]);
		strcat(cmd, " -a ");
		if (settings.qcfail)    { strcat(cmd, " -f "); }
		if (!settings.compress) { strcat(cmd, " -u "); }
		if (settings.output)    { strcat(cmd, " -o "); strcat(cmd, settings.output); }
		if (settings.filter)    { strcat(cmd, " -F "); strcat(cmd, settings.filter); }
		strcat(cmd, " ");
		if (settings.in_bam_file) strcat(cmd, settings.in_bam_file);
		if (strlen(cmd) > 2047) die("Command line too big");
	}

	if (!settings.quiet) display("Cmd: %s\n", cmd);
	settings.cmdline = cmd;

        if (override_intensity_dir) {
                settings.intensity_dir = get_real_path_name(override_intensity_dir);
                if (NULL == settings.intensity_dir) {
                        die("ERROR: can't process intensity dir: %s\n", override_intensity_dir);
                }
        }

	/* preserve starting directory */
	settings.working_dir = getcwd(NULL,0);
	if (NULL == settings.working_dir) {
		die("ERROR: can't obtain working directory: %s\n", strerror(errno));
	}

	/* dump the alignments */
	if (settings.dump) dumpBAM(&settings);

	/* Dump the filter file */
	if (dumpFilter) dumpFilterFile(settings.filter);

	/* calculate the filter */
        if (settings.calculate) {
                /* set the number of regions */
                setRegions(&settings);

                /* read the snp_file */
                if (NULL != settings.snp_file) {
                        settings.snp_hash = readSnpFile(settings.snp_file);
                        if (NULL == settings.snp_hash) {
                                die("ERROR: reading snp file %s\n", settings.snp_file);
                        }
                }

                calculateFilter(&settings);
        }

	/* apply the  filter */
	if (settings.apply) applyFilter(&settings);

	if (NULL != settings.working_dir) free(settings.working_dir);

	return EXIT_SUCCESS;

}
