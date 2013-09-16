//*  -*- mode: c; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 8; -*- */

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
 * TILEVIZ
 *   make tileviz PNG's when calculating filter
 *   to be consistent with perl tileviz need to include reads that fail QC or are not properly paired
 */

#define QC_FAIL
#define PROPERLY_PAIRED
//#define TILEVIZ

#ifdef HAVE_CONFIG_H
#include "pb_config.h"
#endif

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <search.h>
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
#include <gd.h>
#include <gdfonts.h>

/* To turn off assert for a small speed gain, uncomment this line */
/* #define NDEBUG */

/* Hack to stop io_lib from trying to include its own config.h */
#ifdef HAVE_CONFIG_H
#undef HAVE_CONFIG_H
#endif
#include <io_lib/misc.h>
#include <io_lib/pooled_alloc.h>
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

#define REGION_MIN_COUNT            0      // minimum coverage when setting region state
#define REGION_MISMATCH_THRESHOLD   0.016  // threshold for setting region mismatch state
#define REGION_INSERTION_THRESHOLD  0.016  // threshold for setting region insertion state
#define REGION_DELETION_THRESHOLD   0.016  // threshold for setting region deletion state

#define TILE_REGION_THRESHOLD  0.75  // threshold for setting region state at tile level

#define TILE_WIDTH   2048   // default tile width,  used to set initial size of region hash table
#define TILE_HEIGHT  10000  // default tile height, used to set initial size of region hash table

#define MAX_NREGIONS_X         100  // max IX is MAX_NREGIONS_X-1
#define MAX_NREGIONS_Y         1000 // max IY is MAX_NREGIONS_Y-1
#define MAX_REGION_KEY_LENGTH  7    // key is IX:IY so maximum length key is 99:999
    
#define REGION_STATE_MASK  (REGION_STATE_INSERTION | REGION_STATE_DELETION)  // region mask used to filter reads

enum images { IMAGE_COVERAGE,
              IMAGE_DELETION,
              IMAGE_INSERTION,
              IMAGE_MISMATCH,
              IMAGE_QUALITY1,
              IMAGE_QUALITY2,
              N_IMAGES };

char *image_names[] = {"cov",
                       "del",
                       "ins",
                       "mma",
                       "qua",
                       "oqu"};

#define IMAGE_COLUMN_GAP 3
#define IMAGE_LABEL_HEIGHT 25
    
enum colours { COLOUR_LEVEL_0,
               COLOUR_LEVEL_1,
               COLOUR_LEVEL_2,
               COLOUR_LEVEL_3,
               COLOUR_LEVEL_4,
               COLOUR_LEVEL_5,
               COLOUR_LEVEL_6,
               COLOUR_LEVEL_7,
               COLOUR_LEVEL_8,
               COLOUR_LEVEL_9,
               COLOUR_LEVEL_10,
               COLOUR_LEVEL_11,
               COLOUR_TEXT,
               COLOUR_QC_FAIL,
               COLOUR_ZERO_QUAL,
               COLOUR_LOW_QUAL,
               COLOUR_MEDIUM_QUAL,
               COLOUR_HIGH_QUAL,
               N_COLOURS };

int *colour_table = NULL;

typedef struct {
	char *cmdline;
	char *filter;
	char *snp_file;
	char *in_bam_file;
	HashTable *snp_hash;
	int *tileArray;
	char *working_dir;
	char *output;
	int read_length[3];
	int dump;
	int calculate;
	char *tileviz;
	int apply;
	int qcfail;
	int quiet;
	HashTable *region_hash;
	int region_size;
	int nregions_x;
	int nregions_y;
	int nregions;
	int *regions;
	int compress;
	int region_min_count;
	float region_mismatch_threshold;
	float region_insertion_threshold;
	float region_deletion_threshold;
} Settings;

void freeRTS(Settings *s, int ntiles, RegionTable ***rts)
{
        int itile, read, cycle;
        for(itile=0;itile<ntiles;itile++)
            for(read=0;read<N_READS;read++)
            {
                if( NULL == rts[itile*N_READS+read]) continue;
                for(cycle=0;cycle<s->read_length[read];cycle++)
                    free(rts[itile*N_READS+read][cycle]);
                free(rts[itile*N_READS+read]);
            }
        free(rts);
}

/*
 * initialise tileviz image
*/

static gdImagePtr initImage(int width, int height, char *base, int type, int read, int cycle)
{
    gdImagePtr im = gdImageCreate(width, height);

    if( NULL == colour_table ){
        colour_table = smalloc(N_COLOURS * sizeof(int));
    }

    gdImageColorAllocate(im, 0, 0, 0); // black - the background colour

    // white + graduated shades of blue from light to dark
    colour_table[COLOUR_LEVEL_0]  = gdImageColorAllocate(im, 255, 255, 255);
    colour_table[COLOUR_LEVEL_1]  = gdImageColorAllocate(im, 211, 222, 235);   
    colour_table[COLOUR_LEVEL_2]  = gdImageColorAllocate(im, 189, 206, 225);
    colour_table[COLOUR_LEVEL_3]  = gdImageColorAllocate(im, 167, 190, 215);
    colour_table[COLOUR_LEVEL_4]  = gdImageColorAllocate(im, 145, 174, 205);
    colour_table[COLOUR_LEVEL_5]  = gdImageColorAllocate(im, 124, 157, 195);
    colour_table[COLOUR_LEVEL_6]  = gdImageColorAllocate(im, 102, 141, 185);
    colour_table[COLOUR_LEVEL_7]  = gdImageColorAllocate(im,  80, 125, 175);
    colour_table[COLOUR_LEVEL_8]  = gdImageColorAllocate(im,  58, 109, 165);
    colour_table[COLOUR_LEVEL_9]  = gdImageColorAllocate(im,  36,  93, 155);
    colour_table[COLOUR_LEVEL_10] = gdImageColorAllocate(im,  15,  77, 146);
    colour_table[COLOUR_LEVEL_11] = gdImageColorAllocate(im,   0,  61, 136);

    // specific colours
    colour_table[COLOUR_TEXT]        = gdImageColorAllocate(im, 239, 239, 239); // light grey
    colour_table[COLOUR_QC_FAIL]     = gdImageColorAllocate(im, 255,   0,   0); // red
    colour_table[COLOUR_ZERO_QUAL]   = gdImageColorAllocate(im, 255,   0,   0); // red
    colour_table[COLOUR_LOW_QUAL]    = gdImageColorAllocate(im, 244, 211,  71); // yellow
    colour_table[COLOUR_MEDIUM_QUAL] = gdImageColorAllocate(im,  21,  58, 144); // dark blue
    colour_table[COLOUR_HIGH_QUAL]   = gdImageColorAllocate(im, 185, 212, 246); // light blue
    
    char str[1024];

    if( NULL != base ) gdImageString(im, gdFontSmall, 3, 1, (unsigned char *)base, colour_table[COLOUR_TEXT]);

    if( cycle < 0 ){
        sprintf(str, "%c_%s", (read == 1 ? 'F' : 'R'), image_names[type]);
    }else{
        sprintf(str, "%03d%c_%s", cycle, (read == 1 ? 'F' : 'R'), image_names[type]);
    }
    gdImageString(im, gdFontSmall, 3, 11, (unsigned char *)str, colour_table[COLOUR_TEXT]);

    return im;
}

/*
 * generate tileviz plots
*/

static void tileviz(Settings *s, int ntiles, RegionTable ***rts)
{
    int num_surfs = 2;
    int num_cols = 2;
    int num_rows = 16;
    int image_width = s->nregions_x * num_cols * num_surfs + IMAGE_COLUMN_GAP;
    int image_height = (s->nregions_y + 1) * num_rows + IMAGE_LABEL_HEIGHT;
    gdImagePtr im[N_IMAGES];
#if 0
    int oqu = 0;
#else
    int oqu = 1;
#endif
    char *base;
    int filename_sz;
    char *filename;
    FILE *fp;
    int image, ix, iy, read, itile, cycle;

	if (0 >= ntiles)
		return;

    base = strrchr(s->tileviz, '/');
    if( NULL == base )
        base = s->tileviz;
    else
        base++;

    display("Writing tileviz plots to %s\n", s->tileviz);

    filename_sz = (NULL == s->tileviz ? 0 : strlen(s->tileviz)) + 100;
    filename = smalloc(filename_sz);

    for (image=0; image<N_IMAGES; image++)
        im[image] = NULL;

    // create a summary RT
    for (read = 0; read < N_READS; read++) {
        RegionTable summary_rts[ntiles][s->nregions_x][s->nregions_y];
        for (itile=0; itile<ntiles; itile++) {
            int iregion = 0;
            for (ix = 0; ix < s->nregions_x; ix++) {
                for (iy = 0; iy < s->nregions_y; iy++) {
                    RegionTable *summary_rt = &summary_rts[itile][ix][iy];
                    summary_rt->align     = 0;
                    summary_rt->mismatch  = 0;
                    summary_rt->insertion = 0;
                    summary_rt->deletion  = 0;
                    summary_rt->quality1  = 1.0E+6;
                    summary_rt->quality2  = 1.0E+6;
                    if (s->regions[iregion] >= 0) {
                        for (cycle = 0; cycle < s->read_length[read]; cycle++) {
                            RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                            int n = rt->align + rt->insertion + rt->deletion + rt->soft_clip + rt->known_snp;
                            // coverage should be the same for all cycles
                            summary_rt->align = n;
                            if (0 == summary_rt->align) continue;
                            // for mismatch, insertion and deletion convert to a percentage and bin 0(0%), 1(10%), .. 10(100%)
                            rt->mismatch = 10.0 * rt->mismatch / n;
                            rt->insertion = 10.0 * rt->insertion / n;
                            rt->deletion = 10.0 * rt->deletion / n;
                            // for quality values calculate an average value
                            rt->quality1 /= n;
                            rt->quality2 /= n;
                            // ignore the first/last cycles which have a high error rate and lower quality values
                            if( cycle == 0 || cycle == (s->read_length[read]-1) ) continue;
                            // for mismatch, insertion and deletion take the maximum over all cycles
                            summary_rt->mismatch = max(summary_rt->mismatch, rt->mismatch);
                            summary_rt->insertion = max(summary_rt->insertion, rt->insertion);
                            summary_rt->deletion = max(summary_rt->deletion, rt->deletion);
                            // for quality values take the minimum over all cycles
                            summary_rt->quality1 = min(summary_rt->quality1, rt->quality1);
                            summary_rt->quality2 = min(summary_rt->quality2, rt->quality2);
                        }
                    }
                    iregion++;
                }
			}
		}

        // create summary plots

        im[IMAGE_COVERAGE]  = initImage(image_width, image_height, base, IMAGE_COVERAGE,  read, -1);
        im[IMAGE_DELETION]  = initImage(image_width, image_height, base, IMAGE_DELETION,  read, -1);
        im[IMAGE_INSERTION] = initImage(image_width, image_height, base, IMAGE_INSERTION, read, -1);
        im[IMAGE_MISMATCH]  = initImage(image_width, image_height, base, IMAGE_MISMATCH,  read, -1);
        im[IMAGE_QUALITY1]  = initImage(image_width, image_height, base, IMAGE_QUALITY1,  read, -1);
        if(oqu) im[IMAGE_QUALITY2]  = initImage(image_width, image_height, base, IMAGE_QUALITY2,  read, -1);

        for (itile=0; itile<ntiles; itile++) {
            int tile = s->tileArray[itile];
            int surf = tile / 1000;
            int col = (tile - 1000 * surf) / 100;
            int row = tile % 100;

            for (ix = 0; ix < s->nregions_x; ix++) {
                for (iy = 0; iy < s->nregions_y; iy++) {
                    RegionTable *rt = &summary_rts[itile][ix][iy];
                    int x = (surf-1) * (s->nregions_x * num_cols + IMAGE_COLUMN_GAP) + (col-1) * s->nregions_x + ix;
                    int y = IMAGE_LABEL_HEIGHT + (row-1) * (s->nregions_y + 1) + iy;
                    if (0 == rt->align) continue;
                    int colour = (rt->align > COLOUR_LEVEL_11 ? COLOUR_LEVEL_11 : rt->align);
                    gdImageSetPixel(im[IMAGE_COVERAGE],  x, y, colour_table[colour]);
                    colour = rt->deletion;
                    gdImageSetPixel(im[IMAGE_DELETION],  x, y, colour_table[colour]);
                    colour = rt->insertion;
                    gdImageSetPixel(im[IMAGE_INSERTION], x, y, colour_table[colour]);
                    colour = rt->mismatch;
                    gdImageSetPixel(im[IMAGE_MISMATCH],  x, y, colour_table[colour]);
                    if (rt->quality1 > 30) {
                        colour = COLOUR_HIGH_QUAL;
                    } else if (rt->quality1 > 15) {
                        colour = COLOUR_MEDIUM_QUAL;
                    } else if (rt->quality1 < 5) {
                        colour = COLOUR_ZERO_QUAL;
                    } else {
                        colour = COLOUR_LOW_QUAL;
                    }
                    gdImageSetPixel(im[IMAGE_QUALITY1],  x, y, colour_table[colour]);
                    if( NULL!= im[IMAGE_QUALITY2] ) {
                        if (rt->quality2 > 30) {
                            colour = COLOUR_HIGH_QUAL;
                        } else if (rt->quality2 > 15) {
                            colour = COLOUR_MEDIUM_QUAL;
                        } else if (rt->quality2 < 5) {
                            colour = COLOUR_ZERO_QUAL;
                        } else {
                            colour = COLOUR_LOW_QUAL;
                        }
                        gdImageSetPixel(im[IMAGE_QUALITY2],  x, y, colour_table[colour]);
                    }
                }
            }
        }
        
        for (image=0; image<N_IMAGES; image++) {
            if( NULL == im[image] ) continue;
            sprintf(filename, "%s_%c_%s.png", s->tileviz, (read == 1 ? 'F' : 'R'), image_names[image]);
            fp = fopen(filename, "w+");
            if (!fp) die("Can't open tileviz file %s: %s\n", filename, strerror(errno));
            gdImagePng(im[image], fp);
            fclose(fp);
            gdImageDestroy(im[image]);
            im[image] = NULL;
        }
        
        // create cycle by cycle plots
        for (cycle = 0; cycle < s->read_length[read]; cycle++) {

            im[IMAGE_DELETION]  = initImage(image_width, image_height, base, IMAGE_DELETION,  read, cycle);
            im[IMAGE_INSERTION] = initImage(image_width, image_height, base, IMAGE_INSERTION, read, cycle);
            im[IMAGE_MISMATCH]  = initImage(image_width, image_height, base, IMAGE_MISMATCH,  read, cycle);
            im[IMAGE_QUALITY1]  = initImage(image_width, image_height, base, IMAGE_QUALITY1,  read, cycle);
            if(oqu) im[IMAGE_QUALITY2]  = initImage(image_width, image_height, base, IMAGE_QUALITY2,  read, cycle);

            for (itile=0; itile<ntiles; itile++) {
                int tile = s->tileArray[itile];
                int surf = tile / 1000;
                int col = (tile - 1000 * surf) / 100;
                int row = tile % 100;
                int iregion = 0;

                for (ix = 0; ix < s->nregions_x; ix++) {
                    for (iy = 0; iy < s->nregions_y; iy++) {
                        if (s->regions[iregion] >= 0) {
                            RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                            int x = (surf-1) * (s->nregions_x * num_cols + IMAGE_COLUMN_GAP) + (col-1) * s->nregions_x + ix;
                            int y = IMAGE_LABEL_HEIGHT + (row-1) * (s->nregions_y + 1) + iy;
                            if (0 == rt->align) continue;
                            int colour = rt->deletion;
                            gdImageSetPixel(im[IMAGE_DELETION],  x, y, colour_table[colour]);
                            colour = rt->insertion;
                            gdImageSetPixel(im[IMAGE_INSERTION], x, y, colour_table[colour]);
                            colour = rt->mismatch;
                            gdImageSetPixel(im[IMAGE_MISMATCH],  x, y, colour_table[colour]);
                            if (rt->quality1 > 30) {
                                colour = COLOUR_HIGH_QUAL;
                            } else if (rt->quality1 > 15) {
                                colour = COLOUR_MEDIUM_QUAL;
                            } else if (rt->quality1 < 5) {
                                colour = COLOUR_ZERO_QUAL;
                            } else {
                                colour = COLOUR_LOW_QUAL;
                            }
                            gdImageSetPixel(im[IMAGE_QUALITY1],  x, y, colour_table[colour]);
                            if( NULL!= im[IMAGE_QUALITY2] ) {
                                if (rt->quality2 > 30) {
                                    colour = COLOUR_HIGH_QUAL;
                                } else if (rt->quality2 > 15) {
                                    colour = COLOUR_MEDIUM_QUAL;
                                } else if (rt->quality2 < 5) {
                                    colour = COLOUR_ZERO_QUAL;
                                } else {
                                    colour = COLOUR_LOW_QUAL;
                                }
                                gdImageSetPixel(im[IMAGE_QUALITY2],  x, y, colour_table[colour]);
                            }
                        }
                    }
                    iregion++;
                }
            }
            
            for (image=0; image<N_IMAGES; image++) {
                if( NULL == im[image] ) continue;
                sprintf(filename, "%s_%03d%c_%s.png", s->tileviz, cycle+1, (read == 1 ? 'F' : 'R'), image_names[image]);
                fp = fopen(filename, "w+");
                if (!fp) die("Can't open tileviz file %s: %s\n", filename, strerror(errno));
                gdImagePng(im[image], fp);
                fclose(fp);
                gdImageDestroy(im[image]);
                im[image] = NULL;
            }
        }
    }

	return;
}

/*
 * initialise the region table
 */

void initialiseRegionTable(RegionTable *rt) {
    rt->align     = 0;
    rt->mismatch  = 0;
    rt->insertion = 0;
    rt->deletion  = 0;
    rt->soft_clip = 0;
    rt->known_snp = 0;
    rt->quality1  = 0;
    rt->quality2  = 0;
    rt->state     = 0;
}

/*
 * setup a mapping between each potential region and the observed regions
 */

void regionMapping(Settings *s) {
    int iregion = 0, ix, iy;
    if( NULL != s->regions) free(s->regions);
    s->regions = smalloc(s->nregions * sizeof(int));
    for (ix = 0; ix < s->nregions_x; ix++) {
        for (iy = 0; iy < s->nregions_y; iy++) {
            char key[100];
            char *cp;
            HashItem *hi;
            cp = append_int(key, ix);
            cp = append_char(cp, ':');
            cp = append_int(cp, iy);
            *cp = 0;
            if( NULL != (hi = HashTableSearch(s->region_hash, key, strlen(key))) ){
                s->regions[iregion++] = hi->data.i;
            }else{
                s->regions[iregion++] = -1;
            }
        }
    }

    return;
}

/*
 * increase the region size until the average region counts exceeds the minimum region count
*/

static void resizeRegions(Settings *s, int ntiles, size_t nreads, RegionTable ***rts)
{
	int iregion, ix, iy, itile, read, cycle;

	if (0 >= ntiles)
		return;

    // if the minimum region count is not set, set it so that at least 2 reads are required to pass all thresholds
    if( 0 == s->region_min_count ){
        int region_min_count = s->region_min_count;
        if ((region_min_count * s->region_mismatch_threshold) < 2.0 )
            region_min_count = ceil(2.0 / s->region_mismatch_threshold);
        if ((region_min_count * s->region_insertion_threshold) < 2.0 )
            region_min_count = ceil(2.0 / s->region_insertion_threshold);
        if ((region_min_count * s->region_deletion_threshold) < 2.0 )
            region_min_count = ceil(2.0 / s->region_deletion_threshold);
        if (region_min_count != s->region_min_count) {
            display("Setting region_min_count to %d\n", region_min_count);
            s->region_min_count = region_min_count;
        }
    }

    int region_size = s->region_size;
    int nregions_x = s->nregions_x;
    int nregions_y = s->nregions_y;
    int nregions = s->nregions;

    // what is the average #reads per region, assume coverage is reasonably uniform over the whole lane
    int region_count = (float)nreads / (float)(ntiles * nregions);
    display("nregions_x=%d nregions_y=%d nregions=%d region_size=%d region_count=%d region_min_count=%d\n",
            nregions_x, nregions_y, nregions, region_size, region_count, s->region_min_count);

    // increase the region size until at the average region count exceeds the minimum region count
    int scale_factor = 1;
    while ( region_count < s->region_min_count ){
        scale_factor++;
        region_size = scale_factor * s->region_size;
        nregions_x = s->nregions_x / scale_factor + 1;
        nregions_y = s->nregions_y / scale_factor + 1;
        nregions = nregions_x * nregions_y;
        region_count = (float)nreads / (float)(ntiles * nregions);
        display("nregions_x=%d nregions_y=%d nregions=%d region_size=%d region_count=%d region_min_count=%d\n",
                nregions_x, nregions_y, nregions, region_size, region_count, s->region_min_count);
    }

    // nothing to do if we didn't change the region size
    if (scale_factor == 1) return;

    /* setup the new region_hash and a mapping between old and new regions */
    HashTable *region_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE|HASH_FUNC_JENKINS3);
    int *new_regions = smalloc(s->nregions * sizeof(int));
    iregion = 0;
    for (ix = 0; ix < s->nregions_x; ix++) {
        for (iy = 0; iy < s->nregions_y; iy++) {
            if (s->regions[iregion] >= 0) {
                int new_ix = ix / scale_factor;
                int new_iy = iy / scale_factor;
                char key[100];
                char *cp;
                HashItem *hi;
                cp = append_int(key, new_ix);
                cp = append_char(cp, ':');
                cp = append_int(cp, new_iy);
                *cp = 0;
                if( NULL != (hi = HashTableSearch(region_hash, key, strlen(key))) ){
                    new_regions[iregion] = hi->data.i;
                }else{
                    HashData hd;
                    hd.i = region_hash->nused;
                    if( NULL == HashTableAdd(region_hash, key, strlen(key), hd, NULL) ) {
                        fprintf(stderr, "ERROR: building new rts hash table\n");
                        exit(EXIT_FAILURE);
                    }
                    new_regions[iregion] = hd.i;
                }
            }else{
                new_regions[iregion] = -1;
            }
            iregion++;
        }
    }

    // make the new region tables
    for (itile=0; itile<ntiles; itile++) {
        for (read = 0; read < N_READS; read++) {
            for (cycle = 0; cycle < s->read_length[read]; cycle++) {
                RegionTable *new_rts = smalloc(nregions * sizeof(RegionTable));
                for (iregion=0;iregion<nregions;iregion++) {
                    RegionTable *new_rt = new_rts + iregion;
                    initialiseRegionTable(new_rt);
                }
                for (iregion=0;iregion<s->nregions;iregion++) {
                    if (s->regions[iregion] >= 0) {
                        RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                        RegionTable *new_rt = new_rts + new_regions[iregion];
                        new_rt->align     += rt->align;
                        new_rt->mismatch  += rt->mismatch;
                        new_rt->insertion += rt->insertion;
                        new_rt->deletion  += rt->deletion;
                        new_rt->soft_clip += rt->soft_clip;
                        new_rt->known_snp += rt->known_snp;
                        new_rt->quality1  += rt->quality1;
                        new_rt->quality2  += rt->quality2;
                    }
                }
                free(rts[itile*N_READS+read][cycle]);
                rts[itile*N_READS+read][cycle] = new_rts;
            }
        }
    }

    free(new_regions);

    HashTableDestroy(s->region_hash, 0);
    s->region_hash = region_hash;

    s->region_size = region_size;
    s->nregions_x = nregions_x;
    s->nregions_y = nregions_y;
    s->nregions = nregions;

    /* reset the mapping between each potential region and the observed regions */
    regionMapping(s);

	return;
}

/*
 * identify bad tiles quality < mean quality - filter * stdev quality
 * filter is typically 2.
*/

static void findBadRegions(Settings *s, int ntiles, size_t nreads, RegionTable ***rts)
{
	int itile, read, cycle, iregion;

	if (0 >= ntiles)
		return;

    for (itile=0; itile<ntiles; itile++) {
		int tile = s->tileArray[itile];
        for (read = 0; read < N_READS; read++) {
			for (cycle = 0; cycle < s->read_length[read]; cycle++) {
				// set the state for each region
                for (iregion=0; iregion<s->nregions; iregion++) {
                    if (s->regions[iregion] >= 0) {
                        RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                        rt->state = 0;
                        // coverage
                        int n = rt->align + rt->insertion + rt->deletion + rt->soft_clip + rt->known_snp;
                        // coverage - mark spare bins
                        if (n < s->region_min_count) rt->state |= REGION_STATE_COVERAGE;
                        // correct for sparse bins by assuming ALL bins have atleast region_min_count clusters
                        n = max(n, s->region_min_count);
                        // mismatch - mark bins with maximum mismatch rate > threshold
                        if (((float)rt->mismatch  / (float)n) >= s->region_mismatch_threshold)  rt->state |= REGION_STATE_MISMATCH;
                        // insertion - mark bins with maximum insertion rate > threshold
                        if (((float)rt->insertion / (float)n) >= s->region_insertion_threshold) rt->state |= REGION_STATE_INSERTION;
                        // deletion - mark bins with maximum deletion rate > threshold
                        if (((float)rt->deletion  / (float)n) >= s->region_deletion_threshold)  rt->state |= REGION_STATE_DELETION;
                    }
                }
				// ignoring low coverage, if all regions with a non-zero state have the same state and the
                // fraction of regions with this state exceeds a theshold set the state for the whole tile
				int tile_state = -1, nregions = 0;
                for (iregion=0; iregion<s->nregions; iregion++) {
                    if (s->regions[iregion] >= 0) {
                        RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                        int state = rt->state & ~REGION_STATE_COVERAGE;
                        if (!state) continue;
                        if (tile_state == -1) tile_state = state;
                        if (state != tile_state) break;
                        nregions++;
                    }
                }
				if (iregion == s->nregions && (((float)nregions/(float)s->nregions) >= TILE_REGION_THRESHOLD)) {
                    for (iregion=0; iregion<s->nregions; iregion++) {
                        if (s->regions[iregion] >= 0) {
                            RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                            rt->state = tile_state | (rt->state & REGION_STATE_COVERAGE);
                        }
                    }
                }
                if (!s->quiet) {
                    int mismatch = 0, insertion = 0, deletion = 0, soft_clip = 0;
                    long quality_bases = 0, quality_errors = 0;
                    for (iregion=0; iregion<s->nregions; iregion++) {
                        if (s->regions[iregion] >= 0) {
                            RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                            if (rt->state & REGION_STATE_MISMATCH)  mismatch++;
                            if (rt->state & REGION_STATE_INSERTION) insertion++;
                            if (rt->state & REGION_STATE_DELETION)  deletion++;
                            if (rt->state & REGION_STATE_SOFT_CLIP) soft_clip++;
                            quality_bases  += rt->align;
                            quality_errors += rt->mismatch;
                        }
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

/*
 * identify bad tiles metric < mean metric - filter * stdev metric
 *
 * where metric is one of mismatch, insertion, deletion, soft_clip or known_snp
 *
 * filter is typically 2.
*/

static void findBadTiles(Settings *s, int ntiles, size_t nreads, RegionTable ***rts)
{
	int itile, read, cycle, iregion;

	if (0 >= ntiles)
		return;

    for (itile=0; itile<ntiles; itile++) {
		int tile = s->tileArray[itile];
        for (read = 0; read < N_READS; read++) {
			for (cycle = 0; cycle < s->read_length[read]; cycle++) {
				// set the state for each region
                for (iregion=0; iregion<s->nregions; iregion++) {
                    if (s->regions[iregion] >= 0) {
                        RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                        rt->state = 0;
                        // coverage
                        int n = rt->align + rt->insertion + rt->deletion + rt->soft_clip + rt->known_snp;
                        // coverage - mark spare bins
                        if (n < s->region_min_count) rt->state |= REGION_STATE_COVERAGE;
                        // correct for sparse bins by assuming ALL bins have atleast region_min_count clusters
                        n = max(n, s->region_min_count);
                        // mismatch - mark bins with maximum mismatch rate > threshold
                        if (((float)rt->mismatch  / (float)n) >= s->region_mismatch_threshold)  rt->state |= REGION_STATE_MISMATCH;
                        // insertion - mark bins with maximum insertion rate > threshold
                        if (((float)rt->insertion / (float)n) >= s->region_insertion_threshold) rt->state |= REGION_STATE_INSERTION;
                        // deletion - mark bins with maximum deletion rate > threshold
                        if (((float)rt->deletion  / (float)n) >= s->region_deletion_threshold)  rt->state |= REGION_STATE_DELETION;
                    }
                }
				// ignoring low coverage, if all regions with a non-zero state have the same state and the
                // fraction of regions with this state exceeds a theshold set the state for the whole tile
				int tile_state = -1, nregions = 0;
                for (iregion=0; iregion<s->nregions; iregion++) {
                    if (s->regions[iregion] >= 0) {
                        RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                        int state = rt->state & ~REGION_STATE_COVERAGE;
                        if (!state) continue;
                        if (tile_state == -1) tile_state = state;
                        if (state != tile_state) break;
                        nregions++;
                    }
                }
				if (iregion < s->nregions&& (((float)nregions/(float)s->nregions) >= TILE_REGION_THRESHOLD)) {
                    for (iregion=0; iregion<s->nregions; iregion++) {
                        if (s->regions[iregion] >= 0) {
                            RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                            rt->state = tile_state | (rt->state & REGION_STATE_COVERAGE);
                        }
                    }
                }
                if (!s->quiet) {
                    int mismatch = 0, insertion = 0, deletion = 0, soft_clip = 0;
                    long quality_bases = 0, quality_errors = 0;
                    for (iregion=0; iregion<s->nregions; iregion++) {
                        if (s->regions[iregion] >= 0) {
                            RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                            if (rt->state & REGION_STATE_MISMATCH)  mismatch++;
                            if (rt->state & REGION_STATE_INSERTION) insertion++;
                            if (rt->state & REGION_STATE_DELETION)  deletion++;
                            if (rt->state & REGION_STATE_SOFT_CLIP) soft_clip++;
                            quality_bases  += rt->align;
                            quality_errors += rt->mismatch;
                        }
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

void printFilter(Settings *s, int ntiles, RegionTable ***rts) 
{
	FILE *fp;
	int itile, read, cycle, iregion;
	Header hdr;

    fp = fopen(s->filter, "w+");
	if (!fp) die("Can't open filter file %s: %s\n", s->filter, strerror(errno));

	hdr.region_magic = strdup(REGION_MAGIC);
	hdr.coord_shift = COORD_SHIFT; 
	hdr.coord_factor = COORD_FACTOR;
	hdr.ntiles = ntiles;
	hdr.tileArray = s->tileArray;
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
        
	free(hdr.region_magic);
	free(hdr.cmdLine);

    for (itile=0; itile<ntiles; itile++) {
		for (read = 0; read < N_READS; read++) {
			for (cycle = 0; cycle < s->read_length[read]; cycle++) {
                for (iregion=0; iregion<s->nregions; iregion++) {
                    int state = 0;
                    if (s->regions[iregion] >= 0) {
                        RegionTable *rt = rts[itile*N_READS+read][cycle] + s->regions[iregion];
                        state = rt->state;
                    }
					fputc(state, fp);
                }
			}
		}
	}

	fclose(fp);
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
			if (BAM_FREVERSE & bam->core.flag)
				oq = reverse_seq(oq);
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
			printf("\t%s", oq);
		}
		printf("\n");

		nreads_bam++;
	}

	bam_destroy1(bam);

	*nreads = nreads_bam;

	return 0;
}

static int findRegion(Settings *s, RegionTable ***rts, int ntiles, int x, int y)
{
    int iregion = -1;
    int ix = x2region(x, s->region_size);
    int iy = x2region(y, s->region_size);
    char key[100];
    char *cp;
    HashItem *hi;

    cp = append_int(key, ix);
    cp = append_char(cp, ':');
    cp = append_int(cp, iy);
    *cp = 0;
    if( NULL == (hi = HashTableSearch(s->region_hash, key, strlen(key))) ){
        HashData hd;
        iregion = s->region_hash->nused;
        hd.i = iregion;
        if( NULL == HashTableAdd(s->region_hash, key, strlen(key), hd, NULL) ) {
            fprintf(stderr, "ERROR: building rts hash table\n");
            exit(EXIT_FAILURE);
        }
        int nregions_x = (ix >= s->nregions_x ? (ix + 1) : s->nregions_x);
        int nregions_y = (iy >= s->nregions_y ? (iy + 1) : s->nregions_y);
        int nregions = nregions_x * nregions_y;
        if( nregions > s->nregions ){
            int itile, read, cycle;
            s->nregions_x = nregions_x;
            s->nregions_y = nregions_y;
            s->nregions = nregions;
            for(itile=0;itile<ntiles;itile++)
                for(read=0;read<N_READS;read++)
                {
                    if( NULL == rts[itile*N_READS+read]) continue;
                    for(cycle=0;cycle<s->read_length[read];cycle++)
                    {
                        int new_iregion;
                        rts[itile*N_READS+read][cycle] = srealloc(rts[itile*N_READS+read][cycle], s->nregions * sizeof(RegionTable));
                        for (new_iregion=iregion;new_iregion<s->nregions;new_iregion++) {
                            RegionTable *rt = rts[itile*N_READS+read][cycle] + new_iregion;
                            initialiseRegionTable(rt);
                        }
                    }
                }
        }
        iregion = hd.i;
    }else{
        iregion = hi->data.i;
    }

     return iregion;
}

static void updateRegionTable(Settings *s, RegionTable ***rts, int read, int iregion, int *read_qual, char *oq, int *read_mismatch)
{
    /* update region table */
	int cycle;
    for (cycle = 0; cycle < s->read_length[read]; cycle++) {
        RegionTable *rt = rts[read][cycle] + iregion;
        if (read_mismatch[cycle] & BASE_INSERTION) rt->insertion++;
        if (read_mismatch[cycle] & BASE_DELETION) rt->deletion++;
        if (read_mismatch[cycle] & BASE_SOFT_CLIP) rt->soft_clip++;
        if (read_mismatch[cycle] & BASE_KNOWN_SNP) { 
            rt->known_snp++;
        } else {
            if (read_mismatch[cycle] & BASE_ALIGN) rt->align++;
            if (read_mismatch[cycle] & BASE_MISMATCH) rt->mismatch++;
        }
        rt->quality1 += read_qual[cycle];
        if ( NULL != oq ) {
            rt->quality2 += (oq[cycle] - PHRED_QUAL_OFFSET);
        }
    }

	return;
}

/*
 * create an ordered array of tiles and re-order the RegionTable by tile
 */
RegionTable ***orderRegionTableByTile(Settings *s, int *tiles, int ntiles, RegionTable ***rts)
{
	int itile, read;

    // create a sorted array of tiles
	s->tileArray = smalloc(ntiles * sizeof(int));
	for (itile=0;itile<ntiles;itile++)
        s->tileArray[itile] = tiles[itile];
    qsort(s->tileArray, ntiles, sizeof(int), int_sort);

	// re-order the region table by tile
    RegionTable ***new_rts = smalloc(ntiles * N_READS * sizeof(RegionTable **));
	for (itile=0;itile<ntiles;itile++) {
	    int tile = s->tileArray[itile];
        size_t nelem = ntiles;
        void *pitile = lfind(&tile, tiles, &nelem, sizeof(int), &int_cmp);
        int old_itile = ((int*)pitile - tiles);
   	    for (read=0;read<N_READS;read++)
            new_rts[itile*N_READS+read] = rts[old_itile*N_READS+read];
    }

    free(rts);

    return new_rts;
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
RegionTable ***makeRegionTable(Settings *s, samfile_t *fp_bam, int *bam_ntiles, size_t *bam_nreads)
{
    RegionTable ***rts = NULL;

    int *tiles = NULL;
    int ntiles = 0;

	size_t nreads = 0;

	int lane = -1;

	static const int bam_read_buff_size = 1024;
	char bam_read_seq[bam_read_buff_size];
	int bam_read_qual[bam_read_buff_size];
	int bam_read_mismatch[bam_read_buff_size];

	bam1_t *bam = bam_init1();

	/* loop over reads in the bam file */
	while (1) {
		int bam_lane = -1, bam_tile = -1, bam_read = -1, bam_x = -1, bam_y = -1, read_length;

		if (parse_bam_readinfo(fp_bam, bam, &bam_lane, &bam_tile, &bam_x, &bam_y, &bam_read, NULL)) {
			break;	/* break on end of BAM file */
		}

		if (BAM_FUNMAP & bam->core.flag) continue;
#ifndef TILEVIZ
		if (BAM_FQCFAIL & bam->core.flag) continue;
		if (BAM_FPAIRED & bam->core.flag) {
			if (0 == (BAM_FPROPER_PAIR & bam->core.flag)) {
				continue;
			}
		}
#endif

		parse_bam_alignments(fp_bam, bam, bam_read_seq, bam_read_qual, NULL, bam_read_mismatch,
                                                  bam_read_buff_size, s->snp_hash);

		char *name = bam1_qname(bam);
		uint8_t *oq_ptr;
		char *oq = NULL;

		oq_ptr = bam_aux_get(bam, "OQ");
		if (NULL != oq_ptr) {
			oq = bam_aux2Z(oq_ptr);
			if (NULL == oq) {
				die("ERROR: Invalid original qualities %s for read: \"%s\"\n", oq, name);
			}
			if (BAM_FREVERSE & bam->core.flag)
				oq = reverse_seq(oq);
		}

		read_length = strlen(bam_read_seq);
		if (0 == s->read_length[bam_read]) {
			s->read_length[bam_read] = read_length;
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

        // lookup itile from tile in tile array
        size_t nelem = ntiles;
	    void *pitile;
	    int itile;
        pitile = lfind(&bam_tile, tiles, &nelem, sizeof(int), &int_cmp);
	    if (NULL == pitile) {
            int read;
            itile = ntiles;
            ntiles++;
            tiles = srealloc(tiles, ntiles * sizeof(int));
            tiles[itile] = bam_tile;
            rts = srealloc(rts, ntiles * N_READS * sizeof(RegionTable **));
            for (read=0;read<N_READS;read++)
                rts[itile*N_READS+read] = NULL;
            if (!s->quiet) fprintf(stderr, "Processing tile %i (%lu)\n", bam_tile, nreads);
        }else{
            itile = ((int*)pitile - tiles);
        }

        if (NULL == rts[itile*N_READS+bam_read]) {
            int cycle, iregion;
            rts[itile*N_READS+bam_read] = smalloc(read_length * sizeof(HashTable *));
            for(cycle=0;cycle<read_length;cycle++) {
                rts[itile*N_READS+bam_read][cycle] = smalloc(s->nregions * sizeof(RegionTable));
                for(iregion=0;iregion<s->nregions;iregion++) {
                    RegionTable *rt = rts[itile*N_READS+bam_read][cycle] + iregion;
                    initialiseRegionTable(rt);
                }
            }
        }

        int iregion = findRegion(s, rts, ntiles, bam_x, bam_y);
        updateRegionTable(s, &rts[itile*N_READS], bam_read, iregion, bam_read_qual, oq, bam_read_mismatch);

        nreads++;
	}

	bam_destroy1(bam);
	if (!s->quiet) display("nregions_x=%d nregions_y=%d nregions=%d\n", s->nregions_x, s->nregions_y, s->nregions);

    /* re-order the RegionTable by tile */
	rts = orderRegionTableByTile(s, tiles, ntiles, rts);

    /* setup a mapping between each potential region and the observed regions */
    regionMapping(s);

    *bam_ntiles = ntiles;
	*bam_nreads = nreads;

    return rts;
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
	fprintf(usagefp, "      -s --snp_file file\n");
	fprintf(usagefp, "                 set of snps to be removed\n");
	fprintf(usagefp, "                 file in Reference Ordered Data (ROD) format\n");
	fprintf(usagefp, "      --region_size\n");
	fprintf(usagefp, "                 default %d\n", REGION_SIZE);
	fprintf(usagefp, "      --region_mismatch_threshold\n");
	fprintf(usagefp, "                 threshold for setting region mismatch state\n");
	fprintf(usagefp, "                 default %-6.4f\n", REGION_MISMATCH_THRESHOLD);
	fprintf(usagefp, "      --region_insertion_threshold\n");
	fprintf(usagefp, "                 threshold for setting region insertion state\n");
	fprintf(usagefp, "                 default %-6.4f\n", REGION_INSERTION_THRESHOLD);
	fprintf(usagefp, "      --region_deletion_threshold\n");
	fprintf(usagefp, "                 threshold for setting region deletion state\n");
	fprintf(usagefp, "                 default %-6.4f\n", REGION_DELETION_THRESHOLD);
	fprintf(usagefp, "      --region_min_count\n");
	fprintf(usagefp, "                 minimum coverage when setting region state\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "                 If region_min_count = 0 an average region count will be\n");
    fprintf(usagefp, "                 calculated and the region_size increased until at least 2\n");
    fprintf(usagefp, "                 reads of each state are required for the region to pass the\n");
    fprintf(usagefp, "                 corresponding state threshold\n");
#ifdef TILEVIZ
	fprintf(usagefp, "      -t prefix\n");
	fprintf(usagefp, "                 generate tileviz files with this prefix\n");
#endif
	fprintf(usagefp, "\n");
	fprintf(usagefp, "    apply filter:\n");
	fprintf(usagefp, "      -o         output\n");
	fprintf(usagefp, "                 Output bam file name\n");
	fprintf(usagefp, "                 default: stdout\n");
	fprintf(usagefp, "      -u         do not compress the output bam file\n");
	fprintf(usagefp, "                 default: compress\n");
	fprintf(usagefp, "      -f         mark filtered reads as QCFAIL\n");
	fprintf(usagefp, "                 default: do not output filtered reads\n");
	fprintf(usagefp, "\n");

	exit(code);
}


void calculateFilter(Settings *s)
{
	samfile_t *fp_input_bam;
	int ntiles = 0;
	size_t nreads = 0;

	RegionTable ***rts = NULL;
    
	fp_input_bam = samopen(s->in_bam_file, "rb", 0);
	if (NULL == fp_input_bam) {
		die("ERROR: can't open bam file %s: %s\n", s->in_bam_file, strerror(errno));
	}

    s->region_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE|HASH_FUNC_JENKINS3);

	rts = makeRegionTable(s, fp_input_bam, &ntiles, &nreads);

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

    if( 0 == s->region_min_count ){
		display("Resizing regions\n");
        resizeRegions(s, ntiles, nreads, rts);
    }
    

    if( NULL != s->tileviz) {
	    tileviz(s, ntiles, rts);
    }
    
    if( s->calculate) {
        findBadRegions(s, ntiles, nreads, rts);

        if( NULL == s->filter) {
            display("Writing filter to stdout\n");
            s->filter = "/dev/stdout";
        }

        printFilter(s, ntiles, rts);
    }

    free(s->tileArray);
    HashTableDestroy(s->region_hash, 0);
    free(s->regions);
	freeRTS(s, ntiles, rts);
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
		die("ERROR: can't open bam file %s: %s\n", s->in_bam_file, strerror(errno));
	}

	out_bam_header = bam_header_dup(fp_input_bam->header);
	char concat_cmd[2048];
	strcpy(concat_cmd, hdr.cmdLine);
	strcat(concat_cmd, " ; ");
	strcat(concat_cmd, s->cmdline);
	bam_header_add_pg("spf", "spatial_filter", "A program to apply a spatial filter", concat_cmd, out_bam_header);

	out_bam_file = (NULL == s->output ? aprintf("/dev/stdout") : ((s->output)[0] == '/' ? aprintf("%s", s->output) : aprintf("%s/%s", s->working_dir, s->output)));
	fp_output_bam = samopen(out_bam_file, out_mode, out_bam_header);
	if (NULL == fp_output_bam) {
		die("ERROR: can't open bam file %s: %s\n", out_bam_file, strerror(errno));
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
			die("ERROR: can't open bam file %s: %s\n", s->in_bam_file, strerror(errno));
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

void dumpFilterFile(char *filename, int quiet)
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

	if (!quiet) {
        int tile, read, cycle, region;
   	    readFilterData(fp, &hdr);
        for (tile=0; tile<hdr.ntiles; tile++)
            for (read=0; read<hdr.nreads; read++)
                for (cycle=0; cycle<hdr.readLength[read]; cycle++)
                    for (region=0; region<hdr.nregions; region++) {
                        char state = getFilterData(hdr.tileArray[tile], read, cycle, region);
                        if (state & REGION_STATE_MASK)
                            printf("filtering tile=%d read=%d cycle=%d region=%d\n", hdr.tileArray[tile], read, cycle, region);
                    }
    }

	fclose(fp);
}


int main(int argc, char **argv)
{
	Settings settings;
	int dumpFilter = 0;
	char *cmd = NULL;

	settings.filter = NULL;
	settings.quiet = 0;
	settings.dump = 0;
	settings.tileviz = NULL;
	settings.calculate = 0;
	settings.apply = 0;
	settings.qcfail = 0;
	settings.snp_file = NULL;
	settings.snp_hash = NULL;
	settings.tileArray = NULL;
	settings.output = NULL;
	settings.in_bam_file = NULL;
	settings.read_length[0] = 0;
	settings.read_length[1] = 0;
	settings.read_length[2] = 0;
	settings.working_dir = NULL;
	settings.region_size = REGION_SIZE;
	settings.nregions_x = 0;
	settings.nregions_y = 0;
	settings.regions = NULL;
	settings.nregions = 0;
	settings.compress = 1;
	settings.region_min_count = REGION_MIN_COUNT;
	settings.region_mismatch_threshold = REGION_MISMATCH_THRESHOLD;
	settings.region_insertion_threshold = REGION_INSERTION_THRESHOLD;
	settings.region_deletion_threshold = REGION_DELETION_THRESHOLD;

	static struct option long_options[] = {
                   {"snp_file", 1, 0, 's'},
                   {"snp-file", 1, 0, 's'},
                   {"help", 0, 0, 'h'},
                   {"filter", 1, 0, 'F'},
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
#ifdef TILEVIZ
            case 't':	settings.tileviz = optarg; ncmd++; break;
#endif                
			case 'a':	settings.apply = 1; ncmd++; break;
			case 'f':	settings.qcfail = 1;		break;
			case 'u':	settings.compress = 0;		break;
			case 'o':	settings.output = optarg;	break;
			case 's':	settings.snp_file = optarg;	break;
			case 'F':	settings.filter = optarg;	break;
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
   	    if (settings.region_size < 1) die("Error: invalid region size\n");
    }

	// create pseudo command line
	if (settings.calculate) {
		char arg[64];
		cmd = smalloc(2048);
		strcat(cmd, argv[0]);
		strcat(cmd, " -c ");
		if (settings.snp_file)                   { strcat(cmd, " -s "); strcat(cmd, settings.snp_file); }
		if (settings.filter)                     { strcat(cmd, " -F "); strcat(cmd, settings.filter); }
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

	/* preserve starting directory */
	settings.working_dir = getcwd(NULL,0);
	if (NULL == settings.working_dir) {
		die("ERROR: can't obtain working directory: %s\n", strerror(errno));
	}

	/* dump the alignments */
	if (settings.dump) dumpBAM(&settings);

	/* Dump the filter file */
	if (dumpFilter) dumpFilterFile(settings.filter, settings.quiet);

	/* calculate the filter */
    if (settings.calculate || NULL != settings.tileviz) {
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
	if (NULL != settings.cmdline) free(settings.cmdline);

	return EXIT_SUCCESS;

}
