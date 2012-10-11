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
 * Parts of this code have been edited by Illumina.  These edits have been
 * made available under the following license:
 *
 * Copyright (c) 2008-2009, Illumina Inc.
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
 *     * Neither the name of the Illumina Inc. nor the
 *       names of its contributors may be used to endorse or promote
 *       products derived from this software without specific prior
 *       written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY ILLUMINA ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL ILLUMINA BE LIABLE FOR ANY DIRECT,
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
 * This code generates a purity based calibration table
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
 * CALDATA
 *   generate caldata.txt files cf. original version of pb_cal
 *
 * CHECK_BASECALL
 *   this option checks the base call corresponds to channel with
 *   the maximum intensity
 */

#define QC_FAIL
#define PROPERLY_PAIRED
//#define CALDATA
//#define CHECK_BASECALL

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

#define PBP_VERSION PACKAGE_VERSION

const char * pbs_c_rev = "$Revision$";

#define MAX_CIF_CHUNK_BYTES 4194304

/* if we split data by state, using a filter file rather than tile, use tile as a place holder for state */
#define N_STATES 2

#define BASE_ALIGN      (1<<0)
#define BASE_MISMATCH   (1<<1)
#define BASE_INSERTION  (1<<2)
#define BASE_DELETION   (1<<3)
#define BASE_SOFT_CLIP  (1<<4)
#define BASE_KNOWN_SNP  (1<<5)

#define ST_STATUS_GOOD  (1<<0)
#define ST_STATUS_BAD   (1<<1)

typedef struct {
    int         tile;
    int         read;
    int         cycle;
    int         nbins;
    float       offset;
    float       delta;
    float       scale;
    float       *purity;
    long        *num_bases;
    long        *num_errors;
    long        total_bases;
    float       quality;
    int         status;
} SurvTable;

typedef struct {
    int         tile;
    int         read;
    int         cycle;
    int         nbins;
    float       *purity;
    long        *num_bases;
    long        *num_errors;
    float       *frac_bases;
    float       *error_rate;
    float       *quality;
} CalTable;

typedef struct {
    long        lane;
    long        cycle;
    long        read;
    const char *dir;
} CifDir;

typedef struct {
    char *cmdline;
    char *prefix;
    char *intensity_dir;
    char *snp_file;
    HashTable *snp_hash;
    char *working_dir;
    CifDir *cif_dirs;
    size_t  n_cif_dirs;
    size_t *cif_lane_index;
    size_t  n_cif_lanes;
    int read_length[3];
    int cstart[3];
    int quiet;
    int filter_bad_tiles;
    int n_bins_left;
    int n_bins_right;
    int width;
    int height;
    int nregions_x;
    int nregions_y;
    int nregions;
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

static char *complement_table = NULL;

int xy2region(Settings *s, int x, int y)
{
	float x_coord = (float)(x - COORD_SHIFT) / (float)COORD_FACTOR;
	float y_coord = (float)(y - COORD_SHIFT) / (float)COORD_FACTOR;
	return (int)(x_coord / REGION_SIZE) * s->nregions_y + (int)(y_coord / REGION_SIZE);
}

static void setRegions(Settings * s)
{
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

	s->nregions_x = 1 + (int)(s->width / REGION_SIZE);
	s->nregions_y = 1 + (int)(s->height / REGION_SIZE);

	s->nregions = s->nregions_x * s->nregions_y;

        if (!s->quiet) display("nregions_x=%d nregions_y=%d nregions=%d\n", s->nregions_x, s->nregions_y, s->nregions);
        
	return;
}

static void readSnpFile(Settings *s)
{
    FILE *fp;
    HashTable *snp_hash;
    static const int line_size = 8192;
    char line[line_size];
    size_t count = 0;
    char last_chrom[100] = "";
    
    fprintf(stderr, "reading snp file %s\n", s->snp_file);

    fp = fopen(s->snp_file, "rb");
    if (NULL == fp) {
        fprintf(stderr, "ERROR: can't open known snp file %s: %s\n",
                s->snp_file, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if (NULL == (snp_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE|HASH_FUNC_JENKINS3))) {
        fprintf(stderr, "ERROR: creating snp hash table\n");
        exit(EXIT_FAILURE);
    }

    while( fgets(line, line_size, fp) )
    {
        char key[100];
        HashData hd;
        int bin, start, end;
        char chrom[100];
            
        if( 4 != sscanf(line, "%d\t%s\t%d\t%d", &bin, chrom, &start, &end) )
        {
            fprintf(stderr, "ERROR: reading snp file\n%s\n", line);
            exit(EXIT_FAILURE);
        }
        
        /* N.B rod start is 0 based */
        snprintf(key, sizeof(key), "%s:%d", chrom, start);
        hd.i = 0;
	if( NULL == HashTableAdd(snp_hash, key, strlen(key), hd, NULL) )
        {
            fprintf(stderr, "ERROR: building snp hash table\n");
            exit(EXIT_FAILURE);
	}

        if( strcmp(chrom, last_chrom) ){
            /*if( count ) fprintf(stderr, "read %lu snps on chrom %s\n", count, last_chrom);*/
            strcpy(last_chrom, chrom);
            count = 0;
        }

        count++;
    }

    /*if( count ) fprintf(stderr, "read %lu snps on chrom %s\n", count, last_chrom);*/

    fclose(fp);

    s->snp_hash = snp_hash;
}

static void initialiseSurvTable(Settings *s, SurvTable *st, int tile, int read, int cycle)
{
    int i;

    st->tile  = tile;
    st->read  = read;
    st->cycle = cycle;

    st->nbins = 76;

    st->offset = 0.25;
    st->delta  = 0.01;
    st->scale  = 1.0 / st->delta;
    st->purity = smalloc(st->nbins * sizeof(float));
    for (i=0;i<st->nbins;i++)
        st->purity[i] = st->offset + i * st->delta;

    st->num_bases  = smalloc(st->nbins * sizeof(long));
    st->num_errors = smalloc(st->nbins * sizeof(long));
    for (i=0;i<st->nbins;i++) {
        st->num_bases[i] = 0;
        st->num_errors[i] = 0;
    }

    st->total_bases = 0;

    st->status = ST_STATUS_GOOD;
}

static void freeSurvTable(Settings *s, SurvTable **sts)
{
    int itile, read, cycle;
    for(itile=0;itile<=N_TILES;itile++)
        for(read=0;read<N_READS;read++)
        {
            if( NULL == sts[itile*N_READS+read]) continue;
            for(cycle=0;cycle<s->read_length[read];cycle++)
            {
                SurvTable *st = sts[itile*N_READS+read] + cycle;
                if (st->nbins) {
                    free(st->purity);
                    free(st->num_bases);
                    free(st->num_errors);

                    st->nbins = 0;
                }
            }
        }
}

static void completeSurvTable(Settings *s, SurvTable **sts)
{
    float ssc = 1.0;
    int read, cycle, i;

    for(read=0;read<N_READS;read++)
    {
        if( NULL == sts[read]) continue;
        for(cycle=0;cycle<s->read_length[read];cycle++)
        {
            SurvTable *st = sts[read] + cycle;
            long quality_bases = 0;
            long quality_errors = 0;

            for(i=0;i<st->nbins;i++)
            {
                st->total_bases += st->num_bases[i];

                // bases with purity=0.25 are called as N and explicitly get a quality of 0
                if( st->purity[i] <= 0.25 )
                    continue;

                quality_bases  += st->num_bases[i];
                quality_errors += st->num_errors[i];
            }

            st->quality = -10.0 * log10((quality_errors + ssc)/(quality_bases + ssc));
        }
    }
}

/*
 * identify bad tiles quality < mean quality - filter * stdev quality
 * filter is typically 2.
*/
static void findBadTiles(Settings *s, int ntiles, SurvTable **sts)
{
    int filter = s->filter_bad_tiles;
    int read, cycle, itile;

    if (0 == filter || 0 >= ntiles)
        return;
    
    for(read=0;read<N_READS;read++)
        for(cycle=0;cycle<s->read_length[read];cycle++)
        {
            float qsum = 0.0;
            float q2sum = 0.0;
            float qmax = -1.0;
            int nq = 0;
            float qavg, q2avg, qstd, qmin;

            for(itile=0;itile<ntiles;itile++) {
                SurvTable *st = sts[itile*N_READS+read] + cycle;
                qmax = (st->quality > qmax ? st->quality : qmax);
            }

            /* it could be that all tiles for this cycle are all N */
            if( qmax < 0.0 )
                continue;

            /* flag tiles with q < qmax-20, i.e. error rate 100 * min error rate */
            qmin = max(qmax - 20.0, 0.0);

            for(itile=0;itile<ntiles;itile++) {
                SurvTable *st = sts[itile*N_READS+read] + cycle;
                /* skip tiles which are not good */
                if( st->status != ST_STATUS_GOOD )
                    continue;
                if( st->quality < qmin ) {
                    st->status = ST_STATUS_BAD;
                } else {
                    qsum += st->quality;
                    q2sum += st->quality * st->quality;
                    nq++;
                }
            }

            assert(nq > 0);

            qavg = qsum / nq;
            q2avg = q2sum / nq;
            qstd = sqrt(q2avg - qavg * qavg);
            qmin = qavg - filter * qstd;

            for(itile=0;itile<ntiles;itile++) {
                SurvTable *st = sts[itile*N_READS+read] + cycle;
                /* skip tiles which are not good */
                if (st->status != ST_STATUS_GOOD)
                    continue;
                if (st->quality < qmin)
                    st->status = ST_STATUS_BAD;
            }

#if 1
            fprintf(stderr, "read=%1d cycle=%-3d qavg=%.2f qstd=%.2f qmax=%.2f qmin=%.2f ", read, cycle, qavg, qstd, qmax, qmin);
            for(itile=0;itile<ntiles;itile++) {
                SurvTable *st = sts[itile*N_READS+read] + cycle;
                if (st->status == ST_STATUS_GOOD )
                    fprintf(stderr, "%4d  %5.2f ", st->tile, st->quality);
                else
                    // mark poor or bad tiles
                    fprintf(stderr, "%4d* %5.2f ", st->tile, st->quality);
            }
            fprintf(stderr, "\n");
#endif
        }

    return;
}

static void makeGlobalSurvTable(Settings *s, int ntiles, SurvTable **sts)
{
    int read, read_length, cycle, itile, i;

    if (0 >= ntiles)
        return;

    findBadTiles(s, ntiles, sts);

    for(read=0;read<N_READS;read++)
    {
        read_length = s->read_length[read];
        if (0 == read_length) continue;
        sts[ntiles*N_READS+read] = smalloc(read_length * sizeof(SurvTable));

        for(cycle=0;cycle<read_length;cycle++)
        {
            SurvTable *st = sts[ntiles*N_READS+read] + cycle;

            initialiseSurvTable(s, st, -1, read, cycle);

            for(itile=0;itile<ntiles;itile++)
            {
                SurvTable *tile_st = sts[itile*N_READS+read] + cycle;
                if (tile_st->status == ST_STATUS_GOOD )
                {
                    for(i=0;i<st->nbins;i++)
                    {
                        st->num_bases[i] += tile_st->num_bases[i];
                        st->num_errors[i] += tile_st->num_errors[i];
                    }
                    tile_st->total_bases = 0;
                }
            }
        }
    }

    completeSurvTable(s, &sts[ntiles*N_READS]);

    return;
}
    
static void outputSurvTable(Settings *s, SurvTable **sts)
{
    FILE *fp;
    int filename_sz;
    char *filename;
    int itile, read, cycle, i;

    filename_sz = (NULL == s->prefix ? 0 : strlen(s->prefix)) + 100;
    filename = smalloc(filename_sz);

    sprintf(filename, "%s_purity_cycle_surv.txt", s->prefix);
    fp = fopen(filename, "w");
    if( NULL == fp )
    {
        fprintf(stderr, "ERROR: can't open survival table file %s: %s\n",
                filename, strerror(errno));
        exit(EXIT_FAILURE);
    }

    free(filename);

    for(itile=0;itile<=N_TILES;itile++)
        for(read=0;read<N_READS;read++)
        {
            if( NULL == sts[itile*N_READS+read]) continue;
            for(cycle=0;cycle<s->read_length[read];cycle++)
            {
                SurvTable *st = sts[itile*N_READS+read] + cycle;

                // skip st with no data
                if( 0 == st->total_bases ) continue;

                for(i=0;i<st->nbins;i++)
                    fprintf(fp, "%.2f\t%d\t%d\t%d\t%ld\t%ld\n",
                            st->purity[i], st->read, st->cycle, st->tile,
                            st->num_bases[i], st->num_errors[i]);
            }
        }

    fclose(fp);
}

static void initialiseCalTable(SurvTable *st, CalTable *ct)
{
    int i;

    ct->tile  = st->tile;
    ct->read  = st->read;
    ct->cycle = st->cycle;

    ct->nbins = st->nbins;

    ct->purity = smalloc(ct->nbins * sizeof(float));
    for(i=0;i<ct->nbins;i++)
        ct->purity[i] = st->purity[i];

    ct->num_bases  = smalloc(ct->nbins * sizeof(long));
    ct->num_errors = smalloc(ct->nbins * sizeof(long));
    ct->frac_bases = smalloc(ct->nbins * sizeof(float));
    ct->error_rate = smalloc(ct->nbins * sizeof(float));
    ct->quality    = smalloc(ct->nbins * sizeof(float));
    for(i=0;i<ct->nbins;i++)
    {
        ct->num_bases[i]  = st->num_bases[i];
        ct->num_errors[i] = st->num_errors[i];
        ct->frac_bases[i] = 0.0;
        ct->error_rate[i] = 0.0;
        ct->quality[i]    = 0.0;
    }

    return;
}

static void freeCalTable(Settings *s, CalTable **cts)
{
    int itile, read, cycle;
    for(itile=0;itile<=N_TILES;itile++)
        for(read=0;read<N_READS;read++)
        {
            if( NULL == cts[itile*N_READS+read]) continue;
            for(cycle=0;cycle<s->read_length[read];cycle++)
            {
                CalTable *ct = cts[itile*N_READS+read] + cycle;
                if(ct->nbins) {
                    free(ct->purity);
                    free(ct->num_bases);
                    free(ct->num_errors);
                    free(ct->frac_bases);
                    free(ct->error_rate);
                    free(ct->quality);

                    ct->nbins = 0;
                }
            }
        }
}

static void optimisePurityBins(Settings *s, SurvTable *st, CalTable *ct)
{
    float ssc = 1.0;
    long cum_bases[76];
    long cum_errors[76];
    float purity_bins[76];
    int npurity_bins = 0;

    int i, j;
    int iqmax = 0, ipmax = 0;
    float max_diff = 0.0;
    int ipopt = 0, ipbin;
    float pinc;

    /* calc cum_bases and cum_errors */
    for(i=0;i<76;i++)
    {
        cum_bases[i] = 0;
        cum_errors[i] = 0;
        for(j=i;j<76;j++)
        {
            cum_bases[i]  += st->num_bases[j];
            cum_errors[i] += st->num_errors[j];
        }
    }
    /* find optimal purity by maximising diff_frac */
    for(i=0;i<76;i++)
    {
        float frac_bases, frac_errors, diff_frac;
        frac_bases = (float)cum_bases[i]/(float)cum_bases[0];
        frac_errors = (float)cum_errors[i]/(float)cum_errors[0];
        diff_frac = (frac_bases - frac_errors);
        if (diff_frac > max_diff) {
            max_diff = diff_frac;
            ipopt = i;
        }
    }

    /* find maximum (integer) quality */
    for(i=0;i<76;i++)
    {
        float error_rate, quality;
        int iq;
        error_rate = (cum_errors[i] + ssc)/(cum_bases[i] + ssc);
        quality = -10.0 * log10(error_rate);
        iq = (int)(quality + 0.5);
        if (iq > iqmax) {
            iqmax = iq;
            ipmax = i;
        }
    }

    /* if ipmax < ipopt reset ipmax to ipopt */
    if(ipmax < ipopt)
    {
        ipmax = ipopt;
        fprintf(stderr, "Resetting maximal bin to optimal bin read=%d cycle=%d tile=%d\n", st->read, st->cycle, st->tile);
    }

    /* first bin */
    ipbin = 0;
    purity_bins[npurity_bins++] = st->purity[ipbin];
#if 0
    fprintf(stderr,"1: npurity_bins=%d ipbin=%d purity=%f\n", npurity_bins, ipbin, st->purity[ipbin]);
#endif

    if (ipopt > 0)
    {
        /* add n_bins_left bins between first bin and popt */
        pinc = (st->purity[ipopt] - st->purity[0]) / (s->n_bins_left + 1);
        for(i=0,j=1;i<s->n_bins_left;i++,j++)
        {
            float p = st->purity[0] + j * pinc;
            int ip = st->scale * (p - st->offset) + 0.5;
            if (ip == ipopt)
                break;
            if (ip > ipbin) {
                ipbin = ip;
                purity_bins[npurity_bins++] = st->purity[ipbin];
#if 0
                fprintf(stderr,"2: npurity_bins=%d ipbin=%d purity=%f\n", npurity_bins, ipbin, st->purity[ipbin]);
#endif
            }
        }

        /* add the optimal bin */
        ipbin = ipopt;
        purity_bins[npurity_bins++] = st->purity[ipbin];
#if 0
        fprintf(stderr,"3: npurity_bins=%d ipbin=%d purity=%f\n", npurity_bins, ipbin, st->purity[ipbin]);
#endif
    }

    if (ipmax > ipopt)
    {
        /* add n_bins_right bins between popt and pmax */
        pinc = (st->purity[ipmax] - st->purity[ipopt]) / (s->n_bins_right + 1);
        for(i=0,j=1;i<s->n_bins_right;i++,j++)
        {
            float p = st->purity[ipopt] + j * pinc;
            int ip = st->scale * (p - st->offset) + 0.5;
            if (ip == ipmax)
                break;
            if (ip > ipbin) {
                ipbin = ip;
                purity_bins[npurity_bins++] = st->purity[ipbin];
#if 0
                fprintf(stderr,"4: npurity_bins=%d ipbin=%d purity=%f\n", npurity_bins, ipbin, st->purity[ipbin]);
#endif
            }
        }

        /* add pmax */
        ipbin = ipmax;
        purity_bins[npurity_bins++] = st->purity[ipbin];
#if 0
        fprintf(stderr,"5: npurity_bins=%d ipbin=%d purity=%f\n", npurity_bins, ipbin, st->purity[ipbin]);
#endif
    }

    /* add another bin if pmax is not the maximum possible purity value */
    if (ipmax < 75)
        purity_bins[npurity_bins++] = st->purity[75];
#if 0
    if (ipmax < 75)
        fprintf(stderr,"6: npurity_bins=%d ipbin=%d purity=%f\n", npurity_bins, 75, st->purity[75]);
#endif

    /* move data into the new set of bins */

    ct->nbins = npurity_bins;

    for(i=0;i<ct->nbins;i++)
    {
        ct->purity[i] = purity_bins[i];
        ct->num_bases[i]  = 0;
        ct->num_errors[i] = 0;
    }

    for(i=0;i<st->nbins;i++)
    {
        for(j=0;j<ct->nbins;j++)
        {
            if(st->purity[i] <= ct->purity[j])
            {
                ct->num_bases[j]  += st->num_bases[i];
                ct->num_errors[j] += st->num_errors[i];
                break;
            }
        }
    }
}

static int makeCalTable(Settings *s, SurvTable **sts, CalTable **cts)
{
    int nct = 0;
    float ssc = 1.0;
    int itile, read, read_length, cycle, i;

    for(itile=0;itile<=N_TILES;itile++)
        for(read=0;read<N_READS;read++)
        {
            cts[itile*N_READS+read] = NULL;
            if( NULL == sts[itile*N_READS+read]) continue;

            read_length = s->read_length[read];
            cts[itile*N_READS+read] = smalloc(read_length * sizeof(CalTable));
            nct += read_length;
            
            for(cycle=0;cycle<read_length;cycle++)
            {
                SurvTable *st = sts[itile*N_READS+read] + cycle;
                CalTable *ct  = cts[itile*N_READS+read] + cycle;

                // set number of bins in ct to 0
                ct->nbins = 0;

                // skip st with no data
                if( 0 == st->total_bases ) continue;

                initialiseCalTable(st, ct);

                optimisePurityBins(s, st, ct);

                for(i=0;i<ct->nbins;i++)
                {
                    ct->frac_bases[i] = ((float)ct->num_bases[i]) / ((float)st->total_bases);
                    ct->error_rate[i] = (ct->num_errors[i] + ssc)/(ct->num_bases[i] + ssc);
                    ct->quality[i] = -10.0 * (log10(ct->error_rate[i]));
                }
            }
        }

    return nct;
}

static void outputCalTable(Settings *s, CalTable **cts)
{
    FILE *fp;
    int filename_sz;
    char *filename;
    int itile, read, cycle, i;

    filename_sz = (NULL == s->prefix ? 0 : strlen(s->prefix)) + 100;
    filename = smalloc(filename_sz);

    sprintf(filename, "%s_purity_cycle_caltable.txt", s->prefix);
    fp = fopen(filename, "w");
    if (NULL == fp) {
        fprintf(stderr, "ERROR: can't open CT file %s: %s\n",
                filename, strerror(errno));
        exit(EXIT_FAILURE);
    }

    free(filename);

    for(itile=0;itile<=N_TILES;itile++)
        for(read=0;read<N_READS;read++)
        {
            if( NULL == cts[itile*N_READS+read]) continue;
            for(cycle=0;cycle<s->read_length[read];cycle++)
            {
                CalTable *ct = cts[itile*N_READS+read] + cycle;

                // skip ct with no bins
                if( 0 == ct->nbins ) continue;

                for(i=0;i<ct->nbins;i++)
                {
                    /* exclude empty bins except for purity=0.25 and purity=1.0 */
                    if( ct->num_bases[i] == 0 && ( ct->purity[i] > 0.25 && ct->purity[i] < 1.0 ))
                        continue;

                    fprintf(fp, "%f\t%d\t%d\t%d\t%5ld\t%5ld\t%f\t%f\t%f\n",
                            ct->purity[i], ct->read, ct->cycle, ct->tile,
                            ct->num_bases[i], ct->num_errors[i], ct->frac_bases[i],
                            ct->error_rate[i], ct->quality[i]);
                }
            }
        }

    fclose(fp);
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
        if(Imax <= data[i])
            Imax = data[i];
        Isum += data[i];
    }
    if (Imin <= 0)
    {
        Imin--;
        Imax -= Imin;
        Isum -= 4 * Imin;
    }

    return (float)Imax/(float)Isum;
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

/**
 * Return a character representing the complement-base of the supplied
 * parameter.  The mapping is as follows: a->t, c->g, g->c, t|u->a, [->],
 * ]->[, -->-, all others->n. There is a single shot initalisation
 * of a static lookup table to represent this mapping. In addition the
 * case of the supplied parameter is preserved, Uppercase->Uppercase and
 * vice-versa.
 *
 * @param c is the character representing a base. Mapped values include:
 * {a,c,g,t,u,[,],-} all other inputs get a default mapping of 'n'.
 *
 * @returns the character representation of the bilogical compliment of the
 * supplied base. 
 */
char
complement_base(char c)
{

    if (!complement_table) {
	int x;
	complement_table = (char *) calloc(256, sizeof(char)) + 127;

	for (x = -127; x < 128; x++) {
	    if (x == 'a')
		complement_table[x] = 't';
	    else if (x == 'c')
		complement_table[x] = 'g';
	    else if (x == 'g')
		complement_table[x] = 'c';
	    else if (x == 't' || x == 'u')
		complement_table[x] = 'a';
	    else if (x == 'n')
		complement_table[x] = 'n';
	    else if (x == 'A')
		complement_table[x] = 'T';
	    else if (x == 'C')
		complement_table[x] = 'G';
	    else if (x == 'G')
		complement_table[x] = 'C';
	    else if (x == 'T' || x == 'U')
		complement_table[x] = 'A';
	    else if (x == 'N')
		complement_table[x] = 'N';
	    else
		complement_table[x] = x;
	}
    }

    return complement_table[(int) c];
}

/**
 * Convert a string containing a string representation of a base sequence
 * into its biological complement. See complement_base for the mapping. This
 * does not reverses the direction of the string.
 *
 * @see complement_base
 *
 * @param seq is the input sequence. <b>NB.</b> this is destructively modified.
 *
 * @returns a pointer to the original string storage but with the contents
 * modified to have complement characters.
 */
char *
complement_seq(char *seq)
{
    char *s = seq;

    while (*s) {
        *s = complement_base(*s);

        s++;
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
                    int *bam_read,
                    char* read_seq,
                    int* read_qual,
                    char* read_ref,
                    int* read_mismatch,
                    const int read_buff_size) {

    char *name;
    int32_t pos;
    uint32_t *cigar;
    uint8_t *seq, *qual, *m_ptr, *ci_ptr;
    char *mismatch;
    const char *sep = ":#/", *sep2 = "ACGTN^" ;
    const char *cp, *cp2;
    int lane, tile, x, y, offset, read, i, j;
    int skip = 0, clip, insert, delete;
    HashItem *hi;

    if(0 == read_buff_size) return 1;

    read_seq[0] = 0;

    if( 0 > samread(fp, bam)) return 1;
    
    if (BAM_FUNMAP & bam->core.flag)
        return 0;

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

#ifdef QC_FAIL
    if(BAM_FQCFAIL & bam->core.flag)
        return 0;
#endif
    
#ifdef PROPERLY_PAIRED
    if(BAM_FPAIRED & bam->core.flag)
        if(0 == (BAM_FPROPER_PAIR & bam->core.flag))
            return 0;
#endif
    
    pos   = bam->core.pos;
    cigar = bam1_cigar(bam);
    seq   = bam1_seq(bam);
    qual  = bam1_qual(bam);
    m_ptr = bam_aux_get(bam, "MD");

    if(NULL == m_ptr){
        fprintf(stderr,"ERROR: No mismatch for read: \"%s\"\n",name);
        exit(EXIT_FAILURE);
    }else{
        mismatch = bam_aux2Z(m_ptr);
        if(NULL == mismatch){
            fprintf(stderr,"ERROR: Invalid mismatch %s for read: \"%s\"\n",mismatch,name);
            exit(EXIT_FAILURE);
        }
    }

    memset(read_mismatch, 0, bam->core.l_qseq * sizeof(int));

    for (i = 0; i < bam->core.l_qseq; i++) {
        read_seq[i] =  bam_nt16_rev_table[bam1_seqi(seq, i)];
        read_qual[i] = qual[i];
    }
    read_seq[i] = 0;
    read_ref[i] = 0;
    
    j = 0;
    for (i=0; i<bam->core.n_cigar; i++) {
        int l = cigar[i] >> 4, op = cigar[i] & 0xf, k;
        switch(op) {
        case BAM_CMATCH:
            // CIGAR: alignment match;
            for(k=0; k<l; j++, k++)
            {
                read_mismatch[j] |= BASE_ALIGN;
                read_ref[j] = read_seq[j];
            }
            break;
        case BAM_CDEL:
            // CIGAR: deletion from the reference
            if(j == bam->core.l_qseq){
                fprintf(stderr,"ERROR: Trailing deletion for read: %s\n", name);
                exit(EXIT_FAILURE);
            }
            read_mismatch[j] |= BASE_DELETION;
            break;
        case BAM_CINS:
            // CIGAR: insertion to the reference 
            for(k=0; k<l; j++, k++)
            {
                read_mismatch[j] |= BASE_INSERTION;
                read_ref[j] = 'I';
            }
            break;
        case BAM_CSOFT_CLIP:
            // CIGAR: clip on the read with clipped sequence present in qseq
            for(k=0; k<l; j++, k++)
            {
                read_mismatch[j] |= BASE_SOFT_CLIP;
                read_ref[j] = 'S';
            }
            break;
        default:
            fprintf(stderr,"ERROR: Unexpected CIGAR operation: %c\n", op);
            exit(EXIT_FAILURE);
        }
    }
    if(j != bam->core.l_qseq){
        fprintf(stderr,"ERROR: Inconsistent cigar string %d > %d for read: \"%s\"\n", j, bam->core.l_qseq, name);
        exit(EXIT_FAILURE);
    }

    /* clipped sequence is missing from MD */
    for(i=0,clip=0;i<bam->core.l_qseq;i++,clip++)
        if( 0 == (read_mismatch[i] & BASE_SOFT_CLIP))
            break;
    if(0 == (read_mismatch[i] & BASE_ALIGN)){
        fprintf(stderr,"ERROR: Inconsistent cigar string expect alignment after soft clip for read: \"%s\"\n", name);
        exit(EXIT_FAILURE);
    }
    if(clip )
        skip += clip;

    cp2 = mismatch;
    while( NULL != (cp2 = parse_next_int(cp2,&j,sep2)) ){
        /* skip matching bases, exclude insertions which are missing from MD */
        for(insert=0;j>0;i++)
            if(read_mismatch[i] & BASE_INSERTION)
                insert++;
            else
                j--;
        if(insert)
            skip += insert;

        if(0 == strlen(cp2))
            /* reached end of MD string */
            break;

        /* skip insertions which are missing from MD */
        for(insert=0;i<bam->core.l_qseq;i++,insert++)
            if(0 == (read_mismatch[i] & BASE_INSERTION))
                break;
        if (i == bam->core.l_qseq)
        {
            fprintf(stderr,"ERROR: Invalid MD string %s for read: \"%s\"\n", mismatch, name);
            exit(EXIT_FAILURE);
        }
        if(insert)
            skip += insert;

        switch(*cp2) {
        case '^':
            /* deletions missing from read_seq */
            delete = 0;
            while( NULL != strchr("ACGTN", *(++cp2)) )
                delete++;
            if(delete)
                skip -= delete;
            break;
        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'N':
            /* mismatch */
            if(0 == (read_mismatch[i] & BASE_ALIGN)){
                fprintf(stderr,"ERROR: Inconsistent cigar string expect alignment at mismatch for read: \"%s\"\n", name);
                exit(EXIT_FAILURE);
            }
            read_ref[i] = *cp2;
            read_mismatch[i] |= BASE_MISMATCH;
            
            pos = bam->core.pos + (i-skip);

            if( NULL != s->snp_hash ){
                char *chrom = fp->header->target_name[bam->core.tid];
                char key[100];
                /* N.B bam->core.pos is 0 based */
                snprintf(key, sizeof(key), "%s:%d", chrom, pos);
                if (NULL != (hi = HashTableSearch(s->snp_hash, key, strlen(key))))
                {
                    hi->data.i++;
                    read_mismatch[i] |= BASE_KNOWN_SNP;
                }
            }
            i++;
            break;
        default:
            fprintf(stderr,"ERROR: Invalid mismatch %s(%c)\n", mismatch, *cp2);
            exit(EXIT_FAILURE);
        }
    }

    /* clipped sequence is missing from MD */
    for(clip=0;i<bam->core.l_qseq;i++,clip++)
        if(0 == (read_mismatch[i] & BASE_SOFT_CLIP))
            break;
    if(clip )
        skip += clip;
    if (i != bam->core.l_qseq)
    {
        fprintf(stderr,"ERROR: Inconsistent MD string %d != %d for read: \"%s\"\n", i, bam->core.l_qseq, name);
        exit(EXIT_FAILURE);
    }

    if (BAM_FREVERSE & bam->core.flag) {
        read_seq = reverse_seq(read_seq);
        read_seq = complement_seq(read_seq);
        read_qual = reverse_int(read_qual, bam->core.l_qseq);
        read_ref = reverse_seq(read_ref);
        read_ref = complement_seq(read_ref);
        read_mismatch = reverse_int(read_mismatch, bam->core.l_qseq);
    }
    
    *bam_lane = lane;
    *bam_tile = tile;
    *bam_x = x;
    *bam_y = y;
    *bam_offset = offset;
    *bam_read = read;

    return 0;
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

    if (NULL == cif_data->cycles)
        return;
    
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
    
    free(cif_data);
}

static int updateSurvTable(Settings *s, SurvTable **sts, CifData *cif_data,
                           size_t spot_num, int tile, int x, int y, int read, int *read_mismatch,
                           char *read_seq, int *read_qual, char *read_ref, samfile_t *fp, bam1_t *bam, FILE *fp_caldata) {

    int cstart = s->cstart[read];
    int read_length = s->read_length[read];
    int iregion = -1;
    int c, b;

    assert((cstart + read_length) <= cif_data->ncycles);

    if (s->nregions) iregion = xy2region(s, x, y);

#ifdef CALDATA
    {
        int dir = BAM_FREVERSE & bam->core.flag ? 1 : 0;
        char *chrom = fp->header->target_name[bam->core.tid];
        int pos = bam->core.pos;

        fprintf(fp_caldata, "%d\t%d\t%d\t%s\t%d\t", x, y, dir, chrom, pos);
    }
#endif

    /* update survival table */
    for (c = cstart, b = 0; b < read_length; c++, b++) {
        CifCycleData *cycle = cif_data->cycles + c;
        SurvTable *st;
        int channel;
        int bin[cycle->num_channels];
        float purity = -1.0;
        int ibin;

        /* set cycle ct */
        if (s->nregions) {
            int state = (getFilterData(tile, read, b, iregion) & REGION_STATE_MISMATCH) ? 1 : 0;
            st = sts[state*N_READS+read] + b;
        }else{
            st = sts[read] + b;
        }

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

        ibin = st->scale * (purity - st->offset) + 0.5;

        if( read_mismatch[b] & BASE_KNOWN_SNP ) {
            // don't count these
        } else{
            if( read_mismatch[b] & BASE_ALIGN )
                st->num_bases[ibin]++;
            if( read_mismatch[b] & BASE_MISMATCH )
                st->num_errors[ibin]++;
        }

#ifdef CALDATA
        {
            char bases[] = "ACGTNDIS";
            char base = strchr(bases, read_seq[b]);
            int seq_base, ref_base;

            if( NULL == base ){
                fprintf(stderr, "Unknown base %c\n", read_seq[b]);
                exit(EXIT_FAILURE);
            }
            seq_base = base - bases + 1;
            base = strchr(bases, read_ref[b]);
            if( NULL == base ){
                fprintf(stderr, "Unknown ref %c\n", read_ref[b]);
                exit(EXIT_FAILURE);
            }
            ref_base = base - bases + 1;
            fprintf(fp_caldata, "%d %d %d %d %d %d %d %d\t",
                    ref_base, seq_base, read_qual[b], bin[0], bin[1], bin[2], bin[3], read_mismatch[b]);
        }
#endif

#ifdef CHECK_BASECALL
        {
            int Imax = bin[0];
            int i;
            for (i=1; i<cycle->num_channels; i++)
                if(Imax <= bin[i])
                    Imax = bin[i];
            if (Imax > 0 ){
                char bases[] = "ACGTN", flags[] = "    ";
                for (i=0; i<cycle->num_channels; i++)
                {
                    flags[i] = ' ';
                    if(bin[i] == Imax)
                    {
                        flags[i] = '*';
                        if(read_seq[b] == bases[i])
                            break;
                    }
                }
                if( i == cycle->num_channels )
                    fprintf(stderr,"%lu %d %d %d %c %c %d %d%c %d%c %d%c %d%c %f\n", spot_num, x, y, c, read_ref[b], read_seq[b], read_qual[b],
                            bin[0], flags[0], bin[1], flags[1], bin[2], flags[2], bin[3], flags[3], purity);
            }
        }
#endif

    }
    
#ifdef CALDATA
    fprintf(fp_caldata, "\n");
#endif

    return 0;
}

/*
 * Takes the bam file as input and updates the survival table
 * It uses the associated DIF files too to do this.
 *
 * Assumption: within a single input file, all reads are the same length and
 * we're using unclipped data.
 *
 * Returns: 0 written for success
 *	   -1 for failure
 */
int makeSurvTable(Settings *s, samfile_t *fp_bam, SurvTable **sts, int *ntiles, int *nreads) {

    int nst = 0;

    CifData *cif_data = NULL;

    int ncycles_firecrest = -1;

    int lane = -1;
    int tile = -1;

    int tiles[N_TILES];

    int ntiles_bam = 0;
    int nreads_bam = 0;

    static const int bam_read_buff_size = 1024;
    char bam_read_seq[bam_read_buff_size];
    int bam_read_qual[bam_read_buff_size];
    char bam_read_ref[bam_read_buff_size];
    int bam_read_mismatch[bam_read_buff_size];
    FILE *fp_caldata = NULL;

    bam1_t *bam = bam_init1();

    int itile, read;

    checked_chdir(s->intensity_dir);

    for(itile=0;itile<=N_TILES;itile++)
        for(read=0;read<N_READS;read++)
            sts[itile*N_READS+read] = NULL;

    itile = -1;

    /* loop over reads in the bam file */
    while (1){
        int bam_lane = -1, bam_tile = -1, bam_x = -1, bam_y = -1, bam_read = -1, read_length;
        size_t bam_offset = 0;
        int cycle;

        if (0 != parse_bam_file_line(s, fp_bam, bam,
                                     &bam_lane, &bam_tile, &bam_x, &bam_y, &bam_offset, &bam_read,
                                     bam_read_seq, bam_read_qual, bam_read_ref,
                                     bam_read_mismatch, bam_read_buff_size)) {
            break;
        }

        read_length = strlen(bam_read_seq);
        if (0 == read_length) continue;

        if (0 == s->read_length[bam_read]) {
            if(read_length >= N_CYCLES){
                fprintf(stderr,
                        "ERROR: too many cycles for read %d "
                        "in bam file %d > %d.\n",
                        bam_read, read_length, N_CYCLES);
                exit(EXIT_FAILURE);
            }
            s->read_length[bam_read] = read_length;

            if (s->nregions) {
                /* if we have a filter file we allocate memory for each state now */
                int state;
                for(state=0;state<N_STATES;state++) {
                    sts[state*N_READS+bam_read] = smalloc(read_length * sizeof(SurvTable));
                    nst += read_length;
                    for(cycle=0;cycle<read_length;cycle++)
                        initialiseSurvTable(s, sts[state*N_READS+bam_read]+cycle, state, bam_read, cycle);
                }
            }
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
            fprintf(stderr,
                    "Error: Inconsistent lane "
                    "within bam file.\n"
                    "have %d, previously it was %d.\n",
                    bam_lane, lane);
            exit(EXIT_FAILURE);
        }

        if( bam_tile != tile ){
            tile = bam_tile;
            
            if (!s->quiet)fprintf(stderr, "Processing tile %i (%d)\n", tile, nreads_bam);

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

            for (itile=0; itile<ntiles_bam; itile++)
                if (tile == tiles[itile]) {
                    fprintf(stderr,"ERROR: alignments are not sorted by tile.\n");
                    exit(EXIT_FAILURE);
                }

            ntiles_bam++;
            if(ntiles_bam > N_TILES){
                fprintf(stderr,"ERROR: too many tiles %d > %d.\n", ntiles_bam, N_TILES);
                exit(EXIT_FAILURE);
            }

            tiles[itile] = tile;

#ifdef CALDATA
            if (NULL != fp_caldata) fclose(fp_caldata);
            if (!s->quiet)fprintf(stderr, "Opening caldata file %s\n", filename);
            fp_caldata = fopen(filename, "w");
            if (NULL == fp_caldata) {
                fprintf(stderr, "ERROR: can't open caldata file %s: %s\n",
                        filename, strerror(errno));
                exit(EXIT_FAILURE);
            }
#endif            
        }

        if (0 == s->nregions) {
            /* if we don't have a filter file we allocate memory for this tile now */
            if (NULL == sts[itile*N_READS+bam_read]) {
                sts[itile*N_READS+bam_read] = smalloc(read_length * sizeof(SurvTable));
                nst += read_length;
                for(cycle=0;cycle<read_length;cycle++)
                    initialiseSurvTable(s, sts[itile*N_READS+bam_read]+cycle, tile, bam_read, cycle);
            }
        }
        

        if (0 != updateSurvTable(s, (s->nregions ? sts : &sts[itile*N_READS]), cif_data,
                                 bam_offset, bam_tile, bam_x, bam_y, bam_read, bam_read_mismatch,
                                 bam_read_seq, bam_read_qual, bam_read_ref, fp_bam, bam, fp_caldata)) {
            fprintf(stderr,"ERROR: updating quality values for tile %i.\n", tile);
            exit(EXIT_FAILURE);
        }
        
        nreads_bam++;
    }
    
    if (s->nregions) {
        int state;
        for(state=0;state<N_STATES;state++)
            completeSurvTable(s, &sts[state*N_READS]);
    }else{
        for(itile=0;itile<ntiles_bam;itile++)
            completeSurvTable(s, &sts[itile*N_READS]);
    }
    

    if (NULL != cif_data) free_cif_data(cif_data);

    if (NULL != fp_caldata) fclose(fp_caldata);

    bam_destroy1(bam);
    
    *ntiles = ntiles_bam;
    *nreads = nreads_bam;

    return nst;
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
char *
get_real_file_path_name(const char* file) {

    char *tmp_path = sstrdup(file);
    char *cp;
    char *out_path;

    if (NULL != (cp = strrchr(tmp_path, '/'))) {
        *(cp+1) = 0;

        out_path = get_real_path_name(tmp_path);

    } else {
        out_path = alloc_getcwd();
    }

    free(tmp_path);

    return out_path;
}



static
void usage(int code) {
    FILE* usagefp = stderr;

    fprintf(usagefp, "pb_calibration v" PBP_VERSION "\n\n");
    fprintf(usagefp, 
            "Usage: pb_calibration [options] bam_file\n"
            ""
            "  outputs the calibration table to a file called caltable.txt\n"
            "    unless a prefix is specified with the -prefix flag in which case the\n"
            "    file is called prefix_caltable.txt\n"
            "");
    fprintf(usagefp, "  options:\n");
    fprintf(usagefp, "    -snp_file file\n");
    fprintf(usagefp, "               set of snps to be removed be calibration\n");
    fprintf(usagefp, "                 file in Reference Ordered Data (ROD) format\n");
    fprintf(usagefp, "    -p prefix  Output file prefix e.g. prefix_caltable.txt\n");
    fprintf(usagefp, "                 default no prefix\n");
    fprintf(usagefp, "    -intensity-dir dir\n");
    fprintf(usagefp, "               Intensity directory\n");
    fprintf(usagefp, "    -filter_file file\n");
    fprintf(usagefp, "               spatial filter file\n");
    fprintf(usagefp, "    -filter-bad-tiles threshold\n");
    fprintf(usagefp, "               filter tiles with q < qavg - threshold * qstd\n");
    fprintf(usagefp, "                 will output a list of good/bad tiles\n");
    fprintf(usagefp, "                 to good/bad_tiles.lst (or prefix_good/bad_tiles.lst)\n");
    fprintf(usagefp, "    -cstart int\n");
    fprintf(usagefp, "             intensity cycle number of first base of read in single-end bam file\n");
    fprintf(usagefp, "               no default\n");
    fprintf(usagefp, "    -cstart1 int\n");
    fprintf(usagefp, "             intensity cycle number of first base of read 1 in paired-end bam file\n");
    fprintf(usagefp, "               no default\n");
    fprintf(usagefp, "    -cstart2 int\n");
    fprintf(usagefp, "             intensity cycle number of first base of read 2 in paired-end bam file\n");
    fprintf(usagefp, "               no default\n");
    fprintf(usagefp, "    -nL nbins\n");
    fprintf(usagefp, "               number of purity bins between purity 0.25 and optimal cut\n");
    fprintf(usagefp, "                 default 2\n");
    fprintf(usagefp, "    -nR nbins\n");
    fprintf(usagefp, "               number of purity bins between optimal cut and qmax\n");
    fprintf(usagefp, "                 default 2\n");
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
    char *bam_file = NULL;
    samfile_t *fp_bam = NULL;
    const char *override_intensity_dir = NULL;
    const char *filter_file = NULL;
    Header filter_header;
    int ntiles = 0;
    int nreads = 0;
    int nst = 0;
    SurvTable *sts[(N_TILES+1)*N_READS];
    int nct = 0;
    CalTable *cts[(N_TILES+1)*N_READS];

    settings.prefix = NULL;
    settings.quiet = 0;
    settings.filter_bad_tiles = 0;
    settings.n_bins_left = 2;
    settings.n_bins_right = 2;
    settings.cstart[0] = 0;
    settings.cstart[1] = 0;
    settings.cstart[2] = 0;
    settings.intensity_dir = NULL;
    settings.snp_file = NULL;
    settings.snp_hash = NULL;
    settings.read_length[0] = 0;
    settings.read_length[1] = 0;
    settings.read_length[2] = 0;
    settings.working_dir = NULL;
    settings.cif_dirs       = NULL;
    settings.n_cif_dirs     = 0;
    settings.cif_lane_index = NULL;
    settings.n_cif_lanes    = 0;
    settings.width      = 0;
    settings.height     = 0;
    settings.nregions   = 0;
    settings.nregions_x = 0;
    settings.nregions_y = 0;
    
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

	} else if (!strcmp(argv[i], "-snp_file")) {
            if(settings.snp_file != NULL) {
		fprintf(stderr, "ERROR: -snp_file specified multiple times\n");
                usage(1);
            }
            check_arg(i,argc,"-snp_file");
            settings.snp_file = argv[++i];
	} else if (!strcmp(argv[i], "-filter_file")) {
            if(filter_file != NULL) {
		fprintf(stderr, "ERROR: -filter_file specified multiple times\n");
                usage(1);
            }
            check_arg(i,argc,"-filter_file");
            filter_file = argv[++i];

	} else if (!strcmp(argv[i], "-q")) {
	    settings.quiet = 1;
	} else if (!strcmp(argv[i], "-p")) {
            if(settings.prefix != NULL) {
		fprintf(stderr, "ERROR: -p option specified multiple times\n");
                usage(1);
            }
            check_arg(i,argc,"-p");
            settings.prefix = argv[++i];

	} else if (!strcmp(argv[i], "-filter-bad-tiles")){
            settings.filter_bad_tiles = atoi(argv[++i]);
            if(settings.filter_bad_tiles < 1){
                fprintf(stderr,"ERROR: invalid argument to -filter_bad_tiles\n");
                usage(1);
            }
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
	} else if (!strcmp(argv[i], "-nL")){
            check_arg(i,argc,"-nL");
            settings.n_bins_left = atoi(argv[++i]);
            if(settings.n_bins_left < 0){
                fprintf(stderr,"ERROR: invalid argument to -nL\n");
                usage(1);
            }
	} else if (!strcmp(argv[i], "-nR")){
            check_arg(i,argc,"-nR");
            settings.n_bins_right = atoi(argv[++i]);
            if(settings.n_bins_right < 0){
                fprintf(stderr,"ERROR: invalid argument to -nR\n");
                usage(1);
            }

	} else if (!strcmp(argv[i], "-h")) {
	    usage(0);
	} else {
            fprintf(stderr,"ERROR: Unknown option %s\n", argv[i]);
	    usage(1);
	}
    }

    if ((argc-i) < 1)
	usage(0);

    /* preserve starting directory b/c makeSurvTable is going to chdir all over the place */
    settings.working_dir = alloc_getcwd();
    if (NULL == settings.working_dir) {
        fprintf(stderr, "ERROR: can't obtain working directory: %s\n",
                strerror(errno));
        exit(EXIT_FAILURE);
    }

    /* get absolute intensity dir*/
    if (override_intensity_dir) {
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

    /* read the snp_file */
    if (NULL != settings.snp_file) {
        readSnpFile(&settings);
        if (NULL == settings.snp_hash) {
            fprintf(stderr, "ERROR: reading snp file %s\n", settings.snp_file);
            exit(EXIT_FAILURE);
        }
    }

    /* read filter file */
    if (NULL != filter_file) {
        FILE *fp = fopen(filter_file, "rb");
        if (!fp) die("Can't open filter file %s\n", filter_file);
        readHeader(fp, &filter_header);
        readFilterData(fp, &filter_header);

        /* set the number of regions by reading the ImageSize.dat file */
        setRegions(&settings);
        if (0 == settings.nregions) {
                die("ERROR: invalid tile size\n");
        }
    }

    /* Look for CIF directories */
    get_cif_dirs(&settings);

    /* open the bam file */
    bam_file = argv[i++];
    fp_bam = samopen(bam_file, "rb", 0);
    if (NULL == fp_bam) {
        fprintf(stderr, "ERROR: can't open bam file file %s: %s\n",
                bam_file, strerror(errno));
        exit(EXIT_FAILURE);
    }

    /* make the survival table */
    nst = makeSurvTable(&settings, fp_bam, sts, &ntiles, &nreads);
    if (0 == nst) {
        fprintf(stderr,"ERROR: failed to make survival table\n");
        exit(EXIT_FAILURE);
    }

    if (!settings.quiet) {
        fprintf(stderr, "Processed %8d traces\n", nreads);
        if (NULL != settings.snp_hash) {
            size_t nsnps = 0;
            int ibucket;
            for (ibucket=0; ibucket<settings.snp_hash->nbuckets; ibucket++) {
                HashItem *hi;
                for (hi = settings.snp_hash->bucket[ibucket]; hi; hi = hi->next)
                    if (hi->data.i)
                        nsnps += hi->data.i;
            }
            fprintf(stderr, "Ignored %lu snps\n", nsnps);
        }
    }

    /* back to where we belong */
    checked_chdir(settings.working_dir);

    if (0 == settings.nregions) makeGlobalSurvTable(&settings, ntiles, sts);

    outputSurvTable(&settings, sts);

    nct = makeCalTable(&settings, sts, cts);
    if (0 == nct) {
        fprintf(stderr,"ERROR: failed to make calibration table\n");
        exit(EXIT_FAILURE);
    }

    outputCalTable(&settings, cts);

    /* close the bam file */
    samclose(fp_bam);

    freeCalTable(&settings, cts);

    freeSurvTable(&settings, sts);

    if (NULL != settings.cif_dirs) free(settings.cif_dirs);

    if (NULL != settings.working_dir) free(settings.working_dir);

    return EXIT_SUCCESS;

}
