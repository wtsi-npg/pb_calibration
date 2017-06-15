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
 * This code generates a purity or quality based calibration table
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
 *
 * PPPRCN
 *   we build an error profile based on a 4-base sequence context
 *   previous, ref, called and next (PRCN) this option will enable
 *   a 7-base context 4 previous, ref, called and next (PPPPRCN)
 */

#define QC_FAIL
#define PROPERLY_PAIRED
//#define CALDATA
//#define CHECK_BASECALL
//#define PPPPRCN


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
#include <search.h>

/* To turn off assert for a small speed gain, uncomment this line */
/* #define NDEBUG */

/* Hack to stop io_lib from trying to include its own config.h */
#ifdef HAVE_CONFIG_H
#undef HAVE_CONFIG_H
#endif
#include <io_lib/misc.h>
#include <io_lib/hash_table.h>

#include <sam.h>
/* BAM_FSUPPLIMENTARY is new to bwa mem and may not be defined in sam.h */
#ifndef BAM_FSUPPLIMENTARY
#define BAM_FSUPPLIMENTARY 2048
#endif

#include <smalloc.h>
#include <aprintf.h>
#include <die.h>
#include <rts.h>
#include <shared.h>
#include <snp.h>
#include <parse_bam.h>
#include <cif.h>

#include <version.h>

/* if we split data by state, using a filter file rather than tile, use tile as a place holder for state */
#define N_STATES 2

#define ST_MODE_PURITY  (1<<0)
#define ST_MODE_QUALITY (1<<1)

#define LEN_SUBST       2
#define NUM_SUBST       16    // 4 ^ LEN_SUBST
#ifdef PPPPRCN
#define LEN_CNTXT       7
#define NUM_CNTXT       16384 // 4 ^ LEN_CNTXT
#else
#define LEN_CNTXT       4
#define NUM_CNTXT       256   // 4 ^ LEN_CNTXT
#endif

#define ST_HILO_PURITY   0.63
#define ST_HILO_QUALITY  29.5

#define ST_STATUS_GOOD  (1<<0)
#define ST_STATUS_BAD   (1<<1)

typedef struct {
    int         mode;
    int         tile;
    int         read;
    int         cycle;
    int         nbins;
    float       offset;
    float       delta;
    float       scale;
    float       predictor_hilo;
    float       *predictor;
    long        *num_bases;
    long        *num_errors;
    long        *subst[NUM_SUBST];
    long        substH[NUM_SUBST];
    long        substL[NUM_SUBST];
    long        cntxtH[NUM_CNTXT];
    long        cntxtL[NUM_CNTXT];
    long        total_bases;
    long        total_errors;
    float       quality;
    int         status;
} SurvTable;

typedef struct {
    int         mode;
    int         tile;
    int         read;
    int         cycle;
    int         nbins;
    float       *predictor;
    long        *num_bases;
    long        *num_errors;
    float       *frac_bases;
    float       *error_rate;
    float       *quality;
} CalTable;

typedef struct {
    char *cmdline;
    char *prefix;
    char *intensity_dir;
    char *snp_file;
    HashTable *snp_hash;
    char *working_dir;
    int read_length[3];
    int cstart[3];
    int quiet;
    int filter_bad_tiles;
    int n_bins_left;
    int n_bins_right;
} Settings;

static void initialiseSurvTable(Settings *s, SurvTable *st, int tile, int read, int cycle)
{
    int i, j;

    st->mode = (NULL == s->intensity_dir ? ST_MODE_QUALITY : ST_MODE_PURITY);

    st->tile  = tile;
    st->read  = read;
    st->cycle = cycle;

    if (st->mode == ST_MODE_PURITY) {
        st->nbins = 76;

        st->offset = 0.25;
        st->delta  = 0.01;
        st->scale  = 1.0 / st->delta;

        st->predictor_hilo = ST_HILO_PURITY;
    } else {
        st->nbins = 51;

        st->offset = 0.0;
        st->delta  = 1.0;
        st->scale  = 1.0 / st->delta;

        st->predictor_hilo = ST_HILO_QUALITY;
    }

    st->predictor = (float *)smalloc(st->nbins * sizeof(float));
    for (i=0;i<st->nbins;i++)
        st->predictor[i] = st->offset + i * st->delta;

    st->num_bases  = (long *)smalloc(st->nbins * sizeof(long));
    st->num_errors = (long *)smalloc(st->nbins * sizeof(long));
    for (i=0;i<st->nbins;i++) {
        st->num_bases[i] = 0;
        st->num_errors[i] = 0;
    }

    for (j=0;j<NUM_SUBST;j++) {
        st->subst[j] = (long *)smalloc(st->nbins * sizeof(long));
        for (i=0;i<st->nbins;i++) 
            st->subst[j][i] = 0;
        st->substH[j]=0;
        st->substL[j]=0;
    }
    for (j=0;j<NUM_CNTXT;j++) {
        st->cntxtH[j]=0;
   	    st->cntxtL[j]=0;
    }

    st->total_bases = 0;
    st->total_errors = 0;

    st->status = ST_STATUS_GOOD;
}

static void freeSurvTable(Settings *s, int ntiles, SurvTable **sts, int no_cycles)
{
    int itile, read, cycle, i;
    for(itile=0;itile<=ntiles;itile++)
    {
        for(read=0;read<N_READS;read++)
        {
            if( NULL == sts[itile*N_READS+read]) continue;
            for(cycle=0;cycle<(no_cycles ? 1 : s->read_length[read]);cycle++)
            {
                SurvTable *st = sts[itile*N_READS+read] + cycle;
                if (st->nbins) {
                    free(st->predictor);
                    free(st->num_bases);
                    free(st->num_errors);
                    for (i=0;i<NUM_SUBST;i++)
                        free(st->subst[i]);
                    st->nbins = 0;
                }
            }
            free(sts[itile*N_READS+read]);
        }
    }
    free(sts);
}

static int optimalPredictorBin(SurvTable *st, int nbins)
{
    int ipopt = 0;
    float max_diff = 0.0;
    long total_bases = 0;
    long total_errors = 0;
    int i;
    long cum_bases, cum_errors;

    /* find optimal predictor which gives the best separation between bases and errors */

    for(i=0;i<nbins;i++)
    {
        total_bases += st->num_bases[i];
        total_errors += st->num_errors[i];
    }
    cum_bases = total_bases;
    cum_errors = total_errors;

    for(i=0;i<nbins;i++)
    {
        float frac_bases, frac_errors, diff_frac;
        frac_bases = (float)cum_bases/(float)total_bases;
        frac_errors = (float)cum_errors/(float)total_errors;
        diff_frac = (frac_bases - frac_errors);
        if (diff_frac > max_diff) {
            max_diff = diff_frac;
            ipopt = i;
        }
        cum_bases -= st->num_bases[i];
        cum_errors -= st->num_errors[i];
    }

    return ipopt;
}

static int maximumQualityBin(SurvTable *st)
{
    float ssc = 1.0;
    int iqmax = 0, ipmax = 0;
    long cum_bases = st->total_bases;
    long cum_errors = st->total_errors;
    int i;

    /* find maximum (integer) quality */
    for(i=0;i<st->nbins;i++)
    {
        float error_rate, quality;
        int iq;
        error_rate = (cum_errors + ssc)/(cum_bases + ssc);
        quality = -10.0 * log10(error_rate);
        iq = (int)(quality + 0.5);
        if (iq > iqmax) {
            iqmax = iq;
            ipmax = i;
        }
        cum_bases -= st->num_bases[i];
        cum_errors -= st->num_errors[i];
    }

    return ipmax;
}

static void completeSurvTable(Settings *s, SurvTable **sts, int no_cycles)
{
    float ssc = 1.0;
    int read, cycle, i;

    for(read=0;read<N_READS;read++)
    {
        if( NULL == sts[read]) continue;
        for(cycle=0;cycle<(no_cycles ? 1 : s->read_length[read]);cycle++)
        {
            SurvTable *st = sts[read] + cycle;
            long quality_bases = 0;
            long quality_errors = 0;

            for(i=0;i<st->nbins;i++)
            {
                st->total_bases += st->num_bases[i];
                st->total_errors += st->num_errors[i];

                // bases in the first bin (purity=0.25 or quality=0) are called as N and explicitly get a quality of 0
                if( i == 0 )
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
    
    float *quality;
    quality = (float *)smalloc(ntiles *sizeof(float));

    for(read=0;read<N_READS;read++)
        for(cycle=0;cycle<s->read_length[read];cycle++)
        {
            float qsum = 0.0;
            float q2sum = 0.0;
            float qmax = -1.0;
            int nq = 0;
            float qavg, q2avg, qstd, qmin;

            /* calculate a phred quality for each tile */
            memset(quality, 0, ntiles * sizeof(float));
            for(itile=0;itile<ntiles;itile++) {
                SurvTable *st = sts[itile*N_READS+read] + cycle;
                quality[itile] = st->quality;
            }

            /* sort by quality and find the tile with the min error rate (highest quality) */
            qsort(quality, ntiles, sizeof(float), float_sort);
            qmax = quality[ntiles-1];
            
            /* skip this cycle if all bases are error, usually this means all bases are called as N */
            if( qmax < 0.0 )
                continue;

            /* calc qmin = qmax-20, we will exclude tiles with an error rate > 100x the min error rate */
            qmin = max(qmax - 20.0, 0.0);

            /* Warning - for tiles with multiple dropouts the aligner can truncate a large fraction of the reads at the same cycle, as 
               the aligner will not usually allow a mismatch at the last aligned base this can artificially reduce the error rate.
               Check the tile with the min error rate is not an outlier by ensuring we retain at least 75% of the tiles */
            qmin = min(qmin, quality[(int)(0.25*ntiles)]);

            for(itile=0;itile<ntiles;itile++) {
                SurvTable *st = sts[itile*N_READS+read] + cycle;
                /* at this stage ALL the tiles should be good */
                assert(st->status == ST_STATUS_GOOD );
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

            if (!s->quiet) {
                display("read=%1d cycle=%-3d qavg=%.2f qstd=%.2f qmax=%.2f qmin=%.2f ", read, cycle, qavg, qstd, qmax, qmin);
                for(itile=0;itile<ntiles;itile++) {
                    SurvTable *st = sts[itile*N_READS+read] + cycle;
                    if (st->status == ST_STATUS_GOOD )
                        display("%4d  %5.2f ", st->tile, st->quality);
                    else
                        // mark poor or bad tiles
                        display("%4d* %5.2f ", st->tile, st->quality);
                }
                display("\n");
            }
        }

        free(quality);
        return;
}

static SurvTable **makeGlobalSurvTable(Settings *s, int ntiles, SurvTable **sts)
{
    int read, read_length, cycle, itile, i, j;

    if (0 >= ntiles)
        return sts;

    findBadTiles(s, ntiles, sts);

    sts = srealloc(sts, (ntiles+1) * N_READS * sizeof(SurvTable *));
    for(read=0;read<N_READS;read++)
        sts[ntiles*N_READS+read] = NULL;

    for(read=0;read<N_READS;read++)
    {
        read_length = s->read_length[read];
        if (0 == read_length) continue;
        sts[ntiles*N_READS+read] = (SurvTable *)smalloc(read_length * sizeof(SurvTable));

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
                    for(j=0;j<NUM_SUBST;j++)
                    {
                        for(i=0;i<st->nbins;i++)
                            st->subst[j][i] += tile_st->subst[j][i];
                        st->substH[j] += tile_st->substH[j];
                        st->substL[j] += tile_st->substL[j];
                    }
                    
                    for(j=0;j<NUM_CNTXT;j++)
                    {
                        st->cntxtH[j] += tile_st->cntxtH[j];
                        st->cntxtL[j] += tile_st->cntxtL[j];
                    }
                    tile_st->total_bases = 0;
                    tile_st->total_errors = 0;
                }
            }
        }
    }

    completeSurvTable(s, &sts[ntiles*N_READS], 0);

    return sts;
}
    
static void outputErrorTable(Settings *s, int ntiles, SurvTable **sts)
{
    FILE *fp = NULL;
    SurvTable **read_sts;
    int read, cycle, i, j;
    char p;

    // open error file
    for(read=0;read<N_READS;read++)
    {
        SurvTable *st = sts[ntiles*N_READS+read];
        int filename_sz;
        char *filename;
        if( NULL == sts[ntiles*N_READS+read] ) continue;
        filename_sz = (NULL == s->prefix ? 0 : strlen(s->prefix)) + 100;
        filename = smalloc(filename_sz);
        sprintf(filename, "%s_%s_error.txt", s->prefix, (st->mode == ST_MODE_PURITY ? "purity" : "quality"));
        fp = fopen(filename, "w");
        if( NULL == fp )
        {
            fprintf(stderr, "ERROR: can't open Error file %s: %s\n", filename, strerror(errno));
            exit(EXIT_FAILURE);
        }
        free(filename);
        break;
    }
    if( NULL == fp )
    {
        fprintf(stderr, "ERROR: no error data\n");
        exit(EXIT_FAILURE);
    }

    /* generate read summary tables by summing over cycles */

    read_sts = (SurvTable **)smalloc(N_READS * sizeof(SurvTable *));
    for(read=0;read<N_READS;read++)
        read_sts[read] = NULL;

    for(read=0;read<N_READS;read++)
    {
        int read_length = s->read_length[read];
        if (0 == read_length) continue;

        read_sts[read] = (SurvTable *)smalloc(sizeof(SurvTable));
        SurvTable *read_st = read_sts[read];

        initialiseSurvTable(s, read_st, -1, read, -1);

        for(cycle=0;cycle<read_length;cycle++)
        {
            SurvTable *st = sts[ntiles*N_READS+read] + cycle;
            for(j=0;j<NUM_SUBST;j++)
                for(i=0;i<st->nbins;i++)
                    read_st->subst[j][i] += st->subst[j][i];
            for(j=0;j<NUM_SUBST;j++)
                read_st->substH[j] += st->substH[j];
            for(j=0;j<NUM_SUBST;j++)
                read_st->substL[j] += st->substL[j];
            for(j=0;j<NUM_CNTXT;j++)
                read_st->cntxtH[j] += st->cntxtH[j];
            for(j=0;j<NUM_CNTXT;j++)
                read_st->cntxtL[j] += st->cntxtL[j];
        }
    }

    /* substitution RC */

    fprintf(fp, "# Substitution error table. Use `grep ^SET | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per predictor, columns read, predictor followed by substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        SurvTable *st = read_sts[read];
        if( NULL == st) continue;

        for(i=0;i<st->nbins;i++)
        {
            fprintf(fp, "SET\t%d\t%.2f", read, st->predictor[i]);
            for (j=0;j<NUM_SUBST;j++)
            {
                char *subst;
                subst=word2str(j,LEN_SUBST);
                if (subst[0]==subst[1]) continue;
                fprintf(fp, "\t%s\t%ld", subst, st->subst[j][i]);
            }
            fprintf(fp, "\n");
        }
    }
    
    fprintf(fp, "# Mismatch substitutions high predictor. Use `grep ^RCH | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and cycle, columns read, cycle then substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        if( NULL == sts[ntiles*N_READS+read]) continue;
        /* cycle by cycle */
        for(cycle=0;cycle<s->read_length[read];cycle++)
        {
            SurvTable *st = sts[ntiles*N_READS+read] + cycle;
            if( 0 == st->total_bases ) continue;
            fprintf(fp, "RCH\t%d\t%d", read, cycle);
            for (j=0;j<NUM_SUBST;j++)
            {
                char *subst;
                subst=word2str(j,LEN_SUBST);
                if (subst[0]==subst[1]) continue;
                fprintf(fp, "\t%s\t%ld", subst, st->substH[j]);
            }
            fprintf(fp, "\n");
        }
        /* and read summary (cycle = -1) */
        SurvTable *st = read_sts[read];
        fprintf(fp, "RCH\t%d\t%d", read, -1);
        for (j=0;j<NUM_SUBST;j++)
        {
            char *subst;
            subst=word2str(j,LEN_SUBST);
            if (subst[0]==subst[1]) continue;
            fprintf(fp, "\t%s\t%ld", subst, st->substH[j]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "# Mismatch substitutions low predictor. Use `grep ^RCL | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and cycle, columns read, cycle then substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        if( NULL == sts[ntiles*N_READS+read]) continue;
        /* cycle by cycle */
        for(cycle=0;cycle<s->read_length[read];cycle++)
        {
            SurvTable *st = sts[ntiles*N_READS+read] + cycle;
            if( 0 == st->total_bases ) continue;
            fprintf(fp, "RCL\t%d\t%d", read, cycle);
            for (j=0;j<NUM_SUBST;j++)
            {
                char *subst;
                subst=word2str(j,LEN_SUBST);
                if (subst[0]==subst[1]) continue;
                fprintf(fp, "\t%s\t%ld", subst, st->substL[j]);
            }
            fprintf(fp, "\n");
        }
        /* and read summary (cycle = -1) */
        SurvTable *st = read_sts[read];
        fprintf(fp, "RCL\t%d\t%d", read, -1);
        for (j=0;j<NUM_SUBST;j++)
        {
            char *subst;
            subst=word2str(j,LEN_SUBST);
            if (subst[0]==subst[1]) continue;
            fprintf(fp, "\t%s\t%ld", subst, st->substL[j]);
        }
        fprintf(fp, "\n");
    }
    
    /* previous base PRC */

    fprintf(fp, "# Effect of previous base high predictor. Use `grep ^PRCH | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and previous base, columns read then previous base+substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        /* read summary */
        SurvTable *st = read_sts[read];
        if( NULL == st) continue;
#ifdef PPPPRCN
        int cntxt_off = 3;
#else
        int cntxt_off = 0;
#endif
        int len_cntxt = 3;
        int num_cntxt = 64;
        long *count = (long *)smalloc(num_cntxt * sizeof(long));
        for (j=0;j<num_cntxt;j++)
            count[j]=0;
        for (j=0;j<NUM_CNTXT;j++)
        {
            char *cntxt;
            int word;
            cntxt=word2str(j,LEN_CNTXT);
            cntxt+=cntxt_off;
            word=str2word(cntxt,len_cntxt);
            count[word] += st->cntxtH[j];
        }
        for (j=0,p=0;j<num_cntxt;j++)
        {
            char *cntxt;
            cntxt=word2str(j,len_cntxt);
            if (p != cntxt[0]) {
                p = cntxt[0];
                if (j) fprintf(fp, "\n");
                fprintf(fp, "PRCH\t%d", read);
            }
            if (cntxt[1]==cntxt[2]) continue;
            fprintf(fp, "\t%s\t%ld", cntxt, count[j]);
        }
        fprintf(fp, "\n");
        free(count);
    }
    
    fprintf(fp, "# Effect of previous base low predictor. Use `grep ^PRCL | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and previous base, columns read then previous base+substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        /* read summary */
        SurvTable *st = read_sts[read];
        if( NULL == st) continue;
#ifdef PPPPRCN
        int cntxt_off = 3;
#else
        int cntxt_off = 0;
#endif
        int len_cntxt = 3;
        int num_cntxt = 64;
        long *count = (long *)smalloc(num_cntxt * sizeof(long));
        for (j=0;j<num_cntxt;j++)
            count[j]=0;
        for (j=0;j<NUM_CNTXT;j++)
        {
            char *cntxt;
            int word;
            cntxt=word2str(j,LEN_CNTXT);
            cntxt+=cntxt_off;
            word=str2word(cntxt,len_cntxt);
            count[word] += st->cntxtL[j];
        }
        for (j=0,p=0;j<num_cntxt;j++)
        {
            char *cntxt;
            cntxt=word2str(j,len_cntxt);
            if (p != cntxt[0]) {
                p = cntxt[0];
                if (j) fprintf(fp, "\n");
                fprintf(fp, "PRCL\t%d", read);
            }
            if (cntxt[1]==cntxt[2]) continue;
            fprintf(fp, "\t%s\t%ld", cntxt, count[j]);
        }
        fprintf(fp, "\n");
        free(count);
    }
    
    /* previous base + next base PRCN */

    fprintf(fp, "# Effect of previous base and next base high predictor. Use `grep ^PRCNH | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and previous/next base, columns read then previous base+substitution+next base and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        /* read summary */
        SurvTable *st = read_sts[read];
        if( NULL == st) continue;
#ifdef PPPPRCN
        int cntxt_off = 3;
#else
        int cntxt_off = 0;
#endif
        int len_cntxt = 4;
        int num_cntxt = 256;
        long *count = (long *)smalloc(num_cntxt * sizeof(long));
        for (j=0;j<num_cntxt;j++)
            count[j]=0;
        for (j=0;j<NUM_CNTXT;j++)
        {
            char *cntxt;
            int word;
            cntxt=word2str(j,LEN_CNTXT);
            cntxt+=cntxt_off;
            word=str2word(cntxt,len_cntxt);
            count[word] += st->cntxtH[j];
        }
        for (j=0,p=0;j<num_cntxt;j++)
        {
            char *cntxt;
            cntxt=word2str(j,len_cntxt);
            if (p != cntxt[1]) {
                p = cntxt[1];
                if (j) fprintf(fp, "\n");
                fprintf(fp, "PRCNH\t%d", read);
            }
            if (cntxt[1]==cntxt[2]) continue;
            fprintf(fp, "\t%s\t%ld", cntxt, count[j]);
        }
        fprintf(fp, "\n");
        free(count);
    }
    
    fprintf(fp, "# Effect of previous base and next base low predictor. Use `grep ^PRCNL | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and previous/next base, columns read then previous base+substitution+next base and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        /* read summary */
        SurvTable *st = read_sts[read];
        if( NULL == st ) continue;
#ifdef PPPPRCN
        int cntxt_off = 3;
#else
        int cntxt_off = 0;
#endif
        int len_cntxt = 4;
        int num_cntxt = 256;
        long *count = (long *)smalloc(num_cntxt * sizeof(long));
        for (j=0;j<num_cntxt;j++)
            count[j]=0;
        for (j=0;j<NUM_CNTXT;j++)
        {
            char *cntxt;
            int word;
            cntxt=word2str(j,LEN_CNTXT);
            cntxt+=cntxt_off;
            word=str2word(cntxt,len_cntxt);
            count[word] += st->cntxtL[j];
        }
        for (j=0,p=0;j<num_cntxt;j++)
        {
            char *cntxt;
            cntxt=word2str(j,len_cntxt);
            if (p != cntxt[1]) {
                p = cntxt[1];
                if (j) fprintf(fp, "\n");
                fprintf(fp, "PRCNL\t%d", read);
            }
            if (cntxt[1]==cntxt[2]) continue;
            fprintf(fp, "\t%s\t%ld", cntxt, count[j]);
        }
        fprintf(fp, "\n");
        free(count);
    }
    
#ifdef PPPPRCN
    /* previous 4 base homopolymer HRC. */

    fprintf(fp, "# Homopolymer effect high predictor. Use `grep ^HRCH | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and homopolymer, columns read then homopolymer+substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        if( NULL == sts[ntiles*N_READS+read]) continue;
        /* read by read */
        SurvTable *st = read_sts[read];
        int len_homop = 4;
        int cntxt_off = 3;
        int len_cntxt = 3;
        int num_cntxt = 64;
        long *count = (long *)smalloc(num_cntxt * sizeof(long));
        for (j=0;j<num_cntxt;j++)
            count[j]=0;
        for (j=0;j<NUM_CNTXT;j++)
        {
            int k;
            char *cntxt;
            int word;
            cntxt=word2str(j,LEN_CNTXT);
            for (k=0;k<len_homop;k++)
                if (cntxt[k] != cntxt[cntxt_off]) break;
            if (k != len_homop) continue;
            cntxt+=cntxt_off;
            word=str2word(cntxt,len_cntxt);
            count[word] += st->cntxtH[j];
        }
        for (j=0,p=0;j<num_cntxt;j++)
        {
            int k;
            char homop[len_homop];
            char *cntxt;
            cntxt=word2str(j,len_cntxt);
            if (p != cntxt[0]) {
                p = cntxt[0];
                if (j) fprintf(fp, "\n");
                fprintf(fp, "HRCH\t%d", read);
            }
            if (cntxt[1]==cntxt[2]) continue;
            for (k=0;k<(len_homop-1);k++)
                homop[k]=cntxt[0];
            homop[k]=0;
            fprintf(fp, "\t%s%s\t%ld", homop, cntxt, count[j]);
        }
        fprintf(fp, "\n");
        free(count);
    }
    
    fprintf(fp, "# Homopolymer effect low predictor. Use `grep ^HRCL | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and homopolymer, columns read then homopolymer+substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        if( NULL == sts[ntiles*N_READS+read]) continue;
        /* read by read */
        SurvTable *st = read_sts[read];
        int len_homop = 4;
        int cntxt_off = 3;
        int len_cntxt = 3;
        int num_cntxt = 64;
        long *count = (long *)smalloc(num_cntxt * sizeof(long));
        for (j=0;j<num_cntxt;j++)
            count[j]=0;
        for (j=0;j<NUM_CNTXT;j++)
        {
            int k;
            char *cntxt;
            int word;
            cntxt=word2str(j,LEN_CNTXT);
            for (k=0;k<len_homop;k++)
                if (cntxt[k] != cntxt[cntxt_off]) break;
            if (k != len_homop) continue;
            cntxt+=cntxt_off;
            word=str2word(cntxt,len_cntxt);
            count[word] += st->cntxtL[j];
        }
        for (j=0,p=0;j<num_cntxt;j++)
        {
            int k;
            char homop[len_homop];
            char *cntxt;
            cntxt=word2str(j,len_cntxt);
            if (p != cntxt[0]) {
                p = cntxt[0];
                if (j) fprintf(fp, "\n");
                fprintf(fp, "HRCL\t%d", read);
            }
            if (cntxt[1]==cntxt[2]) continue;
            for (k=0;k<(len_homop-1);k++)
                homop[k]=cntxt[0];
            homop[k]=0;
            fprintf(fp, "\t%s%s\t%ld", homop, cntxt, count[j]);
        }
        fprintf(fp, "\n");
        free(count);
    }
    
    /* previous 4 bases + next base PPPPRCN - only the 10 most common */

    fprintf(fp, "# Effect of previous 4 bases and next base high predictor. Use `grep ^PPPPRCNH | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and previous/next base, columns read then previous base+substitution+next base and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        /* read summary */
        SurvTable *st = read_sts[read];
        if( NULL == st ) continue;
	    long *count = (long *)smalloc(NUM_CNTXT * sizeof(long));
        for (j=0;j<NUM_CNTXT;j++)
	        count[j] = st->cntxtH[j];
        qsort(count, NUM_CNTXT, sizeof(long), long_sort);
        int threshold = (count[NUM_CNTXT-10] > 0 ? count[NUM_CNTXT-10] : 1);
        fprintf(fp, "PPPPRCNH\t%d", read);
        for (j=0;j<NUM_CNTXT;j++)
        {
            char *cntxt;
            cntxt=word2str(j,LEN_CNTXT);
            if (st->cntxtH[j] < threshold) continue;
            fprintf(fp, "\t%s\t%ld", cntxt, st->cntxtH[j]);
        }
        fprintf(fp, "\n");
        free(count);
    }

    fprintf(fp, "# Effect of previous 4 bases and next base low predictor. Use `grep ^PPPPRCNL | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and previous/next base, columns read then previous base+substitution+next base and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        /* read summary */
        SurvTable *st = read_sts[read];
        if( NULL == st ) continue;
	    long *count = (long *)smalloc(NUM_CNTXT * sizeof(long));
        for (j=0;j<NUM_CNTXT;j++)
	        count[j] = st->cntxtL[j];
        qsort(count, NUM_CNTXT, sizeof(long), long_sort);
        int threshold = (count[NUM_CNTXT-12] > 0 ? count[NUM_CNTXT-12] : 1);
        fprintf(fp, "PPPPRCNL\t%d", read);
        for (j=0;j<NUM_CNTXT;j++)
        {
            char *cntxt;
            cntxt=word2str(j,LEN_CNTXT);
            if (st->cntxtL[j] < threshold) continue;
            fprintf(fp, "\t%s\t%ld", cntxt, st->cntxtL[j]);
        }
        fprintf(fp, "\n");
        free(count);
    }
#endif

    freeSurvTable(s, 0, read_sts, 1);
    if( NULL != fp) fclose(fp);
}

static void outputSurvTable(Settings *s, int ntiles, SurvTable **sts)
{
    FILE *fp = NULL;
    int itile, read, cycle, i;

    for(itile=0;itile<=ntiles;itile++)
        for(read=0;read<N_READS;read++)
        {
            if( NULL == sts[itile*N_READS+read]) continue;
            for(cycle=0;cycle<s->read_length[read];cycle++)
            {
                SurvTable *st = sts[itile*N_READS+read] + cycle;

                // skip st with no data
                if( 0 == st->total_bases ) continue;

                if( NULL == fp )
                {
                    int filename_sz;
                    char *filename;
                    filename_sz = (NULL == s->prefix ? 0 : strlen(s->prefix)) + 100;
                    filename = smalloc(filename_sz);
                    sprintf(filename, "%s_%s_cycle_surv.txt", s->prefix, (st->mode == ST_MODE_PURITY ? "purity" : "quality"));
                    fp = fopen(filename, "w");
                    if( NULL == fp )
                    {
                        fprintf(stderr, "ERROR: can't open survival table file %s: %s\n",
                                filename, strerror(errno));
                        exit(EXIT_FAILURE);
                    }
                    free(filename);
                }

                for(i=0;i<st->nbins;i++)
                    fprintf(fp, "%.2f\t%d\t%d\t%d\t%ld\t%ld\n",
                            st->predictor[i], st->read, st->cycle, st->tile,
                            st->num_bases[i], st->num_errors[i]);
            }
        }

    if( NULL != fp) fclose(fp);
}

static void initialiseCalTable(SurvTable *st, CalTable *ct)
{
    int i;

    ct->mode  = st->mode;

    ct->tile  = st->tile;
    ct->read  = st->read;
    ct->cycle = st->cycle;

    ct->nbins = st->nbins;

    ct->predictor = (float *)smalloc(ct->nbins * sizeof(float));
    for(i=0;i<ct->nbins;i++)
        ct->predictor[i] = st->predictor[i];

    ct->num_bases  = (long *)smalloc(ct->nbins * sizeof(long));
    ct->num_errors = (long *)smalloc(ct->nbins * sizeof(long));
    ct->frac_bases = (float *)smalloc(ct->nbins * sizeof(float));
    ct->error_rate = (float *)smalloc(ct->nbins * sizeof(float));
    ct->quality    = (float *)smalloc(ct->nbins * sizeof(float));
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

static void freeCalTable(Settings *s, int ntiles, CalTable **cts)
{
    int itile, read, cycle;
    for(itile=0;itile<=ntiles;itile++)
    {
        for(read=0;read<N_READS;read++)
        {
            if( NULL == cts[itile*N_READS+read]) continue;
            for(cycle=0;cycle<s->read_length[read];cycle++)
            {
                CalTable *ct = cts[itile*N_READS+read] + cycle;
                if(ct->nbins) {
                    free(ct->predictor);
                    free(ct->num_bases);
                    free(ct->num_errors);
                    free(ct->frac_bases);
                    free(ct->error_rate);
                    free(ct->quality);

                    ct->nbins = 0;
                }
            }
            free(cts[itile*N_READS+read]);
        }
    }
    free(cts);
}

static void optimisePredictorBins(Settings *s, SurvTable *st, CalTable *ct)
{
    float predictor[st->nbins];
    int nbins = 0;

    int i, j;
    int ipopt = 0, ipmax = 0;
    int ipbin;

    ipopt = optimalPredictorBin(st, st->nbins);
    ipmax = maximumQualityBin(st);

    /* if ipmax < ipopt reset ipmax to ipopt */
    if(ipmax < ipopt)
    {
        ipmax = ipopt;
        fprintf(stderr, "Resetting maximal bin to optimal bin read=%d cycle=%d tile=%d\n", st->read, st->cycle, st->tile);
    }

    /* first bin */
    ipbin = 0;
    predictor[nbins++] = st->predictor[ipbin];
#if 0
    display("1: nbins=%d ipbin=%d predictor=%f\n", nbins, ipbin, st->predictor[ipbin]);
#endif

    if (ipopt > 0)
    {
#ifdef POPT2
        /* find the 2nd optimal bin */
        int ipopt2 = optimalPredictorBin(st, ipopt);
        ipbin = (ipopt2 > 2 ? (ipopt2-1) : 1);
        /* add one bin before popt2 and ALL other bins up to popt */
        for(;ipbin<=ipopt;ipbin++)
        {
            predictor[nbins++] = st->predictor[ipbin];
#if 0
            display("2: nbins=%d ipbin=%d predictor=%f\n", nbins, ipbin, st->predictor[ipbin]);
#endif
        }

#else // POPT2
        /* add n_bins_left bins between first bin and popt */
        float pinc = (st->predictor[ipopt] - st->predictor[0]) / (s->n_bins_left + 1);
        for(i=0,j=1;i<s->n_bins_left;i++,j++)
        {
            float p = st->predictor[0] + j * pinc;
            int ip = st->scale * (p - st->offset) + 0.5;
            if (ip == ipopt)
                break;
            if (ip > ipbin) {
                ipbin = ip;
                predictor[nbins++] = st->predictor[ipbin];
#if 0
                display("2: nbins=%d ipbin=%d predictor=%f\n", nbins, ipbin, st->predictor[ipbin]);
#endif
            }
        }

        /* add the optimal bin */
        ipbin = ipopt;
        predictor[nbins++] = st->predictor[ipbin];
#if 0
        display("3: nbins=%d ipbin=%d predictor=%f\n", nbins, ipbin, st->predictor[ipbin]);
#endif
#endif // POPT2
    }

    if (ipmax > ipopt)
    {
        /* add n_bins_right bins between popt and pmax */
        float pinc = (st->predictor[ipmax] - st->predictor[ipopt]) / (s->n_bins_right + 1);
        for(i=0,j=1;i<s->n_bins_right;i++,j++)
        {
            float p = st->predictor[ipopt] + j * pinc;
            int ip = st->scale * (p - st->offset) + 0.5;
            if (ip == ipmax)
                break;
            if (ip > ipbin) {
                ipbin = ip;
                predictor[nbins++] = st->predictor[ipbin];
#if 0
                display("4: nbins=%d ipbin=%d predictor=%f\n", nbins, ipbin, st->predictor[ipbin]);
#endif
            }
        }

        /* add pmax */
        ipbin = ipmax;
        predictor[nbins++] = st->predictor[ipbin];
#if 0
        display("5: nbins=%d ipbin=%d predictor=%f\n", nbins, ipbin, st->predictor[ipbin]);
#endif
    }

    /* add another bin if pmax is not the maximum possible predictor value */
    if (ipmax < (st->nbins-1))
        predictor[nbins++] = st->predictor[st->nbins-1];
#if 0
    if (ipmax < (st->nbins-1))
        display("6: nbins=%d ipbin=%d predictor=%f\n", nbins, st->nbins-1, st->predictor[st->nbins-1]);
#endif

    /* move data into the new set of bins */

    ct->nbins = nbins;

    for(i=0;i<ct->nbins;i++)
    {
        ct->predictor[i] = predictor[i];
        ct->num_bases[i]  = 0;
        ct->num_errors[i] = 0;
    }

    for(i=0;i<st->nbins;i++)
    {
        for(j=0;j<ct->nbins;j++)
        {
            if(st->predictor[i] <= ct->predictor[j])
            {
                ct->num_bases[j]  += st->num_bases[i];
                ct->num_errors[j] += st->num_errors[i];
                break;
            }
        }
    }
}

static CalTable **makeCalTable(Settings *s, int ntiles, SurvTable **sts)
{
    CalTable **cts = NULL;
    float ssc = 1.0;
    int itile, read, read_length, cycle, i;

    cts = (CalTable **)smalloc((ntiles+1) * N_READS * sizeof(CalTable *));

    for(itile=0;itile<=ntiles;itile++)
        for(read=0;read<N_READS;read++)
        {
            cts[itile*N_READS+read] = NULL;
            if( NULL == sts[itile*N_READS+read]) continue;

            read_length = s->read_length[read];
            cts[itile*N_READS+read] = (CalTable *)smalloc(read_length * sizeof(CalTable));
            
            for(cycle=0;cycle<read_length;cycle++)
            {
                SurvTable *st = sts[itile*N_READS+read] + cycle;
                CalTable *ct  = cts[itile*N_READS+read] + cycle;

                // set number of bins in ct to 0
                ct->nbins = 0;

                // skip st with no data
                if( 0 == st->total_bases ) continue;

                initialiseCalTable(st, ct);

                optimisePredictorBins(s, st, ct);
                
                for(i=0;i<ct->nbins;i++)
                {
                    ct->frac_bases[i] = ((float)ct->num_bases[i]) / ((float)st->total_bases);
                    ct->error_rate[i] = (ct->num_errors[i] + ssc)/(ct->num_bases[i] + ssc);
                    ct->quality[i] = -10.0 * (log10(ct->error_rate[i]));
                }
            }
        }

    return cts;
}

static void outputCalTable(Settings *s, int ntiles, CalTable **cts)
{
    FILE *fp = NULL;
    int itile, read, cycle, i;

    for(itile=0;itile<=ntiles;itile++)
        for(read=0;read<N_READS;read++)
        {
            if( NULL == cts[itile*N_READS+read]) continue;
            for(cycle=0;cycle<s->read_length[read];cycle++)
            {
                CalTable *ct = cts[itile*N_READS+read] + cycle;

                // skip ct with no bins
                if( 0 == ct->nbins ) continue;

                if( NULL == fp )
                {
                    int filename_sz;
                    char *filename;
                    filename_sz = (NULL == s->prefix ? 0 : strlen(s->prefix)) + 100;
                    filename = smalloc(filename_sz);
                    sprintf(filename, "%s_%s_cycle_caltable.txt", s->prefix, (ct->mode == ST_MODE_PURITY ? "purity" : "quality"));
                    fp = fopen(filename, "w");
                    if( NULL == fp )
                    {
                        fprintf(stderr, "ERROR: can't open CT file %s: %s\n",
                                filename, strerror(errno));
                        exit(EXIT_FAILURE);
                    }
                    free(filename);
                }

                for(i=0;i<ct->nbins;i++)
                {
                    /* always output the first and last bin but exclude all other empty bins */
                    if( ct->num_bases[i] == 0 && i > 0 && i < (ct->nbins-1) )
                        continue;

                    fprintf(fp, "%f\t%d\t%d\t%d\t%5ld\t%5ld\t%f\t%f\t%f\n",
                            ct->predictor[i], ct->read, ct->cycle, ct->tile,
                            ct->num_bases[i], ct->num_errors[i], ct->frac_bases[i],
                            ct->error_rate[i], ct->quality[i]);
                }
            }
        }

    if( NULL != fp) fclose(fp);
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


static int updateSurvTable(Settings *s, SurvTable **sts, CifData *cif_data,
                           size_t spot_num, int tile, int x, int y, int read, int *read_mismatch,
                           char *read_seq, int *read_qual, char *read_ref, samfile_t *fp, bam1_t *bam, FILE *fp_caldata) {

    int cstart = s->cstart[read];
    int read_length = s->read_length[read];
    int c, b;

    if (NULL != cif_data) assert((cstart + read_length) <= cif_data->ncycles);

#ifdef CALDATA
    int dir = BAM_FREVERSE & bam->core.flag ? 1 : 0;
    char *chrom = fp->header->target_name[bam->core.tid];
    int pos = bam->core.pos;

    fprintf(fp_caldata, "%d\t%d\t%d\t%d\t%s\t%d\t", read, x, y, dir, chrom, pos);
#endif

    /* update survival table */
    for (c = cstart, b = 0; b < read_length; c++, b++) {
        SurvTable *st = sts[read] + b;
        float predictor = -1.0;
        int ibin;

#ifdef CALDATA
        char bases[] = "ACGTNDIS";
        char *base = strchr(bases, read_seq[b]);
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
#endif

        if (st->mode == ST_MODE_PURITY) {
            CifCycleData *cycle = cif_data->cycles + c;
            int channel;
            int bin[cycle->num_channels];

            read_cif_chunk(cycle, spot_num);
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

            predictor = GetPu(cycle->num_channels, bin);

#ifdef CALDATA
            fprintf(fp_caldata, "%d %d %d %d %d %d %d %d\t",
                    ref_base, seq_base, read_qual[b], bin[0], bin[1], bin[2], bin[3], read_mismatch[b]);
#endif

#ifdef CHECK_BASECALL
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
                            bin[0], flags[0], bin[1], flags[1], bin[2], flags[2], bin[3], flags[3], predictor);
            }
#endif
        } else {
            predictor = read_qual[b];

#ifdef CALDATA
            fprintf(fp_caldata, "%d %d %d %d\t",
                    ref_base, seq_base, read_qual[b], read_mismatch[b]);
#endif
        }
        
        ibin = st->scale * (predictor - st->offset) + 0.5;
        if( ibin >= st->nbins ) ibin=(st->nbins-1);

        if( read_mismatch[b] & BASE_KNOWN_SNP ) {
            // don't count these
        } else {
            if( read_mismatch[b] & BASE_ALIGN )
                st->num_bases[ibin]++;
            if( read_mismatch[b] & BASE_MISMATCH )
                st->num_errors[ibin]++;
            if( read_mismatch[b] & BASE_MISMATCH ){
                char subst[LEN_SUBST+1];
                char cntxt[LEN_CNTXT+1];
                int chr, word;

                chr=0;
                subst[chr++]=read_ref[b];
                subst[chr++]=read_seq[b];
                word=str2word(subst, LEN_SUBST);
                if( word >= 0 ) {
                    st->subst[word][ibin]++;
                    if( predictor >= st->predictor_hilo )
                        st->substH[word]++;
                    else
                        st->substL[word]++;
                }

                chr=0;
#ifdef PPPPRCN
                cntxt[chr++]=(b > 3 ? read_ref[b-4] : 'N');
                cntxt[chr++]=(b > 2 ? read_ref[b-3] : 'N');
                cntxt[chr++]=(b > 1 ? read_ref[b-2] : 'N');
#endif
                cntxt[chr++]=(b > 0 ? read_ref[b-1] : 'N');
                cntxt[chr++]=read_ref[b];
                cntxt[chr++]=read_seq[b];
                cntxt[chr++]=(b < (read_length-1) ? read_ref[b+1] : 'N');
                word=str2word(cntxt, LEN_CNTXT);
                if( word >= 0 ) {
                    if( predictor >= st->predictor_hilo )
                        st->cntxtH[word]++;
                    else
                        st->cntxtL[word]++;
                }
            }
        }
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
 * Returns: +'ve integer for success
 *          0 for failure
 */
SurvTable **makeSurvTable(Settings *s, samfile_t *fp_bam, int *bam_ntiles, size_t *bam_nreads) {
    SurvTable **sts = NULL;

    int *tiles = NULL;
    int ntiles = 0;
    int itile = -1;

    size_t nreads = 0;

    CifData *cif_data = NULL;
    int ncycles_firecrest = -1;

    int lane = -1;
    int tile = -1;

    static const int bam_read_buff_size = 1024;
    char bam_read_seq[bam_read_buff_size];
    int bam_read_qual[bam_read_buff_size];
    char bam_read_ref[bam_read_buff_size];
    int bam_read_mismatch[bam_read_buff_size];

    FILE *fp_caldata = NULL;

    bam1_t *bam = bam_init1();

    if (NULL != s->intensity_dir) checked_chdir(s->intensity_dir);

    /* loop over reads in the bam file */
    while (1){
        int bam_lane = -1, bam_tile = -1, bam_x = -1, bam_y = -1, bam_read = -1, read_length;
        size_t bam_offset = 0;

        if (parse_bam_readinfo(fp_bam, bam, &bam_lane, &bam_tile, &bam_x, &bam_y, &bam_read, (NULL == s->intensity_dir ? NULL : &bam_offset))) {
            break;	/* break on end of BAM file */
	}

        if (BAM_FUNMAP & bam->core.flag) continue;
        if (BAM_FQCFAIL & bam->core.flag) continue;
        if (BAM_FSECONDARY & bam->core.flag) continue;
        if (BAM_FSUPPLIMENTARY & bam->core.flag) continue;
        if (BAM_FPAIRED & bam->core.flag) {
            if (BAM_FMUNMAP & bam->core.flag) continue;
            if (0 == (BAM_FPROPER_PAIR & bam->core.flag)) {
                continue;
            }
        }

        read_length = bam->core.l_qseq;
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

        parse_bam_alignments(fp_bam, bam, bam_read_seq, bam_read_qual, bam_read_ref,
                                     bam_read_mismatch, bam_read_buff_size,s->snp_hash);

        if( bam_tile != tile ){
   	        size_t nelem = ntiles;
	        void *pitile;

            tile = bam_tile;
            
            // lookup itile from tile in tile array
	        pitile = lfind(&tile, tiles, &nelem, sizeof(int), &int_cmp);
	        if (NULL == pitile) {
                int read;
                itile = ntiles;
                ntiles++;
                tiles = srealloc(tiles, ntiles * sizeof(int));
                tiles[itile] = tile;
                sts = srealloc(sts, ntiles * N_READS * sizeof(SurvTable *));
                for(read=0;read<N_READS;read++)
                    sts[itile*N_READS+read] = NULL;
                if (!s->quiet) fprintf(stderr, "Processing tile %i (%lu)\n", tile, nreads);
            }else{
                if( NULL != s->intensity_dir ) {
                    fprintf(stderr,"ERROR: alignments are not sorted by tile.\n");
                    exit(EXIT_FAILURE);
                }
                itile = ((int*)pitile - tiles);
            }
            
            if ( NULL != s->intensity_dir) {
                /* Look for processed trace data */
                if (NULL != cif_data) free_cif_data(cif_data);
                cif_data = load_cif_data(lane, tile, "dif");
            
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
            }

#ifdef CALDATA
            int filename_sz = strlen(s->working_dir) + (NULL == s->prefix ? 0 : strlen(s->prefix)) + 100;
            char *filename;
            filename = smalloc(filename_sz);
            sprintf(filename, "%s/%s_%04d_caldata.txt", s->working_dir, s->prefix, bam_tile);
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

        if (NULL == sts[itile*N_READS+bam_read]) {
            int cycle;
            sts[itile*N_READS+bam_read] = (SurvTable *)smalloc(read_length * sizeof(SurvTable));
            for(cycle=0;cycle<read_length;cycle++)
                initialiseSurvTable(s, sts[itile*N_READS+bam_read]+cycle, tile, bam_read, cycle);
        }
        
        if (0 != updateSurvTable(s, &sts[itile*N_READS], cif_data,
                                 bam_offset, bam_tile, bam_x, bam_y, bam_read, bam_read_mismatch,
                                 bam_read_seq, bam_read_qual, bam_read_ref, fp_bam, bam, fp_caldata)) {
            fprintf(stderr,"ERROR: updating quality values for tile %i.\n", tile);
            exit(EXIT_FAILURE);
        }
        
        nreads++;
    }
    
    for(itile=0;itile<ntiles;itile++)
        completeSurvTable(s, &sts[itile*N_READS], 0);

    if (NULL != tiles) free(tiles);
    if (NULL != cif_data) free_cif_data(cif_data);
    if (NULL != fp_caldata) fclose(fp_caldata);

    bam_destroy1(bam);
    
    *bam_ntiles = ntiles;
    *bam_nreads = nreads;

    return sts;
}

static
void usage(int code) {
    FILE* usagefp = stderr;

    fprintf(usagefp, "pb_calibration %s\n\n", version);
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
    fprintf(usagefp, "               number of bins between minimum predictor and optimal predictor\n");
    fprintf(usagefp, "                 default 2\n");
    fprintf(usagefp, "    -nR nbins\n");
    fprintf(usagefp, "               number of bins between optimal predictor and qmax\n");
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


int main(int argc, char **argv) {

    Settings settings;
    int i;
    char *bam_file = NULL;
    samfile_t *fp_bam = NULL;
    const char *override_intensity_dir = NULL;
    const char *filter_file = NULL;
    int ntiles = 0;
    size_t nreads = 0;
    SurvTable **sts = NULL;
    CalTable **cts = NULL;

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
        fprintf(stderr,"Building purity cycle calibration table\n");
    } else {
        fprintf(stderr,"Building quality cycle calibration table\n");
    }

    /* read the snp_file */
    if (NULL != settings.snp_file) {
        settings.snp_hash = readSnpFile(settings.snp_file);
        if (NULL == settings.snp_hash) {
            fprintf(stderr, "ERROR: reading snp file %s\n", settings.snp_file);
            exit(EXIT_FAILURE);
        }
    }

    /* Look for CIF directories */
    if (NULL != settings.intensity_dir) {
        get_cif_dirs(settings.intensity_dir);
    }

    /* open the bam file */
    bam_file = argv[i++];
    fp_bam = samopen(bam_file, "rb", 0);
    if (NULL == fp_bam) {
        fprintf(stderr, "ERROR: can't open bam file file %s: %s\n",
                bam_file, strerror(errno));
        exit(EXIT_FAILURE);
    }

    /* make the survival table */
    sts = makeSurvTable(&settings, fp_bam, &ntiles, &nreads);
    if (0 == nreads) {
        fprintf(stderr,"WARNING: No data in BAM file\n");
        exit(EXIT_SUCCESS);
    }

    if (!settings.quiet) {
        fprintf(stderr, "Processed %8lu traces\n", nreads);
        if (NULL != settings.snp_hash) {
            HashIter *iter = HashTableIterCreate();
            HashItem *hashItem;
            size_t nsnps = 0;
            while ((hashItem = HashTableIterNext(settings.snp_hash, iter)))
                nsnps += hashItem->data.i;
            fprintf(stderr, "Ignored %lu snps\n", nsnps);
            HashTableIterDestroy(iter);
        }
    }

    /* back to where we belong */
    checked_chdir(settings.working_dir);

    sts = makeGlobalSurvTable(&settings, ntiles, sts);

    outputSurvTable(&settings, ntiles, sts);

    outputErrorTable(&settings, ntiles, sts);

    cts = makeCalTable(&settings, ntiles, sts);
    if (NULL == cts) {
        fprintf(stderr,"ERROR: failed to make calibration table\n");
        exit(EXIT_FAILURE);
    }

    outputCalTable(&settings, ntiles, cts);

    /* close the bam file */
    samclose(fp_bam);

    freeCalTable(&settings, ntiles, cts);
    
    freeSurvTable(&settings, ntiles, sts, 0);

    if (NULL != settings.snp_hash) HashTableDestroy(settings.snp_hash, 0);
    if (NULL != settings.working_dir) free(settings.working_dir);
    if (NULL != settings.cmdline) free(settings.cmdline);
    if (NULL != settings.intensity_dir) free(settings.intensity_dir);

    return EXIT_SUCCESS;

}
