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
 * ST_STRUCTURE
 *   generate error stats
 */

#define QC_FAIL
#define PROPERLY_PAIRED
//#define CALDATA
//#define CHECK_BASECALL
#define ST_STRUCTURE

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
    float       *predictor;
    long        *num_bases;
    long        *num_errors;
#ifdef ST_STRUCTURE
    long        *subst[NUM_SUBST];
    long        *cntxt[NUM_CNTXT];
    long        *homop[NUM_CNTXT];
    long        substH[NUM_SUBST];
    long        substL[NUM_SUBST];
    long        cntxtH[NUM_CNTXT];
    long        cntxtL[NUM_CNTXT];
    long        homopH[NUM_CNTXT];
    long        homopL[NUM_CNTXT];
#endif    
    long        total_bases;
    long        total_errors;
    float       optimal_predictor;
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
    int spatial_filter;
    int region_size;
    int nregions_x;
    int nregions_y;
    int n_bins_left;
    int n_bins_right;
} Settings;

static void initialiseSurvTable(Settings *s, SurvTable *st, int tile, int read, int cycle)
{
    int i;

    st->mode = (NULL == s->intensity_dir ? ST_MODE_QUALITY : ST_MODE_PURITY);

    st->tile  = tile;
    st->read  = read;
    st->cycle = cycle;

    if (st->mode == ST_MODE_PURITY) {
        st->nbins = 76;

        st->offset = 0.25;
        st->delta  = 0.01;
        st->scale  = 1.0 / st->delta;
    } else {
        st->nbins = 51;

        st->offset = 0.0;
        st->delta  = 1.0;
        st->scale  = 1.0 / st->delta;
    }

    st->predictor = smalloc(st->nbins * sizeof(float));
    for (i=0;i<st->nbins;i++)
        st->predictor[i] = st->offset + i * st->delta;

    st->num_bases  = smalloc(st->nbins * sizeof(long));
    st->num_errors = smalloc(st->nbins * sizeof(long));
    for (i=0;i<st->nbins;i++) {
        st->num_bases[i] = 0;
        st->num_errors[i] = 0;
    }
#ifdef ST_STRUCTURE
    int j;
    for (j=0;j<NUM_SUBST;j++) {
        st->subst[j] = (long *)smalloc(st->nbins * sizeof(long));
        for (i=0;i<st->nbins;i++) 
            st->subst[j][i] = 0;
        st->substH[j]=0;
        st->substL[j]=0;
    }
    for (j=0;j<NUM_CNTXT;j++) {
        st->cntxt[j] = (long *)smalloc(st->nbins * sizeof(long));
        for (i=0;i<st->nbins;i++)
            st->cntxt[j][i] = 0;
        st->cntxtH[j]=0;
	st->cntxtL[j]=0;
        st->homop[j] = (long *)smalloc(st->nbins * sizeof(long));
        for (i=0;i<st->nbins;i++)
            st->homop[j][i] = 0;
        st->homopH[j]=0;
	st->homopL[j]=0;
    }
#endif

    st->total_bases = 0;
    st->total_errors = 0;

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
                    free(st->predictor);
                    free(st->num_bases);
                    free(st->num_errors);
#ifdef ST_STRUCTURE
                    int j;
                    for (j=0;j<NUM_SUBST;j++)
                        free(st->subst[j]);
                    for (j=0;j<NUM_CNTXT;j++)
                        free(st->cntxt[j]);
                    for (j=0;j<NUM_CNTXT;j++)
                        free(st->homop[j]);
#endif                    
                    st->nbins = 0;
                }
            }
        }

#ifdef ST_STRUCTURE
    int j;
    for(read=0;read<N_READS;read++)
    {
        SurvTable *st = sts[(N_TILES+1)*N_READS+read];
        if( NULL == st) continue;
        if (st->nbins) {
            free(st->predictor);
            free(st->num_bases);
            free(st->num_errors);
            for (j=0;j<NUM_SUBST;j++)
                free(st->subst[j]);
            for (j=0;j<NUM_CNTXT;j++)
                free(st->cntxt[j]);
            for (j=0;j<NUM_CNTXT;j++)
                free(st->homop[j]);
            st->nbins = 0;
        }
    }
#endif                    

}

static int optimalPredictorBin(SurvTable *st)
{
    int ipopt = 0;
    float max_diff = 0.0;
    long cum_bases = st->total_bases;
    long cum_errors = st->total_errors;
    int i;

    /* find optimal predictor which gives the best separation between bases and errors */
    for(i=0;i<st->nbins;i++)
    {
        float frac_bases, frac_errors, diff_frac;
        frac_bases = (float)cum_bases/(float)st->total_bases;
        frac_errors = (float)cum_errors/(float)st->total_errors;
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
            int ipopt;

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

            ipopt = optimalPredictorBin(st);
            st->optimal_predictor = st->predictor[ipopt];

#ifdef ST_STRUCTURE
            int j;
            for(i=0;i<st->nbins;i++)
            {
                if (st->predictor[i] > st->optimal_predictor)
                {
                    for(j=0;j<NUM_SUBST;j++)
                        st->substH[j] += st->subst[j][i];
                    for(j=0;j<NUM_CNTXT;j++)
                        st->cntxtH[j] += st->cntxt[j][i];
                    for(j=0;j<NUM_CNTXT;j++)
                        st->homopH[j] += st->homop[j][i];
                }
                else
                {
                    for(j=0;j<NUM_SUBST;j++)
                        st->substL[j] += st->subst[j][i];
                    for(j=0;j<NUM_CNTXT;j++)
                        st->cntxtL[j] += st->cntxt[j][i];
                    for(j=0;j<NUM_CNTXT;j++)
                        st->homopL[j] += st->homop[j][i];
                }
            }
#endif                
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
#ifdef ST_STRUCTURE
                    int j;
                    for(j=0;j<NUM_SUBST;j++)
                        for(i=0;i<st->nbins;i++)
                            st->subst[j][i] += tile_st->subst[j][i];
                    for(j=0;j<NUM_CNTXT;j++)
                        for(i=0;i<st->nbins;i++)
                            st->cntxt[j][i] += tile_st->cntxt[j][i];
                    for(j=0;j<NUM_CNTXT;j++)
                        for(i=0;i<st->nbins;i++)
                            st->homop[j][i] += tile_st->homop[j][i];
#endif
                    tile_st->total_bases = 0;
                    tile_st->total_errors = 0;
                }
            }
        }
    }

    completeSurvTable(s, &sts[ntiles*N_READS], 0);

#ifdef ST_STRUCTURE
    int j;
    for(read=0;read<N_READS;read++)
    {
        read_length = s->read_length[read];
        if (0 == read_length) continue;
        sts[(N_TILES+1)*N_READS+read] = smalloc(sizeof(SurvTable));

        SurvTable *st = sts[(N_TILES+1)*N_READS+read];

        initialiseSurvTable(s, st, -1, read, -1);

        for(cycle=0;cycle<read_length;cycle++)
        {
            SurvTable *cycle_st = sts[ntiles*N_READS+read] + cycle;
            for(i=0;i<st->nbins;i++)
            {
                st->num_bases[i] += cycle_st->num_bases[i];
                st->num_errors[i] += cycle_st->num_errors[i];
            }
            for(j=0;j<NUM_SUBST;j++)
                for(i=0;i<st->nbins;i++)
                    st->subst[j][i] += cycle_st->subst[j][i];
            for(j=0;j<NUM_CNTXT;j++)
                for(i=0;i<st->nbins;i++)
                    st->cntxt[j][i] += cycle_st->cntxt[j][i];
            for(j=0;j<NUM_CNTXT;j++)
                for(i=0;i<st->nbins;i++)
                    st->homop[j][i] += cycle_st->homop[j][i];
        }
    }

    completeSurvTable(s, &sts[(N_TILES+1)*N_READS], 1);
#endif
    
    return;
}
    
#ifdef ST_STRUCTURE
static void outputErrorTable(Settings *s, int ntiles, SurvTable **sts)
{
    FILE *fp = NULL;
    int read, cycle, i, j;
    char p;

    // open error file
    for(read=0;read<N_READS;read++)
    {
        if( NULL == sts[ntiles*N_READS+read]) continue;
        for(cycle=0;cycle<s->read_length[read];cycle++)
        {
            SurvTable *st = sts[ntiles*N_READS+read] + cycle;

            // skip st with no data
            if( 0 == st->total_bases ) continue;

            if( NULL == fp )
            {
                int filename_sz;
                char *filename;
                filename_sz = (NULL == s->prefix ? 0 : strlen(s->prefix)) + 100;
                filename = smalloc(filename_sz);
                sprintf(filename, "%s_%s_error.txt", s->prefix, (st->mode == ST_MODE_PURITY ? "purity" : "quality"));
                fp = fopen(filename, "w");
                if( NULL == fp )
                {
                    fprintf(stderr, "ERROR: can't open Error file %s: %s\n",
                            filename, strerror(errno));
                    exit(EXIT_FAILURE);
                }
                free(filename);
            }
        }
    }
    if( NULL == fp )
    {
        fprintf(stderr, "ERROR: no error data\n");
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "# Substitution error table. Use `grep ^SET | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per predictor, columns read, predictor followed by substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        SurvTable *st = sts[(N_TILES+1)*N_READS+read];
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
    
    fprintf(fp, "# Mismatch substitutions high predictor. Use `grep ^MSH | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and cycle, columns read, cycle then substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        if( NULL == sts[ntiles*N_READS+read]) continue;
        for(cycle=0;cycle<s->read_length[read];cycle++)
        {
            SurvTable *st = sts[ntiles*N_READS+read] + cycle;
            if( 0 == st->total_bases ) continue;
            fprintf(fp, "MSH\t%d\t%d", read, cycle);
            for (j=0;j<NUM_SUBST;j++)
            {
                char *subst;
                subst=word2str(j,LEN_SUBST);
                if (subst[0]==subst[1]) continue;
                fprintf(fp, "\t%s\t%ld", subst, st->substH[j]);
            }
            fprintf(fp, "\n");
        }
    }

    fprintf(fp, "# Mismatch substitutions low predictor. Use `grep ^MSL | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and cycle, columns read, cycle then substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        if( NULL == sts[ntiles*N_READS+read]) continue;
        for(cycle=0;cycle<s->read_length[read];cycle++)
        {
            SurvTable *st = sts[ntiles*N_READS+read] + cycle;
            if( 0 == st->total_bases ) continue;
            fprintf(fp, "MSL\t%d\t%d", read, cycle);
            for (j=0;j<NUM_SUBST;j++)
            {
                char *subst;
                subst=word2str(j,LEN_SUBST);
                if (subst[0]==subst[1]) continue;
                fprintf(fp, "\t%s\t%ld", subst, st->substL[j]);
            }
            fprintf(fp, "\n");
        }
    }

    fprintf(fp, "# Substitution profile high predictor. Use `grep ^PRH | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read, columns read then substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        SurvTable *st = sts[(N_TILES+1)*N_READS+read];
        if( NULL == st) continue;

        fprintf(fp, "PRH\t%d", read);
        for (j=0;j<NUM_SUBST;j++)
        {
            char *subst;
            subst=word2str(j,LEN_SUBST);
            if (subst[0]==subst[1]) continue;
            fprintf(fp, "\t%s\t%ld", subst, st->substH[j]);
        }
        fprintf(fp, "\n");
    }
    
    fprintf(fp, "# Substitution profile low predictor. Use `grep ^PRL | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read, columns read then substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        SurvTable *st = sts[(N_TILES+1)*N_READS+read];
        if( NULL == st) continue;

        fprintf(fp, "PRL\t%d", read);
        for (j=0;j<NUM_SUBST;j++)
        {
            char *subst;
            subst=word2str(j,LEN_SUBST);
            if (subst[0]==subst[1]) continue;
            fprintf(fp, "\t%s\t%ld", subst, st->substL[j]);
        }
        fprintf(fp, "\n");
    }
    
    fprintf(fp, "# Effect of previous base high predictor. Use `grep ^P1H | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and previous base, columns read then previous base+substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        SurvTable *st = sts[(N_TILES+1)*N_READS+read];
        if( NULL == st) continue;

        for (j=0,p='\0';j<NUM_CNTXT;j++)
        {
            char *cntxt;
            cntxt=word2str(j,LEN_CNTXT);
            if (p != cntxt[0]) {
                p = cntxt[0];
                if (j) fprintf(fp, "\n");
                fprintf(fp, "P1H\t%d", read);
            }
            if (cntxt[1]==cntxt[2]) continue;
            fprintf(fp, "\t%s\t%ld", cntxt, st->cntxtH[j]);
        }
        fprintf(fp, "\n");
    }
    
    fprintf(fp, "# Effect of previous base low predictor. Use `grep ^P1L | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and previous base, columns read then previous base+substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        SurvTable *st = sts[(N_TILES+1)*N_READS+read];
        if( NULL == st) continue;

        for (j=0, p='\0';j<NUM_CNTXT;j++)
        {
            char *cntxt;
            cntxt=word2str(j,LEN_CNTXT);
            if (p != cntxt[0]) {
                p = cntxt[0];
                if (j) fprintf(fp, "\n");
                fprintf(fp, "P1L\t%d", read);
            }
            if (cntxt[1]==cntxt[2]) continue;
            fprintf(fp, "\t%s\t%ld", cntxt, st->cntxtL[j]);
        }
        fprintf(fp, "\n");
    }
    
    fprintf(fp, "# Homopolymer effect high predictor. Use `grep ^HPH | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and homopolymer, columns read then homopolymer+substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        SurvTable *st = sts[(N_TILES+1)*N_READS+read];
        if( NULL == st) continue;

        for (j=0,p='\0';j<NUM_CNTXT;j++)
        {
            int k;
            char homop[LEN_HOMOP];
            char *cntxt;
            cntxt=word2str(j,LEN_CNTXT);
            if (p != cntxt[0]) {
                p = cntxt[0];
                if (j) fprintf(fp, "\n");
                fprintf(fp, "HPH\t%d", read);
            }
            if (cntxt[1]==cntxt[2]) continue;
            for (k=0;k<(LEN_HOMOP-1);k++)
                homop[k]=cntxt[0];
            homop[k]='\0';
            fprintf(fp, "\t%s%s\t%ld", homop, cntxt, st->homopH[j]);
        }
        fprintf(fp, "\n");
    }
    
    fprintf(fp, "# Homopolymer effect low predictor. Use `grep ^HPL | cut -f 2-` to extract this part\n");
    fprintf(fp, "# One row per read and homopolymer, columns read then homopolymer+substitution and count for 12 substitutions\n");
    for(read=0;read<N_READS;read++)
    {
        SurvTable *st = sts[(N_TILES+1)*N_READS+read];
        if( NULL == st) continue;

        for (j=0,p='\0';j<NUM_CNTXT;j++)
        {
            int k;
            char homop[LEN_HOMOP];
            char *cntxt;
            cntxt=word2str(j,LEN_CNTXT);
            if (p != cntxt[0]) {
                p = cntxt[0];
                if (j) fprintf(fp, "\n");
                fprintf(fp, "HPL\t%d", read);
            }
            if (cntxt[1]==cntxt[2]) continue;
            for (k=0;k<(LEN_HOMOP-1);k++)
                homop[k]=cntxt[0];
            homop[k]='\0';
            fprintf(fp, "\t%s%s\t%ld", homop, cntxt, st->homopL[j]);
        }
        fprintf(fp, "\n");
    }
    
    if( NULL != fp) fclose(fp);
}
#endif

static void outputSurvTable(Settings *s, SurvTable **sts)
{
    FILE *fp = NULL;
    int itile, read, cycle, i;

    for(itile=0;itile<=N_TILES;itile++)
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

    ct->predictor = smalloc(ct->nbins * sizeof(float));
    for(i=0;i<ct->nbins;i++)
        ct->predictor[i] = st->predictor[i];

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
                    free(ct->predictor);
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

static void optimisePredictorBins(Settings *s, SurvTable *st, CalTable *ct)
{
    float predictor[st->nbins];
    int nbins = 0;

    int i, j;
    int ipopt = 0, ipmax = 0;
    int ipbin;
    float pinc;

    ipopt = optimalPredictorBin(st);
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
        /* add n_bins_left bins between first bin and popt */
        pinc = (st->predictor[ipopt] - st->predictor[0]) / (s->n_bins_left + 1);
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
    }

    if (ipmax > ipopt)
    {
        /* add n_bins_right bins between popt and pmax */
        pinc = (st->predictor[ipmax] - st->predictor[ipopt]) / (s->n_bins_right + 1);
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

                optimisePredictorBins(s, st, ct);

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
    FILE *fp = NULL;
    int itile, read, cycle, i;

    for(itile=0;itile<=N_TILES;itile++)
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
    int iregion = -1;
    int c, b;

    if (NULL != cif_data) assert((cstart + read_length) <= cif_data->ncycles);

    if (s->spatial_filter) iregion = xy2region(x, y, s->region_size, s->nregions_x, s->nregions_y);

#ifdef CALDATA
    int dir = BAM_FREVERSE & bam->core.flag ? 1 : 0;
    char *chrom = fp->header->target_name[bam->core.tid];
    int pos = bam->core.pos;

    fprintf(fp_caldata, "%d\t%d\t%d\t%s\t%d\t", x, y, dir, chrom, pos);
#endif

    /* update survival table */
    for (c = cstart, b = 0; b < read_length; c++, b++) {
        SurvTable *st;
        float predictor = -1.0;
        int ibin;

        /* set cycle st */
        if (s->spatial_filter) {
            int state = (getFilterData(tile, read, b, iregion) & REGION_STATE_MISMATCH) ? 1 : 0;
            st = sts[state*N_READS+read] + b;
        }else{
            st = sts[read] + b;
        }

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
#ifdef ST_STRUCTURE
            if( read_mismatch[b] & BASE_MISMATCH ){
                char subst[LEN_SUBST+1];
                char cntxt[LEN_CNTXT+1];
                int word;

                subst[0]=read_ref[b];
                subst[1]=read_seq[b];
                word=str2word(subst,LEN_SUBST);
                if( word >= 0 )
                    st->subst[word][ibin]++;

                cntxt[0]=(b > 0 ? read_ref[b-1] : 'N');
                cntxt[1]=read_ref[b];
                cntxt[2]=read_seq[b];
                word=str2word(cntxt,LEN_CNTXT);
                if( word >= 0 ) 
                    st->cntxt[word][ibin]++;

                if( word >= 0 && b >= LEN_HOMOP ){
                    int l;
                    for(l=1; l<=LEN_HOMOP; l++)
                        if( read_ref[b-l] != cntxt[0] )
                            break;
                    if( l > LEN_HOMOP)
                        st->homop[word][ibin]++;
                }
            }
#endif
        }

#ifdef CALDATA
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
        if( st->mode == ST_MODE_PURITY ){
            fprintf(fp_caldata, "%d %d %d %d %d %d %d %d\t",
                    ref_base, seq_base, read_qual[b], bin[0], bin[1], bin[2], bin[3], read_mismatch[b]);
        } else {
            fprintf(fp_caldata, "%d %d %d %d\t",
                    ref_base, seq_base, read_qual[b], biread_mismatch[b]);
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
 * Returns: +'ve integer for success
 *          0 for failure
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

    int itile, read;

    bam1_t *bam = bam_init1();

    if (NULL != s->intensity_dir) checked_chdir(s->intensity_dir);

    for(itile=0;itile<=N_TILES;itile++)
        for(read=0;read<N_READS;read++)
            sts[itile*N_READS+read] = NULL;
#ifdef ST_STRUCTURE
    for(read=0;read<N_READS;read++)
        sts[(N_TILES+1)*N_READS+read] = NULL;
#endif
    
    itile = -1;

    /* loop over reads in the bam file */
    while (1){
        int bam_lane = -1, bam_tile = -1, bam_x = -1, bam_y = -1, bam_read = -1, read_length;
        size_t bam_offset = 0;
        int cycle;

        if (parse_bam_readinfo(fp_bam, bam, &bam_lane, &bam_tile, &bam_x, &bam_y, &bam_read, (NULL == s->intensity_dir ? NULL : &bam_offset))) {
            break;	/* break on end of BAM file */
		}

        if (BAM_FUNMAP & bam->core.flag) continue;
        if (BAM_FQCFAIL & bam->core.flag) continue;
        if (BAM_FPAIRED & bam->core.flag) {
            if (0 == (BAM_FPROPER_PAIR & bam->core.flag)) {
                continue;
            }
        }

        parse_bam_alignments(fp_bam, bam, bam_read_seq, bam_read_qual, bam_read_ref,
                                     bam_read_mismatch, bam_read_buff_size,s->snp_hash);

        read_length = strlen(bam_read_seq);
        if (0 == s->read_length[bam_read]) {
            s->read_length[bam_read] = read_length;

            if (s->spatial_filter) {
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
   	    size_t nelem = ntiles_bam;
	    void *pitile;

            tile = bam_tile;
            
	    pitile = lfind(&tile, tiles, &nelem, sizeof(int), &int_cmp);
	    if (NULL == pitile) {
                itile = ntiles_bam;
                if(++ntiles_bam > N_TILES){
                    fprintf(stderr,"ERROR: too many tiles %d > %d.\n", ntiles_bam, N_TILES);
                    exit(EXIT_FAILURE);
                }
                tiles[itile] = tile;
                if (!s->quiet)fprintf(stderr, "Processing tile %i (%d)\n", tile, nreads_bam);
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

        if (0 == s->spatial_filter) {
            /* if we don't have a filter file we allocate memory for this tile now */
            if (NULL == sts[itile*N_READS+bam_read]) {
                sts[itile*N_READS+bam_read] = smalloc(read_length * sizeof(SurvTable));
                nst += read_length;
                for(cycle=0;cycle<read_length;cycle++)
                    initialiseSurvTable(s, sts[itile*N_READS+bam_read]+cycle, tile, bam_read, cycle);
            }
        }
        
        if (0 != updateSurvTable(s, (s->spatial_filter ? sts : &sts[itile*N_READS]), cif_data,
                                 bam_offset, bam_tile, bam_x, bam_y, bam_read, bam_read_mismatch,
                                 bam_read_seq, bam_read_qual, bam_read_ref, fp_bam, bam, fp_caldata)) {
            fprintf(stderr,"ERROR: updating quality values for tile %i.\n", tile);
            exit(EXIT_FAILURE);
        }
        
        nreads_bam++;
    }
    
    if (s->spatial_filter) {
        int state;
        for(state=0;state<N_STATES;state++)
            completeSurvTable(s, &sts[state*N_READS], 0);
    }else{
        for(itile=0;itile<ntiles_bam;itile++)
            completeSurvTable(s, &sts[itile*N_READS], 0);
    }
    

    if (NULL != cif_data) free_cif_data(cif_data);

    if (NULL != fp_caldata) fclose(fp_caldata);

    bam_destroy1(bam);
    
    *ntiles = ntiles_bam;
    *nreads = nreads_bam;

    return nst;
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
    int nreads = 0;
    int nst = 0;
#ifdef ST_STRUCTURE
    SurvTable *sts[(N_TILES+2)*N_READS];
#else
    SurvTable *sts[(N_TILES+1)*N_READS];
#endif
    int nct = 0;
    CalTable *cts[(N_TILES+1)*N_READS];

    settings.prefix = NULL;
    settings.quiet = 0;
    settings.filter_bad_tiles = 0;
    settings.spatial_filter = 0;
    settings.region_size = 0;
    settings.nregions_x = 0;
    settings.nregions_y = 0;
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
    nst = makeSurvTable(&settings, fp_bam, sts, &ntiles, &nreads);
    if (0 == nst) {
        fprintf(stderr,"WARNING: No data in BAM file\n");
        exit(EXIT_SUCCESS);
    }

    if (!settings.quiet) {
        fprintf(stderr, "Processed %8d traces\n", nreads);
        if (NULL != settings.snp_hash) {
            HashIter *tileIter = HashTableIterCreate();
            HashItem *hashItem;
            size_t nsnps = 0;
            while ((hashItem = HashTableIterNext(settings.snp_hash, tileIter)))
                nsnps += hashItem->data.i;
            fprintf(stderr, "Ignored %lu snps\n", nsnps);
        }
    }

    /* back to where we belong */
    checked_chdir(settings.working_dir);

    if (!settings.spatial_filter) makeGlobalSurvTable(&settings, ntiles, sts);

    outputSurvTable(&settings, sts);

#ifdef ST_STRUCTURE
    outputErrorTable(&settings, ntiles, sts);
#endif    

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

    if (NULL != settings.working_dir) free(settings.working_dir);

    return EXIT_SUCCESS;

}
