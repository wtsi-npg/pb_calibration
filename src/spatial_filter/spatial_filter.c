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
 *   produce ascii plots for read TILEVIZ like those produced by tileviz
 */

//#define QC_FAIL
//#define PROPERLY_PAIRED
//#define TILEVIZ 1

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

#define PBP_VERSION PACKAGE_VERSION

const char * pbs_c_rev = "$Revision$";

#define N_LANES 8
#define N_TILES 120
#define N_CYCLES 250

#define BASE_ALIGN      (1<<0)
#define BASE_MISMATCH   (1<<1)
#define BASE_INSERTION  (1<<2)
#define BASE_DELETION   (1<<3)
#define BASE_SOFT_CLIP  (1<<4)
#define BASE_KNOWN_SNP  (1<<5)

#define PHRED_QUAL_OFFSET 33  // phred quality values offset

#define COORD_SHIFT       1000.0
#define COORD_FACTOR      10.0
#define REDUCTION_FACTOR  200.0

#define REGION_MISMATCH_THRESHOLD   50.0
#define REGION_INSERTION_THRESHOLD  20.0
#define REGION_DELETION_THRESHOLD   20.0
#define REGION_SOFT_CLIP_THRESHOLD  50.0

#define FILTER_PF         (1<<0)
#define FILTER_COVERAGE   (1<<1)
#define FILTER_MISMATCH   (1<<2)
#define FILTER_INSERTION  (1<<3)
#define FILTER_DELETION   (1<<4)
#define FILTER_SOFT_CLIP  (1<<5)

#define FILTER_MASK  (FILTER_INSERTION | FILTER_DELETION)  // treat reads with ANY indels as a QCFAIL

typedef struct {
    int         lane;
    int         tile;
    int         read;
    int         cycle;
    int         nregions;
    long        *align;
    long        *mismatch;
    long        *insertion;
    long        *deletion;
    long        *soft_clip;
    long        *known_snp;
    long        *no_call;
    float       quality;
} RegionTable;

typedef struct {
    char *cmdline;
    char *prefix;
    char *intensity_dir;
    char *snp_file;
    HashTable *snp_hash;
    char *working_dir;
    char *output;
    int read_length[3];
    int dump;
    int calculate;
    int apply;
    int qcfail;
    int quiet;
    int nregions_x;
    int nregions_y;
    int nregions;
} Settings;

static
void
checked_chdir(const char* dir){
    if(chdir(dir)){
        fprintf(stderr,"ERROR: failed to change directory to: %s\n",dir);
        exit(EXIT_FAILURE);
    }
}

static char *complement_table = NULL;

static void setRegions(Settings *s)
{
    char *file = NULL;
    size_t file_sz = strlen(s->intensity_dir) + 30;
    int fd = -1;
    ssize_t res;
    uint32_t width, height;
    
    file = smalloc(file_sz);
    snprintf(file, file_sz, "%s/../ImageSize.dat", s->intensity_dir);
    fd = open(file, O_RDONLY);
    if (fd < 0) {
        die("Couldn't open %s : %s\n", file, strerror(errno));
    }

    res = read(fd, &width, 4);
    if (4 != res) {
       die("Error reading %s\n", file);
    }
    res = read(fd, &height, 4);
    if (4 != res) {
       die("Error reading %s\n", file);
    }

    close(fd);

    s->nregions_x = 1 + (int)(width / REDUCTION_FACTOR);
    s->nregions_y = 1 + (int)(height / REDUCTION_FACTOR);
    s->nregions = s->nregions_x * s->nregions_y;
    
    fprintf(stderr, "Read tile width=%u height=%u from file %s\n", width, height, file);
    fprintf(stderr, "nregions_x=%d nregions_y=%d nregions=%d\n", s->nregions_x, s->nregions_y, s->nregions);

    free(file);

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

static void initialiseRegionTable(Settings *s, RegionTable *rt, int lane, int tile, int read, int cycle)
{
    int i;

    rt->lane  = lane;
    rt->tile  = tile;
    rt->read  = read;
    rt->cycle = cycle;

    rt->nregions = s->nregions;

    rt->align = smalloc(rt->nregions * sizeof(long));
    rt->mismatch = smalloc(rt->nregions * sizeof(long));
    rt->insertion = smalloc(rt->nregions * sizeof(long));
    rt->deletion = smalloc(rt->nregions * sizeof(long));
    rt->soft_clip = smalloc(rt->nregions * sizeof(long));
    rt->known_snp = smalloc(rt->nregions * sizeof(long));
    for (i=0;i<rt->nregions;i++) {
        rt->align[i] = 0;
        rt->mismatch[i] = 0;
        rt->insertion[i] = 0;
        rt->deletion[i] = 0;
        rt->soft_clip[i] = 0;
        rt->known_snp[i] = 0;
    }
}

static void freeRegionTable(RegionTable *rt)
{
    if (rt->nregions) {
        free(rt->align);
        free(rt->mismatch);
        free(rt->insertion);
        free(rt->deletion);
        free(rt->soft_clip);
        free(rt->known_snp);

        rt->nregions = 0;
    }
}

static uint8_t * load_filter(Settings *s, char *file, size_t *nclusters) {
    int fd = -1;
    ssize_t res;
    uint8_t *filter = NULL;
    
    uint32_t check, version, nfilter;
    
    fd = open(file, O_RDONLY);
    if (fd < 0) {
        die("Couldn't open %s : %s\n", file, strerror(errno));
    }

    res = read(fd, &check, 4);
    if (4 != res) {
       die("Error reading %s\n", file);
    }
    if (0 != check) {
        die("Invalid filter file (%s) check=%ld expected 0\n", file, check);
    }
    res = read(fd, &version, 4);
    if (4 != res) {
       die("Error reading %s\n", file);
    }
    if (3 != version) {
        die("Invalid filter file (%s) version=%ld expected 0\n", file, version);
    }
    res = read(fd, &nfilter, 4);
    if (4 != res) {
       die("Error reading %s\n", file);
    }

    filter = smalloc(nfilter * sizeof(uint8_t));

    res = read(fd, filter, nfilter);
    if (res < 0) {
        die("Error reading %s: %s\n", file, strerror(errno));
    }
    if (res != nfilter) {
        die("Error: Did not get the expected amount of data from %s\n", file);
    }

    *nclusters = nfilter;

    return filter;
}

static void output_filter(Settings *s, int lane, int tile, size_t nclusters, uint8_t *filter) {
    char *file = NULL;
    size_t file_sz = strlen(s->working_dir) + (NULL == s->prefix ? 0 : strlen(s->prefix)) + 100;
    uint32_t check = 0;
    uint32_t version = 3;
    int fd = -1;
    ssize_t res;
    
    file = smalloc(file_sz);
    snprintf(file, file_sz, "%s/%s_%04d.filter", s->working_dir, s->prefix, tile);
    fd = open(file, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) {
        fprintf(stderr,"ERROR: cannot open filter file: %s\n",file);
        exit(EXIT_FAILURE);
    }

    res = write(fd, &check, 4);
    if (4 != res) {
        die("Error writing %s\n", file);
    }
    res = write(fd, &version, 4);
    if (4 != res) {
        die("Error writing %s\n", file);
    }
    res = write(fd, &nclusters, 4);
    if (4 != res) {
        die("Error writing %s\n", file);
    }
    res = write(fd, filter, nclusters);
    if (res != nclusters) {
        die("Error writing %s\n", file);
    }

    return;
}

/*
 * Fast lookup table to convert integer fraction value to a floating point
 * value. This allows floats to be parsed to 2 decimal points with no
 * logs or divisions anywhere. Over-optimising? Maybe, but I have the bit
 * between my teeth and am running with it :-)
 */
static float f_lookup[200] = {
    /* Filler */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90,
    /* Filler */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* Filler */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* Filler */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* Filler */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* Filler */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* Filler */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* Filler */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* Filler */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
    0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19,
    0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,
    0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,
    0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49,
    0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59,
    0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69,
    0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79,
    0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89,
    0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
};

/*
 * simplification of parse_4_float for single float parse
 */
static
const char *
parse_next_float(const char *str, float *val, int *bin)
{
    int minus = 0;
    char c;
    double fval = 0;
    int ival1 = 0, ival2 = 1;
    int is_use_fval = 0;

    *val = 0.;

    c = *str++;
    while (isspace(c))
        c = *str++;
    if (!c)
        return NULL;

    if (c == '-' || c == '+') {
        minus = (c == '-');
        c = *str++;
    }

    while (c >= '0' && c <= '9') {
        ival1 = ival1 * 10 + c-'0';
        c = *str++;
    }

    if (c == '.') {
        c = *str++;
        while (c >= '0' && c <= '9') {
            ival2 = ival2 * 10 + c-'0';
            c = *str++;
        }

        while (ival2 >= 1000)
            ival2 /= 10;

        assert((ival2 >=  10 && ival2 <  20) ||
               (ival2 >= 100 && ival2 < 200));

        fval = ival1 + f_lookup[ival2];

        if (minus)
            fval = -fval;

        is_use_fval = 1;
    }

    switch(c) {
    case '\n': case '\r': case '\t': case ' ': case '\0':
        if (minus)
            ival1 = -ival1;

        if (bin) {
            if (ival1 > 65535)
                bin[65535]++;
            else if (ival1 < -65535)
                bin[-65535]++;
            else
                bin[ival1]++;
        }

        if( is_use_fval ) { *val = fval; }
        else              { *val = ival1; }
        return str;
    }

    /*
     * If we get here we've run out of input before processing float value
     * (the while command terminated) or we've found malformed data and
     * jumped here directly.
     */
    fprintf(stderr, "Error: unexpected character '%c' during parsing "
	    "of string \"%s\"\n", c, str);
    exit(EXIT_FAILURE);
}

/* get next float x,y values from pos file */
static
void
get_next_pos(FILE * fp_pos,
             float* pos_x_ptr,
             float* pos_y_ptr){
    static const int line_size = 1024;
    char line[line_size];
    const char* cp = line;


    if(NULL == fgets(line,line_size,fp_pos)){
        fprintf(stderr,"ERROR: can't parse line from pos.txt file\n");
        exit(EXIT_FAILURE);
    }
    cp = parse_next_float(cp, pos_x_ptr, 0);
    parse_next_float(cp, pos_y_ptr, 0);
}

static int* load_positions(Settings *s, int lane, int tile, size_t nclusters) {
    char *file = NULL;
    size_t file_sz = strlen(s->intensity_dir) + 30;
    FILE *fp = NULL;
    size_t npositions = 0;
    int *iregions = NULL;
    
    file = smalloc(file_sz);
    snprintf(file, file_sz, "%s/s_%d_%d_pos.txt", s->intensity_dir, lane, tile);
        
    if(NULL == (fp = fopen(file, "r"))){
        fprintf(stderr,"ERROR: cannot open pos file: %s\n",file);
        exit(EXIT_FAILURE);
    }

    iregions = smalloc(nclusters * sizeof(int));

    while(npositions < nclusters) {
        float x, y;
        int iregion;
        get_next_pos(fp, &x, &y);
        iregion = (int)(x / REDUCTION_FACTOR) * s->nregions_y + (int)(y / REDUCTION_FACTOR);
        iregions[npositions++] = iregion;
    }

    fclose(fp);
    
    free(file);

    return iregions;
}

static void completeRegionTables(Settings *s, RegionTable *rts)
{
    float ssc = 1.0;
    int read, read_length, irt, cycle, i;

    if( NULL == rts )
        return;

    for(read=0;read<3;read++)
    {
        read_length = s->read_length[read];
        irt = (read > 1 ? 1 : 0) * N_CYCLES;
        for(cycle = 0; cycle < read_length; cycle++, irt++)
        {
            RegionTable *rt = rts + irt;
            long quality_bases = 0;
            long quality_errors = 0;

            for(i=0;i<rt->nregions;i++)
            {
                quality_bases  += rt->align[i];
                quality_errors += rt->mismatch[i];
            }

            rt->quality = -10.0 * log10((quality_errors + ssc)/(quality_bases + ssc));
        }
    }
}

/*
 * identify bad tiles quality < mean quality - filter * stdev quality
 * filter is typically 2.
*/

static void findBadRegions(Settings *s, int ntiles, RegionTable **rts)
{
    char *file = NULL;
    size_t file_sz = strlen(s->intensity_dir) + 100;
    int itile;
    
    if (0 >= ntiles || NULL == rts)
        return;
    
    file = smalloc(file_sz);

    for(itile=0;itile<ntiles;itile++)
    {
        int lane = rts[itile][0].lane;
        int tile = rts[itile][0].tile;
        RegionTable *rt = rts[itile] + 2 * N_CYCLES;
        size_t nclusters;
        uint8_t *filter = NULL;
        int *iregions = NULL;
        int read, read_length, irt, cycle, iregion;
        size_t icluster;

#ifdef TILEVIZ
        // generate tileviz like summary plots for read TILEVIZ

        initialiseRegionTable(s, rt, lane, tile, -1, -1);

        read = TILEVIZ;
        read_length = s->read_length[read];
        irt = (read > 1 ? 1 : 0) * N_CYCLES;
        for(cycle = 0; cycle < read_length; cycle++, irt++)
        {
            RegionTable *cycle_rt = rts[itile] + irt;
            for (iregion=0; iregion<rt->nregions; iregion++)
            {
                int n = cycle_rt->align[iregion]
                      + cycle_rt->insertion[iregion]
                      + cycle_rt->deletion[iregion]
                      + cycle_rt->soft_clip[iregion]
                      + cycle_rt->known_snp[iregion];


                // coverage should be the same for all cycles
                rt->align[iregion] = n;

                // for all other values take the maximum over all cycles
                if( cycle_rt->mismatch[iregion] > rt->mismatch[iregion] )
                    rt->mismatch[iregion] = cycle_rt->mismatch[iregion];
                if( cycle_rt->insertion[iregion] > rt->insertion[iregion] )
                    rt->insertion[iregion] = cycle_rt->insertion[iregion];
                if( cycle_rt->deletion[iregion] > rt->deletion[iregion] )
                    rt->deletion[iregion] = cycle_rt->deletion[iregion];
                if( cycle_rt->soft_clip[iregion] > rt->soft_clip[iregion] )
                    rt->soft_clip[iregion] = cycle_rt->soft_clip[iregion];
                if( cycle_rt->known_snp[iregion] > rt->known_snp[iregion] )
                    rt->known_snp[iregion] = cycle_rt->known_snp[iregion];
            }
        }

        int i, j;

        fprintf(stderr, "tile=%-4d align          mismatch       insertion       deletion      soft_clip\n", tile);
        for (i=0; i<s->nregions_y; i++)
        {
            fprintf(stderr, "          ");
            // truncate coverage to 15 so it displays as a single hex value
            for (j=0; j<s->nregions_x; j++) {
                int k = j * s->nregions_y + i;
                int v = rt->align[k] > 15 ? 15 : rt->align[k];
                fprintf(stderr, "%1x",v);
            }
            fprintf(stderr, "    ");
            // convert all other values to a percentage and bin 0(0%), 1(10%), .. 10(100%)
            // correct for sparse bins by assuming ALL bins have atleast 10 clusters
            for (j=0; j<s->nregions_x; j++) {
                int k = j * s->nregions_y + i;
                int n = rt->align[k] < 10 ? 10 : rt->align[k];
                int v = (int)(10.0 * rt->mismatch[k] / n);
                fprintf(stderr, "%1x",v);
            }
            fprintf(stderr, "    ");
            for (j=0; j<s->nregions_x; j++) {
                int k = j * s->nregions_y + i;
                int n = rt->align[k] < 10 ? 10 : rt->align[k];
                int v = (int)(10.0 * rt->insertion[k] / n);
                fprintf(stderr, "%1x",v);
            }
            fprintf(stderr, "    ");
            for (j=0; j<s->nregions_x; j++) {
                int k = j * s->nregions_y + i;
                int n = rt->align[k] < 10 ? 10 : rt->align[k];
                int v = (int)(10.0 * rt->deletion[k] / n);
                fprintf(stderr, "%1x",v);
            }
            fprintf(stderr, "    ");
            for (j=0; j<s->nregions_x; j++) {
                int k = j * s->nregions_y + i;
                int n = rt->align[k] < 10 ? 10 : rt->align[k];
                int v = (int)(10.0 * rt->soft_clip[k] / n);
                fprintf(stderr, "%1x",v);
            }
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "\n");

        freeRegionTable(rt);
#endif

        initialiseRegionTable(s, rt, lane, tile, -1, -1);
        
        for(read=0;read<3;read++)
        {
            read_length = s->read_length[read];
            irt = (read > 1 ? 1 : 0) * N_CYCLES;
            for(cycle = 0; cycle < read_length; cycle++, irt++)
            {
                RegionTable *cycle_rt = rts[itile] + irt;
#if 1
                fprintf(stderr, "tile=%4d read=%1d cycle=%-3d quality=%.2f\n", tile, cycle_rt->read, cycle_rt->cycle, cycle_rt->quality);
#endif
                for (iregion=0; iregion<rt->nregions; iregion++)
                {
                    int n = cycle_rt->align[iregion]
                          + cycle_rt->insertion[iregion]
                          + cycle_rt->deletion[iregion]
                          + cycle_rt->soft_clip[iregion]
                          + cycle_rt->known_snp[iregion];


                    // coverage should be the same for all cycles
                    rt->align[iregion] = n;

                    // for all other values take the maximum over all cycles
                    if( cycle_rt->mismatch[iregion] > rt->mismatch[iregion] )
                        rt->mismatch[iregion] = cycle_rt->mismatch[iregion];
                    if( cycle_rt->insertion[iregion] > rt->insertion[iregion] )
                        rt->insertion[iregion] = cycle_rt->insertion[iregion];
                    if( cycle_rt->deletion[iregion] > rt->deletion[iregion] )
                        rt->deletion[iregion] = cycle_rt->deletion[iregion];
                    if( cycle_rt->soft_clip[iregion] > rt->soft_clip[iregion] )
                        rt->soft_clip[iregion] = cycle_rt->soft_clip[iregion];
                    if( cycle_rt->known_snp[iregion] > rt->known_snp[iregion] )
                        rt->known_snp[iregion] = cycle_rt->known_snp[iregion];
                }
            }
        }

        // for each region convert value to a simple pass/fail flag
        for (iregion=0; iregion<rt->nregions; iregion++)
        {
            // coverage - correct for sparse bins by assuming ALL bins have atleast 10 clusters
            int n = rt->align[iregion] < 10 ? 10 : rt->align[iregion];
            // align - mark bins with no coverage
            rt->align[iregion] = (rt->align[iregion] ? 0 : 1);
            // mismatch - mark bins with maximum mismatch rate > threshold
            rt->mismatch[iregion] = (100.0 * rt->mismatch[iregion] / n) >= REGION_MISMATCH_THRESHOLD;
            // insertion - mark bins with maximum insertion rate > threshold
            rt->insertion[iregion] = (100.0 * rt->insertion[iregion] / n) >= REGION_INSERTION_THRESHOLD;
            // deletion - mark bins with maximum deletion rate > threshold
            rt->deletion[iregion] = (100.0 * rt->deletion[iregion] / n) >= REGION_DELETION_THRESHOLD;
            // soft_clip - mark bins with maximum soft_clip rate > threshold
            rt->soft_clip[iregion] = (100.0 * rt->soft_clip[iregion] / n) >= REGION_SOFT_CLIP_THRESHOLD;
        }

        snprintf(file, file_sz, "%s/BaseCalls/L00%d/s_%d_%d.filter", s->intensity_dir, lane, lane, tile);
        filter = load_filter(s, file, &nclusters);

        iregions = load_positions(s, lane, tile, nclusters);

        // update filter bitmap for ALL cluster in each region
        for (icluster=0; icluster<nclusters; icluster++)
        {
            iregion = iregions[icluster];
            if( rt->align[iregion] )
                filter[icluster] |= FILTER_COVERAGE;
            if( rt->mismatch[iregion] )
                filter[icluster] |= FILTER_MISMATCH;
            if( rt->insertion[iregion] )
                filter[icluster] |= FILTER_INSERTION;
            if( rt->deletion[iregion] )
                filter[icluster] |= FILTER_DELETION;
            if( rt->soft_clip[iregion] )
                filter[icluster] |= FILTER_SOFT_CLIP;
        }
        
#if 1
        size_t filtered = 0, coverage = 0, mismatch = 0, insertion = 0, deletion = 0, soft_clip = 0, exclude = 0;
        for (icluster=0; icluster<nclusters; icluster++)
        {
            if( filter[icluster] & FILTER_PF )
                filtered++;
            else
                continue;
            if( filter[icluster] & FILTER_COVERAGE )
                coverage++;
            if( filter[icluster] & FILTER_MISMATCH )
                mismatch++;
            if( filter[icluster] & FILTER_INSERTION )
                insertion++;
            if( filter[icluster] & FILTER_DELETION )
                deletion++;
            if( filter[icluster] & FILTER_SOFT_CLIP )
                soft_clip++;
            if( filter[icluster] & (FILTER_INSERTION | FILTER_DELETION) )
                exclude++;
        }
        fprintf(stderr, "tile=%-4d\tfiltered=%-8lu coverage=%-8lu mismatch=%-8lu insertions=%-8lu deletions=%-8lu soft_clip=%-8lu exclude=%-8lu\n",
                tile, filtered, coverage, mismatch, insertion, deletion, soft_clip, exclude);
#endif

        output_filter(s, lane, tile, nclusters, filter);

        if( NULL != iregions)
            free(iregions);
        if( NULL != filter)
            free(filter);
    }

    free(file);
    
    return;
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
parse_bam_file_line_full(Settings *s,
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
                         int* read_mismatch,
                         const int read_buff_size) {

    char *name;
    int32_t pos;
    uint32_t *cigar;
    uint8_t *ci_ptr, *seq, *qual, *m_ptr;
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
    read = -1;

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

    if(BAM_FPAIRED & bam->core.flag){
        if(BAM_FREAD1 & bam->core.flag)
            read = 1;
        if(BAM_FREAD2 & bam->core.flag)
            read = 2;
        if(read == -1){
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
    
    j = 0;
    for (i=0; i<bam->core.n_cigar; i++) {
        int l = cigar[i] >> 4, op = cigar[i] & 0xf, k;
        switch(op) {
        case BAM_CMATCH:
            // CIGAR: alignment match;
            for(k=0; k<l; j++, k++)
                read_mismatch[j] |= BASE_ALIGN;
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
                read_mismatch[j] |= BASE_INSERTION;
            break;
        case BAM_CSOFT_CLIP:
            // CIGAR: clip on the read with clipped sequence present in qseq
            for(k=0; k<l; j++, k++)
                read_mismatch[j] |= BASE_SOFT_CLIP;
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
    read = -1;

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
    *bam_offset = offset;
    *bam_read = read;

    return 0;
}

/*
 * cts - parse bam file line
 *
 * returns 0 on success, 1 on expected failure.
 */
int dump_bam_file(Settings *s, samfile_t *fp_bam, size_t *nreads) {

    size_t nreads_bam = 0;

    static const int bam_read_buff_size = 1024;
    char bam_read_seq[bam_read_buff_size];
    int bam_read_qual[bam_read_buff_size];
    int bam_read_mismatch[bam_read_buff_size];

    bam1_t *bam = bam_init1();

    /* loop over reads in the bam file */
    while (1){
        int bam_lane = -1, bam_tile = -1, bam_read = -1, bam_x = -1, bam_y = -1, read_length;
        size_t bam_offset = 0;

        if (0 != parse_bam_file_line_full(s, fp_bam, bam,
                                          &bam_lane, &bam_tile, &bam_x, &bam_y, &bam_offset, &bam_read,
                                          bam_read_seq, bam_read_qual, bam_read_mismatch, bam_read_buff_size)) {
            break;
        }

        read_length = strlen(bam_read_seq);
        if (0 == read_length) continue;

        char *name = bam1_qname(bam);
        uint8_t *oq_ptr;
        char *oq = NULL;
        int i;

        oq_ptr = bam_aux_get(bam, "OQ");
        if(NULL != oq_ptr){
            oq = bam_aux2Z(oq_ptr);
            if(NULL == oq){
                fprintf(stderr,"ERROR: Invalid original qualities %s for read: \"%s\"\n",oq,name);
                exit(EXIT_FAILURE);
            }
        }
        
        printf("%s\t%d\t%s", name, bam->core.flag, bam_read_seq);
        /* stringify quality values and mismatch bitmap into read_seq */
        for(i=0;i<bam->core.l_qseq;i++)
            bam_read_seq[i] = bam_read_qual[i] + PHRED_QUAL_OFFSET;
        printf("\t%s", bam_read_seq);
        for(i=0;i<bam->core.l_qseq;i++)
            bam_read_seq[i] = bam_read_mismatch[i] + PHRED_QUAL_OFFSET;
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
    char *pn = "spatial_filter";
    char *pp = NULL;
    char *ds = "A program to apply a spatial filter";
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
    
    pgsize = 128 + strlen(pn) + strlen(pn) + strlen(ds) + strlen(PBP_VERSION) + strlen(s->cmdline);
    if( NULL != pp ){
        pgsize += strlen(pp);
        pg = smalloc(pgsize);
        pglen = snprintf(pg, pgsize, "@PG\tID:%s\tPN:%s\tPP:%s\tDS:%s\tVN:" PBP_VERSION "\tCL:%s\n", id, pn, pp, ds, s->cmdline);
    }else{
        pg = smalloc(pgsize);
        pglen = snprintf(pg, pgsize, "@PG\tID:%s\tPN:%s\tDS:%s\tVN:" PBP_VERSION "\tCL:%s\n", pn, pn, ds, s->cmdline);
    }
    assert(pglen < pgsize);
    
    append_header_text(bam_header, pg, pglen);

    free(text);
    free(pg);
}

static int updateRegionTable(Settings *s, RegionTable *rts, int read, int x, int y, int *read_mismatch) {

    float x_coord = (x - COORD_SHIFT) / COORD_FACTOR;
    float y_coord = (y - COORD_SHIFT) / COORD_FACTOR;
    int iregion = (int)(x_coord / REDUCTION_FACTOR) * s->nregions_y + (int)(y_coord / REDUCTION_FACTOR);
    int cycle;

    /* update region table */
    for( cycle=0; cycle<s->read_length[read]; cycle++ ){
        RegionTable *rt = rts + cycle;

        if( read_mismatch[cycle] & BASE_INSERTION )
            rt->insertion[iregion]++;
        if( read_mismatch[cycle] & BASE_DELETION )
            rt->deletion[iregion]++;
        if( read_mismatch[cycle] & BASE_SOFT_CLIP )
            rt->soft_clip[iregion]++;
        if( read_mismatch[cycle] & BASE_KNOWN_SNP ) {
            rt->known_snp[iregion]++;
        } else {
            if( read_mismatch[cycle] & BASE_ALIGN )
                rt->align[iregion]++;
            if( read_mismatch[cycle] & BASE_MISMATCH )
                rt->mismatch[iregion]++;
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
int makeRegionTable(Settings *s, samfile_t *fp_bam, RegionTable **rts, int *ntiles, size_t *nreads) {

    int nrt = 0;

    int lane = -1;
    int tile = -1;

    int ntiles_bam = 0;
    size_t nreads_bam = 0;

    static const int bam_read_buff_size = 1024;
    char bam_read_seq[bam_read_buff_size];
    int bam_read_qual[bam_read_buff_size];
    int bam_read_mismatch[bam_read_buff_size];

    bam1_t *bam = bam_init1();

    int itile = -1;

    /* loop over reads in the bam file */
    while (1){
        int bam_lane = -1, bam_tile = -1, bam_read = -1, bam_x = -1, bam_y = -1, read_length;
        size_t bam_offset = 0;
        int irt, cycle;

        if (0 != parse_bam_file_line_full(s, fp_bam, bam,
                                          &bam_lane, &bam_tile, &bam_x, &bam_y, &bam_offset, &bam_read,
                                          bam_read_seq, bam_read_qual, bam_read_mismatch, bam_read_buff_size)) {
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
            
            if (!s->quiet)fprintf(stderr, "Processing tile %i (%lu)\n", tile, nreads_bam);

            for (itile=0; itile<ntiles_bam; itile++)
                if (tile == rts[itile][0].tile) {
                    fprintf(stderr,"ERROR: alignments are not sorted by tile.\n");
                    exit(EXIT_FAILURE);
                }

            ntiles_bam++;
            if(ntiles_bam >= N_TILES){
                fprintf(stderr,"ERROR: too many tiles %d > %d.\n", ntiles_bam, N_TILES);
                exit(EXIT_FAILURE);
            }

            nrt = 2 * N_CYCLES + 1;
            rts[itile] = smalloc(nrt * sizeof(RegionTable));
            for(irt=0;irt<nrt;irt++)
                rts[itile][irt].nregions = 0;
        }

        irt = (bam_read > 1 ? 1 : 0) * N_CYCLES;
        if (0 == rts[itile][irt].nregions)
            for(cycle = 0; cycle < read_length; irt++, cycle++)
                initialiseRegionTable(s, rts[itile] + irt, lane, tile, bam_read, cycle);
        
        irt = (bam_read > 1 ? 1 : 0) * N_CYCLES;
        if (0 != updateRegionTable(s, rts[itile] + irt, bam_read, bam_x, bam_y, bam_read_mismatch)) {
            fprintf(stderr,"ERROR: updating quality values for tile %i.\n", tile);
            exit(EXIT_FAILURE);
        }
        nreads_bam++;
    }
    
    bam_destroy1(bam);
    
    for (itile=0; itile<ntiles_bam; itile++)
        completeRegionTables(s, rts[itile]);

    *ntiles = ntiles_bam;
    *nreads = nreads_bam;

    return nrt;
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
int filter_bam(Settings *s, samfile_t *fp_in_bam, samfile_t *fp_out_bam, size_t *nreads, size_t *nfiltered) {
    char *file = NULL;
    size_t file_sz = strlen(s->working_dir) + (NULL == s->prefix ? 0 : strlen(s->prefix)) + 100;
    size_t nfilter = 0;
    uint8_t *filter = NULL;

    int tiles[N_TILES];

    int lane = -1;
    int tile = -1;

    int ntiles_bam = 0;
    size_t nreads_bam = 0;
    int nfiltered_bam = 0;

    bam1_t *bam = bam_init1();

    file = smalloc(file_sz);

    /* loop over reads in the input bam file */
    while (1){
        int bam_lane = -1, bam_tile = -1, bam_read = -1, read_length;
        size_t bam_offset = 0;

        if (0 != parse_bam_file_line(s, fp_in_bam, bam, &bam_lane, &bam_tile, &bam_offset, &bam_read)) {
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
            int itile;

            tile = bam_tile;
            
            if (!s->quiet) fprintf(stderr, "Processing tile %i (%lu)\n", tile, nreads_bam);

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
            
            if( NULL != filter) free(filter);
            snprintf(file, file_sz, "%s/%s_%04d.filter", s->working_dir, s->prefix, tile);
            filter = load_filter(s, file, &nfilter);
        }

        assert(bam_offset < nfilter);

        if (filter[bam_offset] & FILTER_MASK) {
            nfiltered_bam++;

            if (s->qcfail) {
                bam->core.flag |= BAM_FQCFAIL;
            }else{
                continue;
            }
        }

        if( 0 > samwrite(fp_out_bam, bam)) {
            fprintf(stderr, "Error: writing bam file\n");
            exit(EXIT_FAILURE);
        }

        nreads_bam++;
    }
    
    *nreads = nreads_bam;
    *nfiltered = nfiltered_bam;

    bam_destroy1(bam);
    
    if (NULL != filter) free(filter);

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

    fprintf(usagefp, "spatial_filter v" PBP_VERSION "\n\n");
    fprintf(usagefp, 
            "Usage: spatial_filter [options] bam_file\n"
            ""
            "  calculate or apply a spatial filter\n"
            "");
    fprintf(usagefp, "  options:\n");
    fprintf(usagefp, "    -d         just dump bam file in 'mismatch' format\n");
    fprintf(usagefp, "    -c         calculate filter files\n");
    fprintf(usagefp, "    -a         apply filter files\n");
    fprintf(usagefp, "    -f         mark filtered reads as QCFAIL\n");
    fprintf(usagefp, "                 default do not output filtered read\n");
    fprintf(usagefp, "    -intensity-dir dir\n");
    fprintf(usagefp, "               Intensity directory\n");
    fprintf(usagefp, "                 no default\n");
    fprintf(usagefp, "    -p prefix\n");
    fprintf(usagefp, "               filter file prefix e.g. prefix_1101.filter\n");
    fprintf(usagefp, "                 no default\n");
    fprintf(usagefp, "    -snp_file file\n");
    fprintf(usagefp, "               set of snps to be removed\n");
    fprintf(usagefp, "                 file in Reference Ordered Data (ROD) format\n");
    fprintf(usagefp, "    -o output\n");
    fprintf(usagefp, "               Output bam file name\n");
    fprintf(usagefp, "                 no default will write output to stdout\n");
    fprintf(usagefp, "    -u         do not compress the output bam file\n");
    fprintf(usagefp, "                 default compress\n");
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
    int i, j;
    int bam_compress = 1;
    char *in_bam_file;
    samfile_t *fp_input_bam;
    char out_mode[5] = "wb";
    bam_header_t *out_bam_header = NULL;
    char *out_bam_file = NULL;
    samfile_t *fp_output_bam;
    const char *override_intensity_dir = NULL;

    settings.prefix = NULL;
    settings.quiet = 0;
    settings.dump = 0;
    settings.calculate = 0;
    settings.apply = 0;
    settings.qcfail = 0;
    settings.intensity_dir = NULL;
    settings.snp_file = NULL;
    settings.snp_hash = NULL;
    settings.output = NULL;
    settings.read_length[0] = 0;
    settings.read_length[1] = 0;
    settings.read_length[2] = 0;
    settings.working_dir = NULL;
    settings.nregions_x = 0;
    settings.nregions_y = 0;
    settings.nregions = 0;
    
    settings.cmdline = get_command_line(argc, argv);

    /* Parse args */
    for (i = 1; i < argc && argv[i][0] == '-'; i++) {
	if (!strcmp(argv[i], "-")) {
	    break;

	} else if (!strcmp(argv[i], "-d")) {
	    settings.dump = 1;

	} else if (!strcmp(argv[i], "-c")) {
	    settings.calculate = 1;

	} else if (!strcmp(argv[i], "-a")) {
	    settings.apply = 1;
	} else if (!strcmp(argv[i], "-f")) {
	    settings.qcfail = 1;
	} else if (!strcmp(argv[i], "-u")) {
            bam_compress = 0;

	} else if (!strcmp(argv[i], "-o")) {
            if(settings.output != NULL) {
		fprintf(stderr, "ERROR: -o option specified multiple times\n");
                usage(1);
            }
            check_arg(i,argc,"-o");
            settings.output = argv[++i];

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

	} else if (!strcmp(argv[i], "-p")) {
            if(settings.prefix != NULL) {
		fprintf(stderr, "ERROR: -p option specified multiple times\n");
                usage(1);
            }
            check_arg(i,argc,"-p");
            settings.prefix = argv[++i];

	} else if (!strcmp(argv[i], "-q")) {
	    settings.quiet = 1;
	} else if (!strcmp(argv[i], "-h")) {
	    usage(0);
	} else {
            fprintf(stderr,"ERROR: Unknown option %s\n", argv[i]);
	    usage(1);
	}
    }

    if ((argc-i) < 1)
	usage(0);

    /* preserve starting directory */
    settings.working_dir = alloc_getcwd();
    if (NULL == settings.working_dir) {
        fprintf(stderr, "ERROR: can't obtain working directory: %s\n",
                strerror(errno));
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

    in_bam_file = argv[i++];

    /* dump the alignments */
    if (settings.dump) {
        size_t nreads = 0;

        fp_input_bam = samopen(in_bam_file, "rb", 0);
        if (NULL == fp_input_bam) {
            fprintf(stderr, "ERROR: can't open bam file file %s: %s\n",
                    in_bam_file, strerror(errno));
            exit(EXIT_FAILURE);
        }

        if (0 != dump_bam_file(&settings, fp_input_bam, &nreads)) {
            fprintf(stderr,"ERROR: failed to dump bam file %s\n", in_bam_file);
            exit(EXIT_FAILURE);
        }

        /* close the bam file */
        samclose(fp_input_bam);

        if (!settings.quiet) {
            fprintf(stderr, "Dumped %8lu traces\n", nreads);
        }
    
        if (NULL != settings.working_dir) free(settings.working_dir);

        return EXIT_SUCCESS;
    }
    
    /* get absolute intensity dir */
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

    if (NULL == settings.prefix) {
        fprintf(stderr,"ERROR: you must specify a prefix\n");
        exit(EXIT_FAILURE);
    }

    /* calculate the filter */
    if (settings.calculate) {
        int ntiles = 0;
        size_t nreads = 0;
        int nrt = 0;
        RegionTable **rts = NULL;

        fp_input_bam = samopen(in_bam_file, "rb", 0);
        if (NULL == fp_input_bam) {
            fprintf(stderr, "ERROR: can't open bam file file %s: %s\n",
                    in_bam_file, strerror(errno));
            exit(EXIT_FAILURE);
        }

        /* set the number of regions */
        setRegions(&settings);
        if (0 > settings.nregions) {
            fprintf(stderr, "ERROR: invalid tile size\n");
            exit(EXIT_FAILURE);
        }

        rts = smalloc(N_TILES * sizeof(RegionTable *));

        nrt = makeRegionTable(&settings, fp_input_bam, rts, &ntiles, &nreads);
        if (0 == nrt) {
            fprintf(stderr,"ERROR: failed to make region table\n");
            exit(EXIT_FAILURE);
        }

        /* close the bam file */
        samclose(fp_input_bam);

        if (!settings.quiet) {
            fprintf(stderr, "Processed %8lu traces\n", nreads);
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

        findBadRegions(&settings, ntiles, rts);

        if (NULL != rts) {
            for (i=0; i<ntiles; i++) {
                for (j=0; j<nrt; j++)
                    freeRegionTable(&rts[i][j]);
                free(rts[i]);
            }
            free(rts);
        }
    }
    
    /* apply the  filter */
    if (settings.apply) {
        size_t nreads = 0;
        size_t nfiltered = 0;

        if (0 == bam_compress)
            strcat(out_mode, "u");

        fp_input_bam = samopen(in_bam_file, "rb", 0);
        if (NULL == fp_input_bam) {
            fprintf(stderr, "ERROR: can't open bam file file %s: %s\n",
                    in_bam_file, strerror(errno));
            exit(EXIT_FAILURE);
        }

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

        if (-1 == filter_bam(&settings, fp_input_bam, fp_output_bam, &nreads, &nfiltered)) {
            fprintf(stderr,"ERROR: failed to filter bam file %s\n", in_bam_file);
            exit(EXIT_FAILURE);
        }

        samclose(fp_input_bam);
        samclose(fp_output_bam);

        if (!settings.quiet) {
            fprintf(stderr, "Processed %8lu traces\n", nreads);
            fprintf(stderr, "Filtered %8lu traces\n", nfiltered);
        }

        /* back to where we belong */
        checked_chdir(settings.working_dir);
    }
    
    if (NULL != settings.working_dir) free(settings.working_dir);

    return EXIT_SUCCESS;

}
