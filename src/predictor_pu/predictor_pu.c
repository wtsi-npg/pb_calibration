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
 * This code applies a purity or quality based calibration table
 *
 */

#ifdef HAVE_CONFIG_H
#include "pb_config.h"
#endif

#ifdef HAVE_PREAD
# define _XOPEN_SOURCE 500 // for pread
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
#include <shared.h>
#include <parse_bam.h>
#include <cif.h>

#include <version.h>

#define CT_MODE_PURITY  (1<<0)
#define CT_MODE_QUALITY (1<<1)

#define PHRED_QUAL_OFFSET 33  // phred quality values offset

typedef struct {
    int         mode;
    int         tile;
    int         read;
    int         cycle;
    float       offset;
    float       delta;
    float       scale;
    int         nbins;
    float       *predictor;
    long        *num_bases;
    long        *num_errors;
    float       *frac_bases;
    float       *error_rate;
    float       *quality;
    int         *ibin;
} CalTable;

typedef struct {
    char *cmdline;
    char *intensity_dir;
    char *output;
    int read_length[3];
    int cstart[3];
    int quiet;
    char *working_dir;
} Settings;

static void initialiseCalTable(CalTable *ct, int mode, int tile, int read, int cycle)
{
    int i;

    ct->mode = mode;

    ct->tile = tile;
    ct->read = read;
    ct->cycle = cycle;

    if (ct->mode == CT_MODE_PURITY) {
        ct->nbins = 76;
        ct->offset = 0.25;
        ct->delta  = 0.01;
        ct->scale  = 1.0 / ct->delta;
    } else {
        ct->nbins = 51;
        ct->offset = 0.0;
        ct->delta  = 1.0;
        ct->scale  = 1.0 / ct->delta;
    }

    ct->predictor = smalloc(ct->nbins * sizeof(float));
    ct->num_bases = smalloc(ct->nbins * sizeof(long));
    ct->num_errors = smalloc(ct->nbins * sizeof(long));
    ct->frac_bases = smalloc(ct->nbins * sizeof(float));
    ct->error_rate = smalloc(ct->nbins * sizeof(float));
    ct->quality = smalloc(ct->nbins * sizeof(float));

    ct->ibin = smalloc(ct->nbins * sizeof(int));
    for(i=0;i<ct->nbins;i++)
        ct->ibin[i] = -1;

    return;
}

static void freeCalTable(CalTable *ct)
{
    if(ct->nbins) {
        free(ct->predictor);
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
    int mode = (NULL == s->intensity_dir ? CT_MODE_QUALITY : CT_MODE_PURITY);
    CalTable *current_ct = NULL;
    float pmin = 1.0, pmax = 0.0;
    int nbins = -1;

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
        float predictor, frac_bases, error_rate, quality;
        char key[100];
        HashItem *hi;
        HashData hd;
        int ibin;
        
        if( 9 != sscanf(line, "%f %d %d %d %d %d %f %f %f", &predictor, &read, &cycle, &tile,
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

        snprintf(key, sizeof(key), "%d:%d:%d", tile, read, cycle);
        if( NULL == (hi = HashTableSearch(ct_hash, key, strlen(key))) ){
            hd.p = smalloc(sizeof(CalTable));
            if( NULL == HashTableAdd(ct_hash, key, strlen(key), hd, NULL) ) {
                fprintf(stderr, "ERROR: building tile cycle hash table\n");
                exit(EXIT_FAILURE);
            }

            if( tile > 0 ) fprintf(stderr, "bad tile=%4d read=%d cycle=%3d\n", tile, read, cycle);
            
            /* current ct is the new ct */
            current_ct = (CalTable *)hd.p;

            initialiseCalTable(current_ct, mode, tile, read, cycle);
            pmin = current_ct->offset;
            pmax = pmin + (current_ct->nbins - 1 ) * current_ct->delta;

            nbins = current_ct->nbins;
            current_ct->nbins = 0;

            nct++;
        }else{
            current_ct = (CalTable *)hi->data.p;
        }

        if( current_ct->nbins == nbins ){
            fprintf(stderr, "ERROR: number of lines in CT file (%s) for tile=%d read=%d cycle=%d exceeds maximum %d\n",
                    calibrationFile, current_ct->tile, current_ct->read, current_ct->cycle, nbins);
            exit(EXIT_FAILURE);
        }

        if( predictor < pmin || predictor > pmax ) {
            fprintf(stderr, "ERROR: predictor %f outside expected range for %s [%f,%f] in CT file (%s)\n",
                    predictor, (current_ct->mode == CT_MODE_PURITY ? "purity" : "quality"), pmin, pmax, calibrationFile);
            exit(EXIT_FAILURE);
        }

        ibin = current_ct->nbins;

        current_ct->predictor[ibin] = predictor;

        current_ct->num_bases[ibin] = num_bases;
        current_ct->num_errors[ibin] = num_errors;
        current_ct->frac_bases[ibin] = frac_bases;
        current_ct->error_rate[ibin] = error_rate;
        current_ct->quality[ibin] = quality;

        current_ct->nbins++;
    }

    fclose(fp);

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

int Get_bin_predictor(CalTable *ct, int value)
{
    int ibin = 0;

    if (ct->ibin[value] < 0) {
        float predictor = ct->offset + ct->delta * value;
        for(ibin=0;ibin<ct->nbins;ibin++)
            if (predictor<=ct->predictor[ibin])
                break;
        if (ibin >= ct->nbins)
            ibin = (ct->nbins-1);
        ct->ibin[value] = ibin;
    }

    ibin = ct->ibin[value];
    return ibin;
}



/*
  Calculate predictor, update quality value and write bam file

  Returns 0 on success
          1 on EOF
  Quits the program if the data in the files is inconsistent.
*/

static int update_bam_qualities(Settings *s, CalTable **cycle_cts, CifData *cif_data,
                                size_t spot_num, int tile, int x, int y, int read,
                                samfile_t *fp_bam, bam1_t *bam) {
    int cstart = s->cstart[read];
    int read_length = s->read_length[read];
    int c, b;
    uint8_t *seq, *qual;
    static const int read_buff_size = 1024;
    uint8_t oq[read_buff_size];
    int read_qual[read_buff_size];

    assert(read_buff_size >= read_length);
    if (NULL != cif_data) assert((cstart + read_length) <= cif_data->ncycles);

    /* copy original qualities to OQ:Z tag, N.B. have to null terminate aux strings */
    qual = bam1_qual(bam);
    for (b = 0; b < read_length; b++)
        oq[b] = qual[b] + PHRED_QUAL_OFFSET;
    oq[b] = 0;
    bam_aux_append(bam, "OQ", 'Z', read_length+1, oq);

    for (b = 0; b < read_length; b++)
        read_qual[b] = qual[b];
    if (BAM_FREVERSE & bam->core.flag)
        reverse_int(read_qual, read_length);

    /* calculate predictor and update qualities */
    for (c = cstart, b = 0; b < read_length; c++, b++) {
        CalTable *ct = cycle_cts[b];
        float predictor = -1;
        int value, ibin, quality;

        if (ct->mode == CT_MODE_PURITY) {
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
        } else {
            predictor = read_qual[b];
        }

        value = ct->scale * (predictor - ct->offset) + 0.5;

        ibin = Get_bin_predictor(ct, value);
        quality = (int)(ct->quality[ibin] + 0.5);
        read_qual[b] = quality;
    }

    if (BAM_FREVERSE & bam->core.flag)
        reverse_int(read_qual, read_length);

    qual = bam1_qual(bam);
    for (b = 0; b < read_length; b++)
        qual[b] = read_qual[b];

    /* if base is N set the quality to 0 otherwise set the quality to 1 if it is 0 */
    seq = bam1_seq(bam);
    for (b = 0; b < read_length; b++) {
        char read_seq = bam_nt16_rev_table[bam1_seqi(seq, b)];
        if (read_seq == 'N')
            qual[b] = 0;
        else if (0 == qual[b])
            qual[b] = 1;
    }

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
int recalibrate_bam(Settings *s, HashTable *ct_hash, samfile_t *fp_in_bam, samfile_t *fp_out_bam, size_t *bam_nreads) {

    HashTable *tile_hash;

    CifData *cif_data = NULL;

    CalTable ***cycle_cts = NULL;

    int ncycles_firecrest = -1;

    int lane = -1;
    int tile = -1;

    size_t nreads = 0;

    bam1_t *bam = bam_init1();

    if (NULL != s->intensity_dir) checked_chdir(s->intensity_dir);

    tile_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE | HASH_FUNC_JENKINS3);
    if (!tile_hash) die("Failed to create tile_hash\n");

    /* loop over reads in the input bam file */
    while (1){
        int bam_lane = -1, bam_tile = -1, bam_x = -1, bam_y = -1, bam_read = -1, read_length;
        size_t bam_offset = 0;

        if (parse_bam_readinfo(fp_in_bam, bam, &bam_lane, &bam_tile, &bam_x, &bam_y, &bam_read, (NULL == s->intensity_dir ? NULL : &bam_offset))) {
            break;	/* break on end of BAM file */
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

        if (bam_tile != tile) {
            tile = bam_tile;
            
            // lookup itile from tile in tile hash
            HashItem *tileItem = HashTableSearch(tile_hash, (char *)&tile, sizeof(tile));
            if (NULL == tileItem) {
  	        HashData hd;
                int read;
                hd.p = smalloc(N_READS * sizeof(CalTable **));
                if( NULL == HashTableAdd(tile_hash, (char *)&tile, sizeof(tile), hd, NULL) ){
                    fprintf(stderr, "ERROR: building tile hash table\n");
                    exit(EXIT_FAILURE);
                }
                cycle_cts = (CalTable ***)hd.p;
                for(read=0;read<N_READS;read++)
                    cycle_cts[read] = NULL;
                if (!s->quiet) fprintf(stderr, "Processing tile %i (%lu)\n", tile, nreads);
            } else {
                if( NULL != s->intensity_dir ) {
                    fprintf(stderr,"ERROR: alignments are not sorted by tile.\n");
                    exit(EXIT_FAILURE);
                }
                cycle_cts = (CalTable ***)tileItem->data.p;
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
        }

        /* reset cycle ct */
        if (NULL == cycle_cts[bam_read]) {
            int cycle;
            cycle_cts[bam_read] = smalloc(read_length * sizeof(CalTable *));
            for(cycle=0;cycle<read_length;cycle++) {
                char key[100];
                HashItem *hi;
                snprintf(key, sizeof(key), "%d:%d:%d", bam_tile, bam_read, cycle);
                if (NULL == (hi = HashTableSearch(ct_hash, key, strlen(key)))) {
                    snprintf(key, sizeof(key), "%d:%d:%d", -1, bam_read, cycle);
                    if (NULL == (hi = HashTableSearch(ct_hash, key, strlen(key)))) {
                        fprintf(stderr,"ERROR: no calibration table for tile=%d read=%d cycle=%d.\n", bam_tile, bam_read, cycle);
                        exit(EXIT_FAILURE);
                    }
                }
                cycle_cts[bam_read][cycle] = (CalTable *)hi->data.p;
            }
        }

        if (0 != update_bam_qualities(s, cycle_cts[bam_read], cif_data,
                                      bam_offset, bam_tile, bam_x, bam_y, bam_read,
                                      fp_out_bam, bam)) {
            fprintf(stderr,"ERROR: updating quality values.\n");
            exit(EXIT_FAILURE);
        }

        nreads++;
    }
    
    *bam_nreads = nreads;

    bam_destroy1(bam);
    
    if (NULL != tile_hash) {
	HashIter *iter = HashTableIterCreate();
	HashItem *tileItem;
	while ((tileItem = HashTableIterNext(tile_hash, iter))) {
            int read;
            cycle_cts = (CalTable ***)tileItem->data.p;
            for(read=0;read<N_READS;read++)
                if (NULL != cycle_cts[read]) free(cycle_cts[read]);
            free(cycle_cts);
        }
        HashTableIterDestroy(iter);
        HashTableDestroy(tile_hash, 0);
    }
    if (NULL != cif_data) free_cif_data(cif_data);

    return 0;
}


static
void usage(int code) {
    FILE* usagefp = stderr;

    fprintf(usagefp, "pb_predictor %s\n\n", version);
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
    settings.cstart[0] = 0;
    settings.cstart[1] = 0;
    settings.cstart[2] = 0;
    settings.intensity_dir = NULL;
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
    if (override_intensity_dir) {
        settings.intensity_dir = get_real_path_name(override_intensity_dir);
        if (NULL == settings.intensity_dir) {
            fprintf(stderr, "ERROR: can't process intensity dir: %s\n",
                    override_intensity_dir);
            exit(EXIT_FAILURE);
        }
        fprintf(stderr,"Assuming purity cycle calibration table\n");
    } else {
        fprintf(stderr,"Assuming quality cycle calibration table\n");
    }

    if (NULL == (ct_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE|HASH_FUNC_JENKINS3))) {
        fprintf(stderr, "ERROR: creating caltable hash\n");
        exit(EXIT_FAILURE);
    }

    // read the callibration table
    nct = restoreCalTable(&settings, ct_filename, ct_hash);
    if( 0 == nct ){
        fprintf(stderr, "ERROR: no rows in calibration table %s\n", ct_filename);
        exit(EXIT_FAILURE);
    }

    /* Look for CIF directories */
    if (NULL != settings.intensity_dir) {
        get_cif_dirs(settings.intensity_dir);
    }

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
    bam_header_add_pg("pb_cal", "predictor_pu", "A program to apply a calibration table", settings.cmdline, out_bam_header);

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
	HashIter *iter = HashTableIterCreate();
	HashItem *hashItem;
	while ((hashItem = HashTableIterNext(ct_hash, iter)))
            freeCalTable((CalTable *)hashItem->data.p);
        HashTableIterDestroy(iter);
        HashTableDestroy(ct_hash, 0);
    }

    if (NULL != settings.working_dir) free(settings.working_dir);
    if (NULL != settings.cmdline) free(settings.cmdline);
    if (NULL != settings.intensity_dir) free(settings.intensity_dir);

    return EXIT_SUCCESS;

}
