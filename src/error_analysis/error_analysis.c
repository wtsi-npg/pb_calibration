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
 */

#define QC_FAIL
#define PROPERLY_PAIRED

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
#include <faidx.h>

#include <smalloc.h>
#include <aprintf.h>
#include <die.h>
#include <rts.h>
#include <shared.h>
#include <snp.h>
#include <parse_bam.h>
#include <cif.h>

#include <version.h>

typedef struct {
    char *cmdline;
    char *snp_file;
    HashTable *snp_hash;
    char *ref_file;
    faidx_t *fai;
    int depth;
    int maf;
    int mad;
    int qmin;
    int quiet;
} Settings;

#define MAX_DEPTH 50

typedef struct {  
    Settings *s;
    size_t npos;
    size_t nbases;
    int ref_tid;
    int ref_len;
    char *ref;
    size_t depth[MAX_DEPTH+1];
    size_t freq[101];
    samfile_t *in;  
} tmpstruct_t;  
  
// callback for sampileup
static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)  
{  
    tmpstruct_t *tmp = (tmpstruct_t*)data;  
    int maf = -1, mad = -1, known_snp = 0, i;
    
    if (NULL != tmp->s->fai && tid != tmp->ref_tid) {
        if( NULL != tmp->ref ) {
            free(tmp->ref);
            tmp->ref = NULL;
        }
        if (!tmp->s->quiet) fprintf(stderr, "fetching target %s\n", tmp->in->header->target_name[tid]);
        tmp->ref = faidx_fetch_seq(tmp->s->fai, tmp->in->header->target_name[tid], 0, 0x7fffffff, &tmp->ref_len);
        if (NULL == tmp->ref) {
            fprintf(stderr, "ERROR: reading fetching target %s\n", tmp->in->header->target_name[tid]);
            exit(EXIT_FAILURE);
        }
        tmp->ref_tid = tid;
    }

    int rb = (NULL != tmp->ref && pos < tmp->ref_len) ? tmp->ref[pos] : 'N';
    int depth = 0;
    int ncalls[16];
    memset(ncalls, 0, 16 * sizeof(int));
    for (i=0; i<n; i++) {
        if (pl[i].is_del) continue;
        int c = bam_nt16_rev_table[bam1_seqi(bam1_seq(pl[i].b), pl[i].qpos)];
        uint8_t *q = bam1_qual(pl[i].b);
        if (q[pl[i].qpos] < tmp->s->qmin) continue;
        if (c == '=') {
            /* if calls that match the reference are a mixture of [=ACGTN] they may be counted separately if we have no reference */
            ncalls[bam_nt16_table[rb]]++;
        }else{
            ncalls[bam_nt16_table[c]]++;
        }
        depth++;
    }
    qsort(ncalls, 16, sizeof(int), int_sort);

    if (depth < tmp->s->depth) return 0;

    tmp->npos++;
    tmp->depth[depth > MAX_DEPTH ? MAX_DEPTH : depth]++;

    /* frequency of the minor allele rounded to an integer */
    maf = 100.0 * ncalls[(rb == 'N' ? 14 : 15)] / depth;

    /* depth of the minor allele */
    mad = ncalls[(rb == 'N' ? 14 : 15)];

    if (NULL != tmp->s->snp_hash && ncalls[14]) {
        char *chrom = tmp->in->header->target_name[tid];
        char key[100];
        HashItem *hi;
        snprintf(key, sizeof(key), "%s:%d", chrom, pos);
        if (NULL != (hi = HashTableSearch(tmp->s->snp_hash, key, strlen(key)))) {
            hi->data.i = n;
            known_snp = 1;
        }
    }
    if (!known_snp && maf >= tmp->s->maf && mad >= tmp->s->mad) tmp->freq[maf]++;
    
    tmp->nbases += depth;

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
int runPileup(Settings *s, samfile_t *fp_bam, size_t *npos, size_t *nbases) {

    int ret = 0;
    int depth, maf;

    tmpstruct_t tmp;  
    tmp.s = s;
    tmp.npos = 0;
    tmp.nbases = 0;
    tmp.ref_tid = -1;
    tmp.ref_len = -1;
    tmp.ref = NULL;
    memset(tmp.depth, 0, (MAX_DEPTH+1) * sizeof(size_t));
    memset(tmp.freq, 0, 101 * sizeof(size_t));
    tmp.in = fp_bam;

    sampileup(tmp.in, -1, pileup_func, &tmp);

    for(depth=0; depth<(MAX_DEPTH+1); depth++) {
        if( 0 == tmp.depth[depth] ) continue;
        fprintf(stderr,"depth=%d count=%lu\n", depth, tmp.depth[depth]);
    }
    for(maf=0; maf<101; maf++) {
        if( 0 == tmp.freq[maf] ) continue;
        fprintf(stderr,"maf=%d count=%lu\n", maf, tmp.freq[maf]);
    }
            
    *npos = tmp.npos;
    *nbases = tmp.nbases;

    return ret;
}

static
void usage(int code) {
    FILE* usagefp = stderr;

    fprintf(usagefp, "pb_calibration %s\n\n", version);
    fprintf(usagefp, 
            "Usage: error_analysis [options] bam_file\n"
            ""
            "  outputs a histogram of allele frequencies\n"
            "");
    fprintf(usagefp, "  options:\n");
    fprintf(usagefp, "    -snp_file file\n");
    fprintf(usagefp, "               set of snps to be removed be calibration\n");
    fprintf(usagefp, "                 file in Reference Ordered Data (ROD) format\n");
    fprintf(usagefp, "    -ref_file file\n");
    fprintf(usagefp, "               reference fasta file\n");
    fprintf(usagefp, "    -depth depth\n");
    fprintf(usagefp, "               minimum depth\n");
    fprintf(usagefp, "    -maf minor allele frequency\n");
    fprintf(usagefp, "               minimum minor allele frequenecy\n");
    fprintf(usagefp, "    -mad minor allele depth\n");
    fprintf(usagefp, "               minimum minor allele depth\n");
    fprintf(usagefp, "    -qmin qual\n");
    fprintf(usagefp, "               minimum quality\n");
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
    size_t npos = 0;
    size_t nbases = 0;
    int ret;

    settings.quiet = 0;
    settings.depth = 4;
    settings.maf = 5;
    settings.mad = 1;
    settings.qmin = 30;
    settings.snp_file = NULL;
    settings.snp_hash = NULL;
    settings.ref_file = NULL;
    settings.fai = NULL;
    
    settings.cmdline = get_command_line(argc, argv);

    /* Parse args */
    for (i = 1; i < argc && argv[i][0] == '-'; i++) {
	if (!strcmp(argv[i], "-")) {
	    break;

	} else if (!strcmp(argv[i], "-snp_file")) {
            if(settings.snp_file != NULL) {
		fprintf(stderr, "ERROR: -snp_file specified multiple times\n");
                usage(1);
            }
            check_arg(i,argc,"-snp_file");
            settings.snp_file = argv[++i];

	} else if (!strcmp(argv[i], "-ref_file")) {
            if(settings.ref_file != NULL) {
		fprintf(stderr, "ERROR: -ref_file specified multiple times\n");
                usage(1);
            }
            check_arg(i,argc,"-ref_file");
            settings.ref_file = argv[++i];

	} else if (!strcmp(argv[i], "-depth")) {
            check_arg(i,argc,"-depth");
	    settings.depth = atoi(argv[++i]);
	} else if (!strcmp(argv[i], "-maf")) {
            check_arg(i,argc,"-maf");
	    settings.maf = atoi(argv[++i]);
	} else if (!strcmp(argv[i], "-mad")) {
            check_arg(i,argc,"-mad");
	    settings.mad = atoi(argv[++i]);
	} else if (!strcmp(argv[i], "-qmin")) {
            check_arg(i,argc,"-qmin");
	    settings.qmin = atoi(argv[++i]);

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

    /* read the snp_file */
    if (NULL != settings.snp_file) {
        settings.snp_hash = readSnpFile(settings.snp_file);
        if (NULL == settings.snp_hash) {
            fprintf(stderr, "ERROR: reading snp file %s\n", settings.snp_file);
            exit(EXIT_FAILURE);
        }
    }

    /* read the reference */
    if (NULL != settings.ref_file) {
        settings.fai = fai_load(settings.ref_file);
        if (0 == settings.fai) {
            fprintf(stderr, "ERROR: reading reference file %s\n", settings.ref_file);
            exit(EXIT_FAILURE);
        }
    }

    /* open the bam file */
    bam_file = argv[i++];
    fp_bam = samopen(bam_file, "rb", 0);
    if (NULL == fp_bam) {
        fprintf(stderr, "ERROR: can't open bam file file %s: %s\n",
                bam_file, strerror(errno));
        exit(EXIT_FAILURE);
    }

    /* run pileup */
    ret = runPileup(&settings, fp_bam, &npos, &nbases);
    if (0 > ret) {
        fprintf(stderr,"Error: running pileup\n");
        exit(EXIT_SUCCESS);
    }

    if (!settings.quiet) {
        fprintf(stderr, "Processed %lu positions\n", npos);
        fprintf(stderr, "Processed %lu bases\n", nbases);
        if (NULL != settings.snp_hash) {
            HashIter *iter = HashTableIterCreate();
            HashItem *hashItem;
            size_t nsnps = 0;
            while ((hashItem = HashTableIterNext(settings.snp_hash, iter)))
                nsnps += hashItem->data.i;
            fprintf(stderr, "Ignored %lu snps\n", nsnps);
        }
    }

    /* close the bam file */
    samclose(fp_bam);


    return EXIT_SUCCESS;

}
