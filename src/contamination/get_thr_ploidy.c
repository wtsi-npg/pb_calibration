/*  File: get_thr_ploidy_vcfa_ch.c // filters 1 1 vcfq file and computes av std Depth  depending on Depth of vcfq file

 * Authors: designed by Irina Abnizova (ia1) and edited by Steve Leonard(srl)
 *
  Last edited:

  14 April- merge with Steve's changes

  March 2016- prepared to read chr names and letters; skip bad lines
  18 Jan 2016 - real std
  11 June: one thr file. output is e.g.  2 2 41 113
  9 June 2014 :n1 from pipe of three:filters 1 1 vcfq file and computes av std Depth  depending on Depth of vcfq file
  28 May : not output filtered vcfq
  22 May  get_filters 1 1 .vcfq  RA_ms.filters
  21 May  filters bad varians: abnormal Depth
  30 April: adds Filteref 1,1 vcfq + distr_1_1
  29 April
  // computes suitable thresholds bases on avDP4=actual (not cov before as earlier!!!

  3 Feb-put thrR thrA into file thrD ;
  - get autom thr, starting with default ones 1,5   still old fields of vcf input

  *-------------------------------------------------------------------
 * Description: computes computes threshold file depending on Depth of vcfq file, and some extra outputs

 * Exported functions:
 * HISTORY:

*/

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <getopt.h>

// ******** predefined user constants,
#define NPAR 3                // Number of arguments
#define THR 1                 // Default reference allele threshold
#define THA 1                 // Default alternative allele threshold
#define NRID 6                // Number of fields (pos,Depth,DP4) in extracted_vcf files
#define MIN_COUNT_AFTERF 200  // Minimum variant count after filtering

// ******** declarations of  functions
void usage(int code); // usage
int GetAf (int []); // computes AF percent for a variant position
int GetMax (int [], int i1, int i2); // percent giving max AAF
int ploidy (int ind); // defines ploidy from results of GetMax ploid=2 if ~50%

////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////
int main (int argc, char *argv[]) {

    // min_depths for Ref and Alternative alleles
    int thR = THR;
    int thA = THA;

    // verbose mode
    int verbose = 0;

    // file names
    char *extract_vcf_file, *thr_file, *filtered_vcf_file;

    // file handles
    FILE *extract_vcfFile, *thrFile, *filtered_vcfFile;

    // to read vcfq
    static const int line_size = 8192; // maximum line size
    char line[line_size];

    // values in vcfq file
    int DP4[4]; // 4xdepths
    int pos, D; // genome pos of error, depth

    int count_beforeF; // variant count before filtering
    int count_afterF; // variant count after filtering

    // arrays
    int histAF[100]; // maf percentages

    // stats
    float sumDP4, sumDP4_2; // sum of depths and sum of square depths used for mu/std calculation

    // to compute for outputs
    int ploid; // polidy
    int avDP4, stdDP4; // mean and std of depth after filtering
    int avDP41; // adjusted mean depth
    int thRR, thAA; // estimated thresholds based on actual depths avDP4 and stdDP4

    int i; // array index

    static struct option long_options[] = 
        { {"thR", 1, 0, 'r'},
          {"thA", 1, 0, 'a'},
          {"verbose", 0, 0, 'v'},
          {"help", 0, 0, 'h'},
          {0, 0, 0, 0}
        };

    char c;
    while ( (c = getopt_long(argc, argv, "r:a:vh?", long_options, 0)) != -1) {
        switch (c) {
            case 'r': thR = atoi(optarg); break;
            case 'a': thA = atoi(optarg); break;
            case 'v': verbose = 1; break;
            case 'h':
            case '?': usage(0); break;
            default: fprintf(stderr, "ERROR: Unknown option %c\n", c);
                     usage(1);
                     break;
        }
    }

    if ( (argc-optind) < NPAR) {
        // not enough parameters
        usage(-1);
    } else {
        extract_vcf_file  = argv[optind+0];
        thr_file          = argv[optind+1];
        filtered_vcf_file = argv[optind+2];
    }

    // open files for read/write
    extract_vcfFile = fopen(extract_vcf_file,"r");
    if (extract_vcfFile == NULL) {
        fprintf(stderr, "cannot open input_vcf_file %s: %s\n", extract_vcf_file, strerror(errno));
        exit(EXIT_FAILURE);
    }
    thrFile = fopen(thr_file,"w");
    if (thrFile == NULL) {
        fprintf(stderr, "cannot open threshold_file %s: %s\n", thr_file, strerror(errno));
        exit(EXIT_FAILURE);
    }
    filtered_vcfFile = fopen(filtered_vcf_file,"w");
    if (filtered_vcfFile == NULL) {
        fprintf(stderr, "cannot open filtered_vcf_file %s: %s\n", filtered_vcf_file, strerror(errno));
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "get_thr_ploidy\n");

    // initialise AAF histogram
    for (i=0;i<100;i++) {
        histAF[i] = 0;
    }
    
    // initialise counts
    count_beforeF = 0;
    count_afterF = 0;

    // initialise sums
    sumDP4 = 0;
    sumDP4_2 = 0;

    // read vcfq file - 6 fields position, depth and 4xdepths (ref/forward, ref/reverse, alt/forward and alt/reverse) 
    while (fgets(line, line_size, extract_vcfFile)) {
        int k = sscanf(line, "%d,%d,%d,%d,%d,%d", &pos, &D, &DP4[0], &DP4[1], &DP4[2], &DP4[3]);
        if (k != NRID) {
            // number of fields read not correct
            fprintf(stderr, "skipping malformed VCF line %s\n", line);
        }
        count_beforeF++;

        // filter for DP4 separately for ref and alternative alleles
        if (DP4[0] >= thR && DP4[1]>= thR && DP4[2]>= thA && DP4[3]>= thA) {
            int AF;
            float sDP4;

            count_afterF++;

            // compute AF and update histo
            AF = GetAf(DP4);
            histAF[AF-1]++;

            // calc sDP4 and update sumDP4, sumDP4_2
            sDP4 = DP4[0] + DP4[1] + DP4[2] + DP4[3];
            sumDP4 +=  sDP4;
            sumDP4_2 +=  sDP4*sDP4;

            // write to filtered vcf file
            if (fprintf(filtered_vcfFile,"%d,%d,%d,%d,%d,%d\n",pos,D, DP4[0], DP4[1],DP4[2], DP4[3]) <= 0) {
                fprintf(stderr, "error writing filteredVCF file: %s\n", strerror(errno));
                exit(EXIT_FAILURE);
            }
        }
    }

    fprintf(stderr, "count of variants before filtering = %d\n", count_beforeF);
    fprintf(stderr, "count of variants after filtering = %d\n", count_afterF);

    // initialise values
    avDP4 = 0; // mean depth
    stdDP4 = 0; // stddev depth
    ploid = 1; // ploidy
    avDP41 = 0; // adjustered mean depth
    thRR = 0; // estimated reference threshold
    thAA = 0; // estimated alternative threshold

    // do we have enough data ?
    if (count_afterF < MIN_COUNT_AFTERF) {
      fprintf(stderr, "Not enough variants after filtering, count_afterF < %d\n", MIN_COUNT_AFTERF);

    } else {
        int ind, ind50;
        float perc_left;
        float ratio;

        // ploidy detection
        ind = GetMax(histAF, 1, 80);// percent giving max peak in hist in interval (1,80)%
        fprintf(stderr, "percent of max AAF= %d\n", ind);

        ind50 = GetMax(histAF, 40, 60);//percent giving max peak in hist in interval (40,60)%
        fprintf(stderr, "percent of max AAF around 50= %d\n", ind50);//

        ploid = ploidy(ind);// automated detection of ploidy: 1 or 2
        fprintf(stderr, "ploidy= %d\n", ploid);

        // compute stats: average depth : default and actual after 1,1 and q25 filtering

        avDP4 = ceil(sumDP4 / count_afterF); // actual average cov  after Filt
        stdDP4 = ceil(sqrt((sumDP4_2 / count_afterF) - avDP4 * avDP4));//
        perc_left = (100.0 * count_afterF / count_beforeF);

        fprintf(stderr, "average depth after filtering= %d\n", avDP4);
        fprintf(stderr, "stdev depth after filtering= %d\n", stdDP4);
        fprintf(stderr, "percent data left filtering= %.2f\n", perc_left);

        // adjust depth to account for bad (very deep) regions

        avDP41 = avDP4;
        ratio = sqrt((sumDP4_2 / count_afterF) - (avDP4 * avDP4)) / (sumDP4 / count_afterF); // stdDP4 / avDP4;
        fprintf(stderr, "cv actual depth = %.2f\n", ratio);

        if (ratio > 1.5) {
            // bad regions
            avDP41 = ceil(0.8*avDP4);
            fprintf(stderr, "moderately deep regions detected\n");
            fprintf(stderr, "adjusted average depth = %d\n", avDP41);
        }

        if (ratio > 3) {
            // very bad regions
            avDP41 = ceil(0.7*avDP4);
            fprintf(stderr, "extremely deep covered detected\n");
            fprintf(stderr, "adjusted average depth = %d\n", avDP41);
        }

        // calculate thresholds using adjusted Depth
        if (avDP41 > 0 && avDP41 < 4 )
        {
            thRR = 100; thAA = 100; // stupidly high ?
            fprintf(stderr, "coverage too low, adjusted average depth= %d\n", avDP41);
            fprintf(stderr, "theshold set to exclude all variants in mixture calculation\n");
        }
        if (avDP41 >= 4 && avDP41 < 10) {
            thRR = 1; thAA = 1;
        }
        if (avDP41 >= 10 && avDP41 < 70) {
            thRR = 2; thAA = 2;
        }
        if (avDP41 >= 70 && avDP41 < 90) {
            thRR = 3; thAA = 2;
        }
        if (avDP41 >= 90) {
            thRR = 3; thAA = 3;
        }
    }

    // write threshold file
    if (fprintf(thrFile, "%d %d %d %d %d\n", thRR, thAA, avDP41, stdDP4, ploid) <= 0) {
        fprintf(stderr, "error writing threshold_file: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    
    // close files
    fclose(extract_vcfFile);
    fclose(thrFile);
    fclose(filtered_vcfFile);

    fprintf(stderr, "done get_thr_ploidy.\n");
    return 0;
}

////////////////////////////////////////////////////
// usage
////////////////////////////////////////////////////
void usage(int code)
{
    FILE *usagefp = stderr;

	fprintf(usagefp, "get_thr_ploidy\n\n");
	fprintf(usagefp,
		"Usage: get_thr_ploidy [options] input_vcf_file threshold_file filtered_vcf_file\n"
		"\n" "  computes minimum depth thresholds for reference and alternative alleles\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "  options:\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "    -r  minimum reference allele depth\n");
	fprintf(usagefp, "        default 1\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "    -a  minimum reference allele depth\n");
	fprintf(usagefp, "        default 1\n");
	fprintf(usagefp, "\n");
	fprintf(usagefp, "    -v  verbose\n");
	fprintf(usagefp, "        default false\n");
	fprintf(usagefp, "\n");

	exit(code);
}

////////////////////////////////////////////////////
// Calculates AF=alternative allele frequency from data=DP4 counts
// input - array (4-vector) of DP4 counts
// output -index into array of percentages 0,..,100
////////////////////////////////////////////////////
int GetAf (int data[])
{
    float sum, AF1;
    int i, AF;

    sum = data[0];
    for (i=1; i<4; i++) {
        sum += data[i];
    }
    AF1 = (data[2] + data[3]) / sum;
    AF = (int) (100*AF1 + 0.5);

    return AF;
}

////////////////////////////////////////////////////
// find position of the max value in a sub-interval (i1<=i<=i2) of an array
////////////////////////////////////////////////////
int GetMax (int data[], int i1, int i2)
{
  int Imax;
  int ind, i;

  Imax = data[i1];
  ind = i1;
  for (i=i1; i<i2; i++) {
      if( Imax < data[i]) {
          Imax = data[i];
          ind = i;
      }
  }

  return ind;
}

////////////////////////////////////////////////////
// identify ploidy
// input - percentage corresponding to max AAF
// output - 2 if max between 40% and 60% otherwise 1
////////////////////////////////////////////////////
int ploidy (int ind)
{
    int ploid = 1;

    if ( ind >40 && ind < 60) {
        ploid = 2;
    }

    return ploid;
}
