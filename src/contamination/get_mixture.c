 /*  File: get_mixture.c // computes and analyses AAF( maf) distribution over variant positions from vcfq files
 * Authors: designed by Irina Abnizova (ia1)edited by Steve Leonard (srl)
 *

  Last edited:
  13 April- added Steven's corrections and outputs
  April - if empty vcfq file; if small data in histo thr_tot=200 variants
  14 March  - three modes and likelihoods: low, middle , high


 *-------------------------------------------------------------------
 * Description: computes AAF distrib over variant positions from vcf extracted 6 fields, stores histogram file 'distr',
   computes possible percentage of contaminated mixture (two modes) and its confidence

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
#define NPAR 4                // Number of arguments
#define NTHR 5                // Number of fields (thR, thA, mu, std and ploidy (1=haploid 2=diploid))
#define NRID 6                // Number of fields (pos,Depth,DP4) in extracted_vcf files
#define MIN_COUNT_AFTERF 200  // Minimum variant count after filtering
#define DEFAULT_DIST 11       // Default bin width for hist to compute skewness

// ******** declarations of  functions
void usage(int code); // usage
int GetAf (int []);// computes AF percentage for a variant position
float GetMu( int data[]);
float GetStd( int data[]);
void skewnesses (int data[], float sk[], float skc[], int st[], int stc[]);
void MixMode(int mix[], int st[], int data[], int dist);
void MixModeComplement(int mixc[], int stc[], int data[], int dist);
int GetMaxPeakInterval(int data[], int n1, int n2);//------max peak from the interval
void Likely (float Lik[], float sk[],float skc[], int mixc[],float avDP4);

////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////
int main    (int argc, char *argv[]) {

    // verbose mode
    int verbose = 0;

    // file names
    char *thr_file, *extract_vcf_file, *mix_file, *distrib_file;

    // file handles
    FILE *thrFile, *extract_vcfFile, *mixFile, *distribFile;

    // values in threshold file
    int muD, stdD;
    int thR, thA, ploid;

    // to read vcfa/vcfq
    static const int line_size = 8192; // maximum line size
    char line[line_size];

    // values in vcfq file
    int DP4[4];
    int pos, D;

    int i;

    // stats to compute
    float sumDP4 ;//across all pos for mu std
    int avDP4; // of actual sDP4 depth after RA & bad regions filtering

    float filtered_left;

    // settings for skewnesses and their likelihoods: hardcoded within functions

    int dist = DEFAULT_DIST; // bin width for hist to compute skewness

    //---------------results to compute:
    int count_beforeF; // variant count before filtering
    int count_afterF; // variant count after filtering

    //---------------------------arrays
    int histAF[100]; // to store mafs percentages
    int st[4], stc [4]; // starts of summing intervals for histo
    float sk[3], skc[3]; // three skewnesses for mode0,1,2
    int mix[3], mixc[3]; // mixtures three modes and their complements
    float Lik[3]; // likelihood or confidences three modes

    static struct option long_options[] = 
        { {"verbose", 0, 0, 'v'},
          {"help", 0, 0, 'h'},
          {0, 0, 0, 0}
        };

    char c;
    while ( (c = getopt_long(argc, argv, "r:a:vh?", long_options, 0)) != -1) {
        switch (c) {
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
        thr_file         = argv[optind+0];
        extract_vcf_file = argv[optind+1];
        mix_file         = argv[optind+2];
        distrib_file     = argv[optind+3];
    }

    // open files for read/write
    thrFile = fopen(thr_file,"r");
    if (thrFile == NULL) {
        fprintf(stderr, "cannot open threshold_file %s: %s\n", thr_file, strerror(errno));
        exit(EXIT_FAILURE);
    }
    extract_vcfFile = fopen(extract_vcf_file,"r");
    if (extract_vcfFile == NULL) {
        fprintf(stderr, "cannot open filtered_vcf_file %s: %s\n", extract_vcf_file, strerror(errno));
        exit(EXIT_FAILURE);
    }
    mixFile = fopen(mix_file,"w");
    if (mixFile == NULL) {
        fprintf(stderr, "cannot open mixture_file %s: %s\n", mix_file, strerror(errno));
        exit(EXIT_FAILURE);
    }
    distribFile = fopen(distrib_file,"w");
    if (distribFile == NULL) {
        fprintf(stderr, "cannot open distribution_file %s: %s\n", distrib_file, strerror(errno));
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "get_mixture\n");

    // initialise AAF histogram
    for (i=0;i<100;i++) {
        histAF[i] = 0;
    }

    // read thr file - 5 fields theshold ref, threshold alt, mean depath, stddev depth and ploidy
    while (fgets(line, line_size, thrFile)) {
        int k = sscanf(line, "%d %d %d %d %d", &thR, &thA, &muD, &stdD, &ploid);
        if (k != NTHR) {
            // number of fields read not correct
            fprintf (stderr, "corrupt threshold file %s\n", line);
            exit(EXIT_FAILURE);
        }
    }
    fprintf(stderr, "ploidy= %d\n", ploid);

    // initialise counts
    count_beforeF = 0;
    count_afterF = 0;

    // initialise sums
    sumDP4 = 0;

    // read vcfq file - 6 fields position, depth and 4xdepth (ref/forward, ref/reverse, alt/forward and alt/reverse) 
    while (fgets(line, line_size, extract_vcfFile)) {
        int k = sscanf(line, "%d,%d,%d,%d,%d,%d", &pos, &D, &DP4[0], &DP4[1], &DP4[2], &DP4[3]);
        if (k != NRID) {
            // number of fields read not correct
            fprintf(stderr, "skipping malformed VCF line %s\n", line);
            continue;
        }
        count_beforeF++;

        // filter for DP4 separately for ref and alternative alleles and abnormaly large Depth
        if( DP4[0] >= thR && DP4[1]>= thR && DP4[2]>= thA && DP4[3]>= thA && D < (muD+2*stdD) ) {
            int AF;
            float sDP4;

            count_afterF++;

            // compute AF and update AAF histo
            AF = GetAf(DP4);
            histAF[AF-1]++;

            // calc sDP4 and update sumDP4
            sDP4 = DP4[0] + DP4[1] + DP4[2] + DP4[3];
            sumDP4 = sumDP4 + sDP4;
        }
    }

    fprintf(stderr, "count of variants before filtering = %d\n", count_beforeF);
    fprintf(stderr, "count of variants after filtering = %d\n", count_afterF);

    // initialise interval sums
    for (i=0;i<4;i++) {
        st[i] = 0;
        stc[i] = 0;
    }
    // initialise skwnesses, mix, conf for modes 0,1,2
    for (i=0;i<3;i++) {
        sk[i] = 0.0;
        skc[i] = 0.0;
        mix[i] = 0;
        mixc[i] = 0;
        Lik[i] = 0.0;
    }

    // initialise variables
    avDP4 = 0; // average depth

    // do we have enough data ?
    if (count_afterF < MIN_COUNT_AFTERF) {
      fprintf(stderr, "Not enough variants after filtering, count_afterF < %d\n", MIN_COUNT_AFTERF);

    } else {
        // stats after filtering
        avDP4 = ceil(sumDP4 / count_afterF); // average depth after filtering
        filtered_left = (100.0 * count_afterF / count_beforeF); // percentage left after filtering
        fprintf(stderr, "percentage left after filtering = %.2f\n", filtered_left);
        fprintf(stderr, "coverage after filtering = %d\n",  avDP4);

        // calculate skewness
        skewnesses(histAF, sk, skc, st, stc);

        // mixtures per mode; confidences per mixture
        MixMode(mix, st, histAF, dist);
        MixModeComplement(mixc, stc, histAF, dist);
        Likely(Lik, sk, skc, mixc, avDP4);

        if (verbose) {
            fprintf(stderr, "Results:\n");
            for (i=0;i<3;i++) {
                fprintf(stderr, "sk = %.2f\n", sk[i]);
            }
            for (i=0;i<3;i++) {
                fprintf(stderr, "skc = %.2f\n", skc[i]);
            }
            for (i=0;i<3;i++) {
                fprintf(stderr, "Mode %d: likelihood of mix = %.2f\n", i, Lik[i]);
            }
        }
    }

    // write distribution file
    for (i=0;i<101;i++) {
        if (fprintf(distribFile,"%d %d\n",i, histAF[i]) <= 0) {
            fprintf(stderr, "error writing distribution_file: %s\n", strerror(errno));
            exit(EXIT_FAILURE);
        }
    }
    
    // write mixture file
    if (fprintf(mixFile,"mix low freq=%d\nconfidence low freq=%.4f\n", mix[0], Lik[0]) <= 0 ) {
        fprintf(stderr, "error writing mixture_file: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    if (fprintf(mixFile,"mix middle freq= %d\nconfidence middle freq=%.4f\n",mix[1], Lik[1]) <= 0 ) {
        fprintf(stderr, "error writing mixture_file: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    if (fprintf(mixFile,"mix high freq= %d\nconfidence high freq=%.4f\n", mix[2], Lik[2]) <= 0 ) {
        fprintf(stderr, "error writing mixture_file: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    if (fprintf(mixFile,"AvActDepth=%d\nmin_depthR=%d\nmin_depthA=%d\n", avDP4, thR, thA) <= 0 ) {
        fprintf(stderr, "error writing mixture_file: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    // close files
    fclose(thrFile);
    fclose(distribFile);
    fclose(extract_vcfFile);
    fclose(mixFile);

    fprintf(stderr, "done get mixture.\n");
    return 0;
}

////////////////////////////////////////////////////
// usage
////////////////////////////////////////////////////
void usage(int code)
{
    FILE *usagefp = stderr;

    fprintf(usagefp, "get_mixture\n\n");
    fprintf(usagefp,
            "Usage: get_mixture [options] threshold_file filtered_vcf_file mixture_file distribution_file\n"
            "\n" "  applies minimum depth thresholds for reference and alternative alleles, calculates a histogram of allele frequencies and estimates mixture likelihoods\n");
    fprintf(usagefp, "\n");
    fprintf(usagefp, "  options:\n");
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
    for (i=1;i<4;i++) {
        sum += data[i];
    }
    AF1 = (data[2] + data[3]) / sum;
    AF = (int) (100*AF1 + 0.5);

    return AF;
}

////////////////////////////////////////////////////
// GetMu.c calculates mu for all non-zero elements of 100-vector
////////////////////////////////////////////////////
float GetMu( int data[])
{
    float mu, sum;
    int i;
    int cnz = 0;

    sum = data[0];
    for (i=1;i<101;i++) {
        if (data[i] > 0) {
           cnz++;
           sum += data[i];
        }
    }
    mu = sum / cnz;

    return mu;
}

////////////////////////////////////////////////////
// standard deviation
////////////////////////////////////////////////////
float GetStd( int data[])
{
    float mu, std, sum, sum2;
    int i;
    int cnz = 0;

    sum = data[0];
    sum2 = data[0] * data[0];
    for (i=1;i<101;i++) {
        if (data[i] > 0) {
            cnz++;
            sum += data[i];
            sum2 += data[i] * data[i];
        }
    }
    mu = sum / cnz;
    std = sqrt( (sum2 / cnz) - (mu * mu) );
    return std;
}

////////////////////////////////////////////////////
// new skewnesses for three modes
////////////////////////////////////////////////////
void skewnesses (int data[], float sk[], float skc[], int st[], int stc[])
{
    int i, j, d, n1, n2;
    float s[4];
    float sc[4];

    d = DEFAULT_DIST;

    //1.  starts, starts complementary
    st[0] = 5;
    stc[0] = 51;
    for (i=1;i<4;i++) {
        st[i] = st[i-1]+d;
        stc[i] = stc[i-1]+d;
    }

    //2. sums between starts: 3 sums and 3 sums complementary
    for(i=0;i<4;i++) {
        n1 = st[i];
        n2 = stc[i];
        s[i] = 0;
        sc[i] = 0;
        for (j=n1;j<n1+d;j++) {
            s[i] = s[i] + data[j];
        }
        for (j=n2;j<n2+d;j++) {
            sc[i] = sc[i] + data[j];
        }
    }

    //3.     // compute skewnesses
    for (i=0;i<3;i++) {
        if (s[i+1]>0) {
            sk[i]=s[i]/s[i+1];
        }
        if (sc[i]>0) {
            skc[i]=sc[i+1]/sc[i];
        }
    }
}

////////////////////////////////////////////////////
// mixture mode
////////////////////////////////////////////////////
void MixMode(int mix[], int st[], int data[], int dist)
{
    int i;
    int n1, n2;
    int m;

    for (i=0;i<3;i++) {
        n1 = st[i];
        n2 = n1 + dist;
        m = GetMaxPeakInterval(data, n1, n2);
        mix[i] = m;
    }
}

////////////////////////////////////////////////////
// mixture mode complement
////////////////////////////////////////////////////
void MixModeComplement(int mixc[], int stc[], int data[], int dist)
{
    int i;
    int n1, n2;
    int m;

    for (i=0;i<3;i++) {
        n1 = stc[i+1];
        n2 = n1 + dist;
        m = GetMaxPeakInterval(data, n1, n2);
        mixc[i] = m;
    }
}

////////////////////////////////////////////////////
// max peak from the interval
////////////////////////////////////////////////////
int GetMaxPeakInterval(int data[], int n1, int n2)
{
    // input data = AAF
    int max_peak;
    int i, mixx, k, j;
    int peak[100], perc_peak[100];

    // initialise mix, in case no peaks exist
    mixx = 0.0;

    // find max peak among peaks in the interval only, data is maf(AAF) vector

    // what if no peaks?
    peak[0] = 0;
    perc_peak[0] = 0;

    // find peaks
    k = 0;
    for (i=n1;i<n2;i++) {
        if ( ((data[i]-data[i-1]) >=0 ) && ((data[i+1]-data[i]) < 0)) {
            peak[k] = data[i];
            perc_peak[k] = i;
            k++;
        }
    }

    if (k == 1) {
        // one peak
        max_peak = peak[0];
        mixx = perc_peak[0];
    } else if (k > 1) {
        // more than one peak
        max_peak = peak[0];
        for (j=1;j<k;j++) {
            if (max_peak < peak[j]) {
                max_peak = peak[j];
            }
        }

        // find percentage giving max_peak, first one if several are max
        for (j=0;j<k;j++) {
            if (peak[j] == max_peak) {
                mixx = perc_peak[j];
                break;
            }
        }
    }

    return mixx;
}

////////////////////////////////////////////////////
// confidences modes 0,1,2
////////////////////////////////////////////////////
void Likely (float Lik[], float sk[],float skc[], int mixc[],float avDP4)
{
    int i;
    float l, lc, lp;

    for (i=0;i<3;i++) {
        lp = 0.0;
        if (mixc[i]>0) {
            lp = skc[2-i];
        }
        l = sk[i];
        lc = skc[2-i];
        Lik[i] = (1-1/avDP4) * (l+lc+lp) / 3;
    }
}






