/*  File: get_thr_ploidy.c // filters 1 1 vcfq file and computes av std Depth  depending on Depth of vcfq file

 * Authors: designed by Irina Abnizova (ia1) and edited by Steve Leonard(srl)
 *
  Last edited:

  27 April output of the Actual Depth distribution

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
#define NPAR 4                // Number of arguments
#define THR 1                 // Default reference allele threshold
#define THA 1                 // Default alternative allele threshold
#define NRID 6                // Number of fields (pos,Depth,DP4) in extracted_vcf files
#define MIN_COUNT_AFTERF 200  // Minimum variant count after filtering

// ******** declarations of  functions
void usage(int code); // usage
int GetAf (int []); // computes AF percent for a variant position
int GetActualDepth (int []);// computes actual depth for a variant position, sumDP4
int GetMaxIndex (int [], int i1, int i2); // percent giving max AAF
int ploidy (int ind); // defines ploidy from results of GetMaxIndex ploid=2 if ~50%

int GetMax (int hist[], int nbins);
int GetMin(int hist[], int nbins);
int EstimateStd (int hist[],int bins[],int nbins, int mu, int height);// returns sd
void GetThresholds(int thr[],int depth_min, int depth_max, int avDP4, int stdDP4, int mode_depth, int sdmo, int bin_width);


////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////
int main (int argc, char *argv[]) {

    // depth_mins for Ref and Alternative alleles
    int thR = THR;
    int thA = THA;

    // verbose mode
    int verbose = 0;

    // file names
    char *extract_vcf_file, *thr_file, *filtered_vcf_file, *depth_distr_file;

    // file handles
    FILE *extract_vcfFile, *thrFile, *filtered_vcfFile, *depth_distrFile;

    // to read vcfq
    static const int line_size = 8192; // maximum line size
    char line[line_size];

    // values in vcfq file
    int DP4[4]; // 4xdepths
    int pos, D; // genome pos of error, depth
    int bin_width;// for depth histo

    int count_beforeF; // variant count before filtering
    int count_afterF; // variant count after filtering

    // arrays
    int histAF[101]; // maf percentages
    int histAAD[101];// distribution of average actual depth DP4 per one percentage of AAF
    int bins_depth[101];
    int hist_depth[101];


    // stats
    float sumDP4, sumDP4_2; // sum of depths and sum of square depths used for mu/std calculation
    float sDP4;


    // to compute for outputs
    int ploid; // ploidy
    int avDP4, stdDP4; // mean and std of depth after filtering
    //int thRR, thAA; // estimated thresholds based on actual depths avDP4 and stdDP4
    int depth_max,depth_min;
    int thr_low,thr_up;


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
        depth_distr_file  = argv[optind+3];
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

    //*depth_distr_file

    depth_distrFile = fopen(depth_distr_file,"w");
        if (depth_distrFile == NULL) {
            fprintf(stderr, "cannot open depth_distr_file %s: %s\n", depth_distr_file, strerror(errno));
            exit(EXIT_FAILURE);
    }

    fprintf(stderr, "get_thr_ploidy\n");

    // initialise depth histograms
        bin_width=10;
        for (i=0;i<101;i++) {
            hist_depth[i] = 0.0;
            bins_depth[i]=1+i*bin_width;
        }

    // initialise AAF and AAD histograma
    for (i=0;i<101;i++) {
        histAF[i] = 0;
        histAAD[i]=0;
    }

    // initialise counts
    count_beforeF = 0;
    count_afterF = 0;

    // initialise sums
    sumDP4 = 0;
    sumDP4_2 = 0;
    sDP4=30.0;

    depth_max=sDP4;
    depth_min=sDP4;


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
            int AF,AAD;//altrenative allele frequency, actual average depth per variant
            //float sDP4;


            count_afterF++;

            // compute AF , AAD and update histos
            AF = GetAf(DP4);
            histAF[AF-1]++;
            AAD =GetActualDepth(DP4);
            histAAD[AF-1]=histAAD[AF-1]+AAD;


            // calc sDP4 and update sumDP4, sumDP4_2
            sDP4 = DP4[0] + DP4[1] + DP4[2] + DP4[3];
            sumDP4 +=  sDP4;
            sumDP4_2 +=  sDP4*sDP4;

            if (depth_max < sDP4)
            depth_max=sDP4;

            if (depth_min > sDP4)
            depth_min=sDP4;

            for (i=0;i<100;i++) {
                if ((sDP4 >= bins_depth[i]) && (sDP4 < bins_depth[i+1]))
                    hist_depth[i]++;
            }

            // write to filtered vcf file
            if (fprintf(filtered_vcfFile,"%d,%d,%d,%d,%d,%d\n",pos,D, DP4[0], DP4[1],DP4[2], DP4[3]) <= 0) {
                fprintf(stderr, "error writing filteredVCF file: %s\n", strerror(errno));
                exit(EXIT_FAILURE);
            }
        }
    }
         fprintf(stderr, "min depth after filtering 1 1= %d\n", depth_min);
          fprintf(stderr, "max depth after filtering 1 1= %d\n", depth_max);

        for(i=0;i<101;i++){
            if (histAF[i] > 0) {
                histAAD[i]=ceil(histAAD[i]/histAF[i]);

            }
        }
        // write distribution file
        for (i=0;i<100;i++) {
            if (fprintf(depth_distrFile,"%d %d\n",bins_depth[i], hist_depth[i]) <= 0) {
                    fprintf(stderr, "error writing depth distribution_file: %s\n", strerror(errno));
                    exit(EXIT_FAILURE);
           }
        }


    fprintf(stderr, "count of variants before filtering = %d\n", count_beforeF);
    fprintf(stderr, "count of variants after filtering = %d\n", count_afterF);

    // initialise values
    avDP4 = 0; // mean depth
    stdDP4 = 0; // stddev depth
    ploid = 1; // ploidy
    //avDP41 = 0; // adjusted mean depth
    //thRR =100; // estimated reference threshold: set up stupidly high in case not enough data
    //thAA =100; // estimated alternative threshold
    thr_low=100;
    thr_up=0;

    // do we have enough data ?
    if (count_afterF < MIN_COUNT_AFTERF) {
      fprintf(stderr, "Not enough variants after filtering, count_afterF < %d\n", MIN_COUNT_AFTERF);

    } else {
        int ind, ind50;
        int mode_depth,ind_mode;
        int height,sdmo;
        float perc_left;
        //float ratio;

        int thr[2];//for output low up thr

        //initialise thr up low
        for (i=0; i<2;i++){
            thr[i]=0;
        }


        // ploidy detection
        ind = GetMaxIndex(histAF, 1, 80);// percent giving max peak in hist in interval (1,80)%
        fprintf(stderr, "percent of max AAF= %d\n", ind);


     // -----------compute params of main mode: mode and its sdmo
        // mode is the bin with its index
            height = GetMax(hist_depth,100);
            ind_mode = GetMaxIndex(hist_depth, 1, 100);
            mode_depth=bins_depth[ind_mode];

            //nbins=100;
            sdmo = EstimateStd(hist_depth,bins_depth,100,mode_depth,height);//for depth hist

         fprintf(stderr, "mode depth= %d\n", mode_depth);//
         fprintf(stderr, "std of mode depth= %d\n", sdmo);//

        ind50 = GetMaxIndex(histAF, 40, 60);//percent giving max peak in hist in interval (40,60)%
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


       //------------------------------ automated working out of thresholds: thr_low thr_upper

        GetThresholds(thr,depth_min,depth_max, avDP4, stdDP4, mode_depth, sdmo, bin_width);
        thr_low=thr[0];
        thr_up=thr[1];
   }// if enough data

    // write threshold file
    if (fprintf(thrFile, "%d %d %d %d %d\n", thr_low, thr_up, avDP4, stdDP4, ploid) <= 0) {
        fprintf(stderr, "error writing threshold_file: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    // close files
    fclose(extract_vcfFile);
    fclose(thrFile);
    fclose(filtered_vcfFile);
    fclose(depth_distrFile);

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
// Calculates actual depth  from data=DP4 counts
// input - array (4-vector) of DP4 counts
// output -sumDP4 into array of percentages 0,..,100
////////////////////////////////////////////////////
int GetActualDepth (int data[])
{
    float sum;
    int i;

    sum = data[0];
    for (i=1; i<4; i++) {
        sum += data[i];
    }

    return sum;
}

////////////////////////////////////////////////////
// find position of the max value in a sub-interval (i1<=i<=i2) of an array
////////////////////////////////////////////////////
int GetMaxIndex (int data[], int i1, int i2)
{
  int data_max;
  int ind, i;

  data_max = data[i1];
  ind = i1;
  for (i=i1; i<i2; i++) {
      if( data_max < data[i]) {
          data_max = data[i];
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

////////////////////////////////////////////////////
// estimate std of a mode
// input - histogram, bins, mode,its height
// output - std of this mode
////////////////////////////////////////////////////
int EstimateStd (int hist[],int bins[],int nbins, int mu, int height)
{
    int n;
    int ma_bin = -1;
    int mi_bin = -1;
    float threshold = 0.5*height;

    int sd;

    for (n=0; n<nbins; n++)
    {
        if (bins[n] >= mu)// to the Right of mode
        {
            if ((hist[n] > threshold) & (hist[n] >= hist[n+1]))// added monotonity condition 29 Jan
            {
                ma_bin = bins[n];
            }
            else
            {
                break;
            }
        }
    }

    for (n=nbins-1; n>=0; n--)
    {
        if (bins[n] <= mu)// to the Left of mode
        {
            if ((hist[n] > threshold) & (hist[n-1] <= hist[n]))// added monotonity condition 29 Jan
            {
                mi_bin = bins[n];
            }
            else
            {
                break;
            }
        }
    }

    sd = ceil(0.5 * (ma_bin - mi_bin));

    if (mu==bins[0])
    sd = (ma_bin - mi_bin);
    if (mu==bins[nbins-1])
    sd = (ma_bin - mi_bin);


    return sd;
}
///////////////////////
//// ---------  maximum value in a vector
/////
int GetMax (int hist[], int nbins)
{
    float Hmax;
    int i;

    Hmax = hist[0];
    for (i=0; i<nbins; i++)
    {
        if( Hmax < hist[i])
            Hmax = hist[i];

    }

    return Hmax;
}
///////////////////////
//// ---------  minimum value in a vector
/////
int GetMin(int hist[], int nbins)
{
    float Hmin;
    int i;

    Hmin = hist[0];
    for (i=0; i<nbins; i++)
    {
        if( Hmin > hist[i])
            Hmin = hist[i];

    }

    return Hmin;
}

///////////////////////
//// ---------  work out threshold up and low from depth_hist and its params: updates thr array
// input int thr[]  = thr_up thr_low
/////
void GetThresholds(int thr[],int depth_min, int depth_max, int avDP4, int stdDP4, int mode_depth, int sdmo, int bin_width)
{
      int i;
      float cvv;
      float thrL,thrU;
      int data[2];

      //initialise thresholds
      thrL=depth_min*0.25;
      thrU=depth_max;

      // depending on depth distribution: its symmetry, range, variation
      cvv=0.0;
      if (avDP4 > 0){
          cvv=1.0*stdDP4/avDP4;
      }
        fprintf(stderr, "CV of depth= %.2f\n", cvv);

        if (cvv>1){
         fprintf(stderr, "some regions are stupidly deep covered\n");
        }


      //1 assymetric distribution

      //1.1 skewed to the left: (avDP4-mode_depth) > bin_width

      if ((avDP4-mode_depth) > bin_width){
          data[0]=depth_min;
          data[1]=mode_depth-sdmo;
          thrU=GetMax(data,2);

          data[0]=depth_max;
          data[1]=avDP4+3*stdDP4;
          thrL=GetMin(data,2);
      }

      //1.2 skewed to the right: (avDP4-mode_depth) < - bin_width

       if ((avDP4-mode_depth) < - bin_width){
                data[0]=depth_min;
                data[1]=avDP4 - 2*stdDP4;
                thrU=GetMax(data,2);

                data[0]=depth_max;
                data[1]=mode_depth+sdmo;
                thrL=GetMin(data,2);
        }

       //1.3 symmetrical: mu+-3 std


         if ( (avDP4-mode_depth) >= - bin_width &&  (avDP4-mode_depth) <= bin_width  ){
                 data[0]=depth_min;
                 data[1]=avDP4 - 2*stdDP4;
                 thrU=GetMax(data,2);

                 data[0]=depth_max;
                 data[1]=avDP4 + 2*stdDP4;;
                 thrL=GetMin(data,2);
           }



        thr[0]=ceil(thrU*0.25);
        thr[1]=ceil(thrL);

     for (i=0;i<2;i++){
     fprintf(stderr, "thr filt low_up = %d\n", thr[i]);

    }

}




















