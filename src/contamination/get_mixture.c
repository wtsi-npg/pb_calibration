 /*  File: get_mixture.c // computes maf distribution over variant positions from vcfq files
 * Authors: designed by Irina Abnizova (ia1)edited by Steve Leonard (srl)
 *

  Last edited:

  28 April-to remove training printfs, to output mixture=0 if not any
   25 April: incorporate info about (100-mix)% to increase proly of detection
  24 April-output in infoAF 4 fields for each variant position, last is its actual coverage sDP4
  22 April: max sginif or first max if not signif  + all 100 AAF to see 'complimentary'
  10 April:   if ( DP4[0] >= THR[0] && DP4[1]>= THR[0] && DP4[2]>= THR[1] && DP4[3]>= THR[1])// 0 ref alt are ok
              thr_sk
  10 March 2014 for skewness to first 25% AAF
  3 feb-reads thr from a file; 29 Jan 2014  - get autom thr, starting with default ones 1,5
  still old fields of vcf input
  21 Nov 2013, for two thresholds

 *-------------------------------------------------------------------
 * Description: computes maf distrib over variant positions from just extracted 6 fileld

 * Exported functions:
 * HISTORY:

skewness to first 25% of AAF  ----------introduced 10 March
computes AF (alternative to Ref allele frequency) for each error base call from vcf

applies minimum depth threshold, both for Ref and Alternative alleles fwd rev (total min depth will be
4*MIN_DEPTH

input1 thrfile  (thRR thAA)
input2: extracted vcf,


output1: info AF.txt=

pos  %AF Q30frac sumDP4 (for this position)
1546 100 0.89  101
4469 100 0.86  90
4532 100 0.78  87
4559  97 0.95  60

output2: distrAF.tsigt=histogram of AF with significance for each bin
% count signif
8 0 0
9 0 0
10 385 1
11 780 3


usage:./get_mixture ra_i_j.thr extracted.vcf .infAF .distrAF

*/

#include <stdio.h>
//#include "conio.h"
#include <time.h>
#include <math.h>
#include <string.h>

// ******** predefined user constants, 22 october

#define NRID 6      // Number of variant ids=fields: pos, chromosome, Depth,DP4,(no PV4+maf) in esigtracted_vcf files
#define NPAR 5      // masig Number of PARameters and arguments
#define Nthr 2      //  Number of thresholds


// ******** declarations of  functions
int   GetAf (int []);// computes AF percent for a variant position
float   GetMu( int data[]);
float   GetStd( int data[]);

int GetMaxPeak2(int data[]);
float ComplimentaryPeaks(int data[], int perc);

//=====================================================MAIN

int main    (int argc, char *argcv[]) {

    // flags
    int firstSixAreOK  = 1;
    int canWriteAF        = 1;
    int canWriteDistrib      = 1;

    FILE *extract_vcfFile, *afFile, *distribFile, *thrFile; //file handles

    int n,i,perc;//, cntZero=0, counter=0;
    //char chrom[50];
    int count_afterF=0;// after filtering with min depths
    int count_beforeF=0;// before filtering with min depths
    int DP4[4];//
    int D,pos,AF;// Depth, genome pos of error
    int sumD=0,sumDP4=0;//across all pos!
    int sumDbefore=0;
    int thR,thA;// from file
    int less50;
    int more25;// sum AF  (25,50)
    int less25;
    int bw25_40;
    int sDP4;//act coverage of a base

    float sDP4f;//act coverage of a base as float
    float mu,std;//mean and std of alternative frequency, AF
    float avDbefore;
    float avD,avDP4,Q30frac;

    float perc_left;
    float sig,conf;//confidence (0,1)
    float sk25;//skewness to first 25% of AAF  ----------introduced 10 March
    float thr_sk;
    float fprolly;

    int histAF[100];// to store mafs percentages
    int THR[2];// define from input1 thr file


    if(argc < NPAR)//four input_output files are submitted
    {
        printf("not enough of parms |input_output files\n");
        printf("usage:./get_signifAF name.thr input2 output1 output2\n");
        return -1;
    }


    thrFile=fopen(argcv[1],"r");
		    if (thrFile == NULL) {
		      printf("cannot open first input _threshold file %s\n", argcv[1]);
		      return -1;
    }

    extract_vcfFile=fopen(argcv[2],"r");
    if (extract_vcfFile == NULL) {
      printf("cannot open first input _vcf.tsigt file %s\n", argcv[2]);
      return -1;
    }

    afFile=fopen(argcv[3],"w");
    if (afFile == NULL) {
    printf("cannot open first output1 infoAF.tsigt file %s for writing\n", argcv[3]);
    return -1;
    }

    distribFile=fopen(argcv[4],"w");
    if (distribFile == NULL) {
    printf("cannot open second output file distAF %s for writing\n", argcv[4]);
    return -1;
    }

    printf("get_peak_thrFile\n");

    // give value for thr_sk=    0.55 for human; 1.0 for pathogens
             thr_sk=0.6;

    // initiate zero vector for histograme
			 for(i=0;i<100;i++)//
	         {
	              histAF[i]=0;
		     }
    // initiate zero vector for threshold
			 for(i=0;i<1;i++)//
			 {
			 	  THR[i]=0;
		     }

		// scan thr file	and assign current precomputed thresholds for filtering
    while( (n = fscanf(thrFile,"%d %d", &thR, &thA )) >= 0)
         // until the end of the input threshold  file
    {
		      if (n != Nthr) // incorrect format, NCOLS=number of columns in pipeCT input
		      {
		      printf ("corrupted input thrFile format\n");
		      return -1;
		      }

		THR[0]=thR;
		THR[1]=thA;
	}

	        printf("recommended precomputed Ref threshold= %d\n", thR);
	        printf("recommended precomputed Alternative threshold= %d\n", thA);

	//6 fields of input file	// no PV4

    while( (n = fscanf(extract_vcfFile,"%d,%d,%d,%d,%d,%d", &pos, &D, &DP4[0], &DP4[1], &DP4[2], &DP4[3])) >= 0 && firstSixAreOK == 1 && canWriteAF == 1)
    // read the Read Title
    {
        if( n != NRID )     // incorrect format
        {
            firstSixAreOK = 0;
            break;
        }
           count_beforeF++;
           sumDbefore=sumDbefore+D;

           // f1 Josie : filter for DP4 separately for ref and alternative alleles

            if ( DP4[0] >= THR[0] && DP4[1]>= THR[0] && DP4[2]>= THR[1] && DP4[3]>= THR[1])// 0 ref alt are ok
           {
                AF = GetAf(DP4);
                histAF[AF-1]++;
                /*for(i=0;i<100;i++)//
			    {
				   if (AF==i+1)
				   histAF[i]=histAF[i]+1;//  sh be +1
		        }
		        */
		        // compute average run depth default (Q13), D across pos,sumD
				// compute average run depth actual after Q30-25, sum(DP4) across pos,sumDP4
				        count_afterF++;
				        sumD=sumD+D;
				        sDP4=DP4[0]+DP4[1]+DP4[2]+DP4[3];
				        sumDP4=sumDP4+sDP4;

				        sDP4f=sDP4;
				        Q30frac=sDP4f/D;

                // output check:4 columns position AF% highQ frac, actual depth=sDP4
                if( fprintf(afFile,"%d %d %.2f %d\n", pos, AF, Q30frac,sDP4) < 0 )
                 {
                 canWriteAF = 0;
                 break;
                 }

             }// end min_depth filter of DP4

        if (canWriteAF == 0) {break;}

        // removing space symbols before carriage return
        do  {
            n = fgetc (extract_vcfFile);
        }
        while ((char)n != '\n');

    } //END of outer while loop across vcfq file

    // to find significat peak of AF histo!
        mu=GetMu(histAF);
        std=GetStd(histAF);


      less50=histAF[0];
      for(i=0;i<101;i++)//we used to consider only 49, because we need only lower freq, esp for diploid organisms
      {
		  sig=0;// significance of histo bin
		  if (std > 0) {
		  sig=(histAF[i]-mu)/std;
          less50+=histAF[i];
	      }

	     if( fprintf(distribFile,"%d %d %.2f\n",i, histAF[i],sig) <= 0 ) {
         canWriteDistrib = 0;
	     }

       }//  for i=100 cycle

       less25=histAF[0];
	   for(i=0;i<25;i++)//only  because we need only lower freq, esp for diploid organisms
	   {
	   		less25+=histAF[i];
	   }
	   bw25_40=histAF[26];
	   for(i=27;i<40;i++)//only  because we need only lower freq, esp for diploid organisms
	   {
	   	   		bw25_40+=histAF[i];
	   }

//compute 25% AF skewness
       //if (bw25_40==0){
	   sk25=0.00;// if no bases for bw25_40
       //}
	   if (bw25_40 >0 )
	   {
	   sk25=1.0*less25/bw25_40;// if proportion of less 25% AF is high (more 1?), more likely a mixture
       }

      // find  first max peak among first 25% AAF

          perc=GetMaxPeak2(histAF);//peaks in first 25% AAF

       // condition on mixture percent:
       //if more than 25, it is likely to be diploid allele  :  17 Nov+ amount AF<25%  skewed (too high)

       if (perc < 25 && sk25 >= thr_sk)// thr_sk was 0.95 for path, 0.62 for human
       {
		   printf("percent mixture= %d\n", perc);
	   }
		if (perc >=25 || sk25 < thr_sk )
		{
		   perc=0.0;
		   printf("percent mixture= %d\n", perc);
		   printf("not likely there is a low percent mixture here at this min depth\n");
        }

       // compute avarege depth :default and after filtering
      avD=sumD/count_afterF;
      avDP4=sumDP4/count_afterF;//actual after RA filter
      avDbefore=sumDbefore/count_beforeF;
      perc_left=(100.0*count_afterF/count_beforeF);
      printf("actual average depth= %.2f\n", avDP4);

       // confidence of mixture estimating
       conf=(1-1/avDP4)*0.7;

       if (sk25>0.8)// larger any how!
       {conf=(1-1/avDP4)*0.9;//(sigsig[perc]+1)*(1-1/avDP4);// 0.8 is MY constant
       }

       fprolly=ComplimentaryPeaks(histAF, perc);

         conf=conf*fprolly;
         if (conf>1)
         {conf=1.0;}

         printf("Confidence of non-zero mixture = %.2f\n", conf);
	     printf("skewness of first 25 perc AF = %.2f\n", sk25);

    fclose(distribFile);
    fclose(extract_vcfFile);
    fclose(afFile);

    if( firstSixAreOK  == 0 ||
        canWriteAF == 0 || canWriteDistrib ==0) {
        printf ("Error during execution. Details: \n");
        printf ("\tfirstSixAreOK %d\n",  firstSixAreOK);
        printf ("\tcanWriteAF %d\n",        canWriteAF);
        printf ("\tcanWriteDistrib %d\n",      canWriteDistrib);
        printf ("Execution aborted\n");
        return -1;
    }

    printf(" done.\n");
    return 0;
}//main

/***************************************************** Functions******************************/


////////////////////////////////////////////////////
// Calculates AF=alternative allele frequebcy from data=DP4 counts
// input - array (4-vector) of DP4 counts
//output -one float number from (0,1)
////////////////////////////////////////////////////
int GetAf (int data[])
{
    float Isum,AF1;
    int i,AF;

    Isum = data[0];
    for (i=1; i<4; i++)
    {
        Isum += data[i];
    }
    AF1=(data[2]+data[3])/Isum;
    AF=(long) (100*AF1+0.5);
    //(long) (sig+0.5)
    return AF;
}
// GetMu.c calculates mu for all non-zero elements of 100-vector
float GetMu( int data[])
 {
	 float mu,Isum;
	 int i;
     int cnz=0;

     Isum = data[0];
     //Isum2 = data[0]*data[0];
	 for (i=1; i<101; i++)
	 {
		 if (data[i] > 0)
		 {
		 cnz++;
	     Isum += data[i];
	     //Isum2 += data[i]*data[i];
	     }
      }
      mu=Isum/cnz;
      //std=Isum2/cnz-mu*mu;
      return mu;
}

// ========================stand deviation
float GetStd( int data[])
 {
	 float mu,std,Isum,Isum2;
	 int i;
     int cnz=0;

     Isum = data[0];
     Isum2 = data[0]*data[0];
	 for (i=1; i<101; i++)
	 {
		 if (data[i] > 0)
		 {
		 cnz++;
	     Isum += data[i];
	     Isum2 += data[i]*data[i];
	     }
      }
      mu=Isum/cnz;
      std=sqrt(Isum2/cnz-mu*mu);
      return std;
}


//----------------------max from the first 25
 int GetMaxPeak2(int data[])
 {

	 int Imax;
     int i,perc,k;
     int peak[25], per_peak[25];

     // initialise perc

     perc=0.0;

// find max peak among peaks, first 25 percent only

     // find peaks
     k=0;
     for (i=1; i<25; i++)
     {
		 if ( ((data[i]-data[i-1]) >0 ) && ((data[i+1]-data[i]) <=0))// && ( Imax < data[i]))//here peak in [i]
         {

             peak[k] = data[i];
             per_peak[k]=i;
             //printf("peak = %d\n", peak[k]);
             //printf("per_peak = %d\n", per_peak[k]);
             k=k+1;

         }
	 }
// find max peak
	 Imax = peak[0];

	 	      for (i=1; i<k; i++)
	 	      {
	 	          if( Imax < peak[i])
	 	              Imax = peak[i];
	 	              //perc=i;
	           }

       //printf("max peak = %d\n", Imax);
   // find percent giving Imax, first one if several are max
       for (i=0; i<k; i++)
  	   {
  	          if( peak[i]==Imax)
  	            { perc=per_peak[i];

  	              break;}
       }

       return perc;

}

//=================================================fihg AAF complimentary peaks, (100-perc)AAF if any
float ComplimentaryPeaks(int data[], int perc)
 {

	 float fprolly;
     int i,k;
     int peak[25], per_peak[25];
     int around;

     around=10;// frequency may vary up to this amount of %
     fprolly=1.0;

// find max peak among peaks, last 75% (because as mixture was first 25 percent only)

     // find peaks
     k=0;
     for (i=75; i<101; i++)
     {
		 if ( ((data[i]-data[i-1]) >0 ) && ((data[i+1]-data[i]) <=0))// && ( Imax < data[i]))//here peak in [i]
         {

             peak[k] = data[i];
             per_peak[k]=i;
             //printf("peak = %d\n", peak[k]);
             //printf("per_peak = %d\n", per_peak[k]);
             k=k+1;

         }
	 }

	 // see if these peaks are complimentary to detected perc=% mixture <25%
	 for (i=0; i<k; i++)
	 {
	     if ( (per_peak[i] > (100-perc-around)) && (per_peak[i] < (100-perc+around)))
	     fprolly=1.1;
	 }
 //&& dp(i,1)>(100-perc-around) & dp(i,1)<(100-perc_ms+m1)


       return fprolly;

}