 /*  File: get_mixture.c // computes and analyses AAF( maf) distribution over variant positions from vcfq files
 * Authors: designed by Irina Abnizova (ia1)edited by Steve Leonard (srl)
 *

  Last edited:
  19 Jan confidences
  mix1   conf1  low freq allele mixture <25%
  mix2   conf2  high freq allele mixture (25% , 50%)


  11 Jan 2016 : look at ploidy, get hard codes inside corresponding functions
  void get_bins_ploidy	(int bins_ploid[], int plo);// // how to consider first n25 depending on ploidy
  void get_skew_ploidy	(float thr_sk[], int plo);// // how to consider skew25, skew75 depending on ploidy

  5 January 2016: picks up automaticaly if it is diploid or haploid
  17 December 2015: two modes 1. less 25%  2. b/w 30 and 40 percentage
  13 June -output sk25,sk75
  12 june-- final mixture file; modified as output name.mix:  mix_percentage=17
possible_mix_percentage= 17
confidence_nonzeroMix=0.64
AvActDepth=52 (avDP4)
min_depthR=2
min_depthA=2

  30 May last 75 bins implemented


 *-------------------------------------------------------------------
 * Description: computes AAF distrib over variant positions from vcf extracted 6 fields, stores histogram file 'distr',
   computes possible percentage of contaminated mixture (two modes) and its confidence

 * Exported functions:
 * HISTORY:

skewness to first 25% of AAF  ----------introduced 10 March
computes AF (alternative to Ref allele frequency) for each error base call from extracted vcf

applies minimum depth threshold, both for Ref and Alternative alleles fwd rev (total min depth will be
4*MIN_DEPTH

input1 thrfile  (thRR thAA  mu  std  ploidy )
                 3     2    82  110     2

input2: extracted vcf: pos, D , DP4 (DP4 is four fields)


output1: mix_conf AF.txt=
mixture confidence
12        0.87

output2: .distrAF=histogram of AF with significance for each bin
% count
8 0
9 0
10 385
11 780


usage:./get_mixture rams.thr extracted.vcf name.mix name.distrAF

*/

#include <stdio.h>
//#include "conio.h"
#include <time.h>
#include <math.h>
#include <string.h>

// ******** predefined user constants, 22 october

#define NRID 6      // Number of variant ids=fields: pos, chromosome, Depth,DP4,(no PV4+maf) in esigtracted_vcf files
#define NPAR 5      // masig Number of PARameters and arguments
#define Nthr 5     //  Number of thresholds and muD1 stdD1 12 June + ploidy=1 (haploid) (or 2-diploid)


// ******** declarations of  functions
int     GetAf (int []);// computes AF percentage for a variant position
void get_bins_ploidy	(int bins_ploid[], int plo);// // how to consider first n25 depending on ploidy
void get_skew_ploidy	(float thr_sk[], int plo);// // how to consider skew25, skew75 depending on ploidy
float   GetMu( int data[]);
float   GetStd( int data[]);
float   GetSkew25( int data[],int n25, int n49);
float   GetSkew75( int data[],int n51, int n75);

float Confidence25(float sk25,float sk75,int mix25,int mixc25,float thr25,float thr75, float avDP4);
float Confidence25_50(float sk25,float sk75,int mix25_50,int mixc25_50,float thr25,float thr75, float avDP4);


int GetMaxPeak25(int data[], int n25);
int GetPeakComplement25(int data[], int n75);

int GetMaxPeak25_50(int data[], int n25, int n49);
int GetPeakComplement25_50(int data[], int n51, int n75);


//=====================================================MAIN

int main    (int argc, char *argcv[]) {

    // flags
    int firstSixAreOK  = 1;
    int canWriteMix        = 1;
    int canWriteDistrib      = 1;

    FILE *thrFile, *extract_vcfFile, *mixFile, *distribFile ; //file handles

        int n,i;
        int dist;// distance bewteen mix and (100-mix_compl): should be 0 ideally

      // ---------- params from input1 rams(now four+1 values thrR thrA muD stdD ploidy)
		 int muD,stdD;
		 int thR,thA, ploid,plo;

		 int THR[2];// compute from input1 thr file
	     int msD[2];
	     int bins_ploid[4];// 4-vector of bins depending on ploidy
		 float thr_sk[2];//2-vector of sk25,sk75 depending on ploidy


    //---------------params from input2 vcfq
	    int DP4[4];//
	    int D,pos;// Depth, genome pos of error

	    //----------stats to compute
	    int avDP4;// of actual sDP4 depth after RA & bad regions filtering
	    int sumDP4=0;//across all pos for mu std
	    int sumDbefore=0;
	    float avDbefore;

		float sDP4;
		float mix25_left;

       // ----------------settings for skewness and their likelihoods: hardcoded within functions
       // mode 1: mix<=25%
       int n25,n49;
       int n51,n75;
       // mode 2 : 30< mix < 45
       int n2_25,n2_49;
       int n2_51,n2_75;
       // skewness thresholds
       float thr25,thr75;


      //---------------results to compute:
       int count_afterF=0;// after filtering with min depths
	   int count_beforeF=0;// before filtering with min depths

       int mix25,mixc25;
       int mix25_update;
       int mix25_50,mixc25_50;
       int AF;
       float conf,conf25,conf25_50;//confidence (0,1)
       float sk25,sk75;//skewness to first 25% of AAF  ----------introduced 10 March
       float likely75,likely25;

       //arrays
       int histAF[100];// to store mafs percentages

//-----------------------------------------open files for read/write
    if(argc < NPAR)//four input_output files are submitted
    {
        printf("not enough of parms |input_output files\n");
        printf("usage:./get_mixture input1.thr input2.vcfq output1.mix output2.distr\n");
        return -1;
    }

    thrFile=fopen(argcv[1],"r");
		    if (thrFile == NULL) {
		      printf("cannot open first input _threshold file %s\n", argcv[1]);
		      return -1;
    }

    extract_vcfFile=fopen(argcv[2],"r");
    if (extract_vcfFile == NULL) {
      printf("cannot open first input _.vcfq file %s\n", argcv[2]);
      return -1;
    }

    mixFile=fopen(argcv[3],"w");
    if (mixFile == NULL) {
    printf("cannot open first output1 mix_conf file %s for writing\n", argcv[3]);
    return -1;
    }

    distribFile=fopen(argcv[4],"w");
    if (distribFile == NULL) {
    printf("cannot open second output file distAF %s for writing\n", argcv[4]);
    return -1;
    }

    printf("get_mixture: parameters\n");

    // ===SETTINGS     parameters for skewness measurements and peak consideretion: if Diploid then more noise around 50%
           dist=10;
    // initiate zero vector for histograme
			 for(i=0;i<100;i++)//
	         {
	              histAF[i]=0;
		     }
    // initiate zero vectors for thresholds and constant for ploidy
			 for(i=0;i<1;i++)//
			 {
			 	  THR[i]=0;
			 	  msD[i]=0;
			 	  thr_sk[i]=0.0;
		     }
		     plo=0;// initial ploidy
		     //initial bins
		     for(i=0;i<4;i++)//
			 {
			 	   bins_ploid[i]=0;
			 }

		// scan thr file	and assign current precomputed thresholds mu sd ploidy for filtering
    while( (n = fscanf(thrFile,"%d %d %d %d %d", &thR, &thA,&muD, &stdD,&ploid)) >= 0)
         // until the end of the input threshold  file
    {
		      if (n != Nthr) // incorrect format input
		      {
		      printf ("corrupted input thrFile format\n");
		      return -1;
		      }
		THR[0]=thR;
		THR[1]=thA;
		msD[0]=muD;
        msD[1]=stdD;
        plo=ploid;
	}
  printf("ploidy= %d\n", plo);

     //===============================get precomputed bins for skewness dep on ploidy
     get_bins_ploidy(bins_ploid,plo);//

     for(i=0;i<4;i++){
     printf("ploid_bin= %d\n", bins_ploid[i]);

     }
    // mode 1: mix<25% diploid
	 		             n25=bins_ploid[0];
	 		             n49=bins_ploid[1];
	 		             n51=bins_ploid[2];
	 		             n75=bins_ploid[3];

	get_skew_ploidy	(thr_sk, plo);
      for(i=0;i<2;i++){
	      printf("ploid_skew= %.2f\n", thr_sk[i]);

      }
       // give value for thr25=    0.55 for human  diploid; 1.0, 0.8 for pathogens-haploid
	                 thr25=thr_sk[0];
	                 thr75=thr_sk[1];// complementary peak's (100-mix) skewness

	  //6 fields of input2 file	// no PV4: do it to compute AAF=AF


    while( (n = fscanf(extract_vcfFile,"%d,%d,%d,%d,%d,%d", &pos, &D, &DP4[0], &DP4[1], &DP4[2], &DP4[3])) >= 0 && firstSixAreOK == 1 && canWriteMix == 1)

    {
        if( n != NRID )     // incorrect format
        {
            firstSixAreOK = 0;
            break;
        }
           count_beforeF++;
           sumDbefore=sumDbefore+D;

         // f1 Josie : filter for DP4 separately for ref and alternative alleles+ filter abnormal Depth

           if ( DP4[0] >= THR[0] && DP4[1]>= THR[0] && DP4[2]>= THR[1] && DP4[3]>= THR[1] && D < (msD[0]+msD[1]))// 0 ref alt are ok
           {
                AF = GetAf(DP4);
                histAF[AF-1]++;

		        // for EACH POSITION, compute depths after thrRA filtering:
                       // average run depth after Q30, sum(DP4) across pos,=  sumDP4
                        count_afterF++;
				        sDP4=DP4[0]+DP4[1]+DP4[2]+DP4[3];
				        sumDP4=sumDP4+sDP4;
			}// end min_depth filter of DP4

			// -------------output stupidly deep positions if (D >= (msD[0]+msD[1])) { fill the file

       // removing space symbols before carriage return
        do  {
            n = fgetc (extract_vcfFile);
        }
        while ((char)n != '\n');

    } //END of outer while loop across vcfq file   AF is computed!

    //====================== output2    write AF
      for(i=0;i<101;i++)
      {
		  if( fprintf(distribFile,"%d %d\n",i, histAF[i]) <= 0 ) {
          canWriteDistrib = 0;
	      }

       }//  for i=100 cycle
//=======================================analysis of AF distribution to find mixture sample

        sk25=GetSkew25(histAF,n25,n49);
         sk75=GetSkew75(histAF,n51,n75);

//===================================mode 1:
// find  what percentage AAF gives first max peak among first 25% AAF
         mix25=GetMaxPeak25(histAF,n25);//peaks in first 25% AAF
         mixc25=GetPeakComplement25(histAF,n75);// complement to 25% peaks (100-mix25)

 //-----==================================mode2
		  mix25_50=GetMaxPeak25_50(histAF,n25,n49);//peaks in (25,49) percentage  AAF
		  mixc25_50=GetPeakComplement25_50(histAF,n51,n75);


// stats after filtering :compute average depth :default and after filtering

          avDP4=ceil(sumDP4/count_afterF);//actual average cov  after Filt
	      mix25_left=(100.0*count_afterF/count_beforeF);// percentage of filteerd out regions ( ? do we need?)

// make decision about on mixture percentage: update it to zero if needed
         // mix25_update=MakeDecisionMix(mix25, sk25, thr25,sk75,thr75);



  // confidences of mixture estimating

           conf25=Confidence25(sk25,sk75,mix25,mixc25, thr25, thr75, avDP4);
           conf25_50=Confidence25_50(sk25,sk75,mix25_50,mixc25_50, thr25, thr75, avDP4);

           printf("Results:\n");
           //printf("updated percentage mixture= %d\n", mix25_update);
           printf("possible low percentage mixture= %d\n", mix25);
           //printf("possible complement low percentage mixture= %d\n", mixc25);

           //if (( mix25-100+mixc25)<dist)
            //   conf25=conf25*1.1;

           printf("confidence of low freq mixture = %.4f\n", conf25);

           printf("possible high percentage mixture= %d\n", mix25_50);
           printf("confidence of high freq mixture = %.4f\n", conf25_50);
            //printf("possible complement high percentage mixture= %d\n", mixc25_50);

           //printf("Confidence of non-zero mixture = %.2f\n", conf);
	       //printf("skewness of first 25 mix25 AF = %.2f\n", sk25);
	       //printf("skewness of last 25 mix25 AF = %.2f\n", sk75);
	      // printf("likelihood of last 25 mix25 AF = %.2f\n", likely75);

	   // ======================MAIN output: mixFile
	  //if (fprintf(mixFile,"mix_percentage=%d\npossible_mix_percentage= %d\nhigh_freq_possible_mix_percentage= %d\nconfidence_nonzeroMix=%.2f\nAvActDepth=%d\nmin_depthR=%d\nmin_depthA=%d\nskew25=%.2f\nskew75=%.2f\n", mix25_update,mix25,mix25_50,conf,avDP4,thR,thA,sk25,sk75)<= 0 )
 	  //{
	  //             canWriteMix = 0;
	                //return -1;
      //}
       if (fprintf(mixFile,"mix lowfreq=%d\nconfidence low freq= %.4f\nmix high freq= %d\nconfidence high freq=%.4f\nAvActDepth=%d\nmin_depthR=%d\nmin_depthA=%d\nskew25=%.2f\nskew75=%.2f\n", mix25,conf25,mix25_50,conf25_50,avDP4,thR,thA,sk25,sk75)<= 0 )
	   	  {
	  	                canWriteMix = 0;
	  	                //return -1;
      }

    fclose(distribFile);
    fclose(extract_vcfFile);
    fclose(mixFile);

    if( firstSixAreOK  == 0 ||
        canWriteMix == 0 || canWriteDistrib ==0) {
        printf ("Error during execution. Details: \n");
        printf ("\tfirstSixAreOK %d\n",  firstSixAreOK);
        printf ("\tcanWriteMix %d\n",        canWriteMix);
        printf ("\tcanWriteDistrib %d\n",      canWriteDistrib);
        printf ("Execution aborted\n");
        return -1;
    }

    printf(" get mixture is performed.\n");
    return 0;
}//main

/***************************************************** Functions******************************/

 ////////////////////////////////////////////////////
// Calculates AF=alternative allele frequency from data=DP4 counts
// input - array (4-vector) of DP4 counts
//output -one float number from (0,1)
////////////////////////////////////////////////////
int GetAf (int data[])
{
    float Isum,AF1;
    int i,AF;

    Isum = data[0];
    for (i=1; i<4; i++)//there are four conts for each base
    {
        Isum += data[i];
    }
    AF1=(data[2]+data[3])/Isum;
    AF=(long) (100*AF1+0.5);
    return AF;
}
// GetMu.c calculates mu for all non-zero elements of 100-vector
float GetMu( int data[])
 {
	 float mu,Isum;
	 int i;
     int cnz=0;

     Isum = data[0];
     for (i=1; i<101; i++)
	 {
		 if (data[i] > 0)
		 {
		 cnz++;
	     Isum += data[i];
	     }
      }
      mu=Isum/cnz;

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
//==================================================================
//----------------------max from the 25-49
//----------------------max from the ( 25, 45)
 int GetMaxPeak25_50(int data[], int n25, int n49)
 {

  // input data=AAF
	 int Imax;
     int i,mix25_50,k;
     int peak[n25], mix25_peak[n25];

     // initialise mix25_50, in case no peaks exist

     mix25_50=0.0;

// find max peak among peaks, first 25 percentage only, data is maf(AAF) vector

//what if no peaks?
     peak[0]=0;

     // find peaks
     k=0;
     for (i=n25; i<n49; i++)
     {
		 if ( ((data[i]-data[i-1]) >0 ) && ((data[i+1]-data[i]) < 0))// && ( Imax < data[i]))//here peak in [i]
         {

             peak[k] = data[i];
             mix25_peak[k]=i;
             k=k+1;

         }
	 }
// found max peak/s if any exist

  if (peak[0]>0)
  {
	 Imax = peak[1];

	 	      for (i=1; i<k; i++)
	 	      {
	 	          if( Imax < peak[i])
	 	              Imax = peak[i];
	 	              //mix25_50=i;
	           }

       printf("Mode 2:max peak 25_50 = %d\n", Imax);

   // find percentage giving Imax, first one if several are max
       for (i=0; i<k; i++)
  	   {
  	          if( peak[i]==Imax)
  	            { mix25_50=mix25_peak[i];

  	              break;}
       }
   }// if (peak[1]>0)-exist any peak <25%
      printf("Mode 2: percentage giving this max = %d\n", mix25_50);
       return mix25_50;

}//==================================================================
//----------------------max from the first 25
 int GetMaxPeak25(int data[], int n25)
 {

  // input data=AAF
	 int Imax;
     int i,mix25,k;
     int peak[n25], mix25_peak[n25];

     // initialise mix25, in case no peaks exist

     mix25=0.0;

// find max peak among peaks, first 25 percentage only, data is maf(AAF) vector

//what if no peaks?
     peak[0]=0;

     // find peaks
     k=0;
     for (i=1; i<n25; i++)
     {
		 if ( ((data[i]-data[i-1]) >0 ) && ((data[i+1]-data[i]) < 0))// && ( Imax < data[i]))//here peak in [i]
         {

             peak[k] = data[i];
             mix25_peak[k]=i;
             k=k+1;

         }
	 }
// find max peak/s if any exist

  if (peak[0]>0)
  {
	 Imax = peak[1];

	 	      for (i=1; i<k; i++)
	 	      {
	 	          if( Imax < peak[i])
	 	              Imax = peak[i];
	 	              //mix25=i;
	           }

       printf("Mode 1: max peak25 = %d\n", Imax);

   // find percentage giving Imax, first one if several are max
       for (i=0; i<k; i++)
  	   {
  	          if( peak[i]==Imax)
  	            { mix25=mix25_peak[i];

  	              break;}
       }
   }// if (peak[1]>0)-exist any peak <25%
      printf("MOde1: percentage giving this max = %d\n", mix25);
       return mix25;

}
//====================================================================================

//===================================
float   GetSkew25( int data[],int n25, int n49)
 {

	 float sk25;
     int i,k;
     int less25;
     int bw25_50;
           // sum up from 0 to n25 %
           less25=data[0];
	 	   for(i=0;i<n25;i++)
	 	   {
	 	   		less25+=data[i];
	 	   }

	 	   bw25_50=data[n25];
	 	   for(i=n25;i<n49;i++)
	 	   {
	 	   	    bw25_50+=data[i];
	 	   }

	   //compute 25% AF skewness, sk25

	 	   sk25=0.00;// if no bases for bw25_50

	 	   if (bw25_50 >0 )
	 	   {
	 	   sk25=1.0*less25/bw25_50;// if proportion of less 25% AF is high (more 1?), more likely to be a mixture
	        }
	  return sk25;
}
//  ---------------------float   GetSkew( int data[],int n51, int n75);
float   GetSkew75( int data[],int n51, int n75)
 {

	 float sk75;
     int i,k;
     int more75;
     int bw75_50;
     int last;

     last=5;
           // sum up from n75 to almost 100%
           more75=data[n75];
	 	   for(i=n75;i<(100-last);i++)
	 	   {
	 	   		more75+=data[i];
	 	   }
           // sum 50 to 75
	 	   bw75_50=data[n51];
	 	   for(i=n51;i<n75;i++)
	 	   {
	 	   	    bw75_50+=data[i];
	 	   }

	   //compute 75% AF skewness, sk75

	 	   sk75=0.00;// if no bases for bw75_50

	 	   if (bw75_50 >0 )
	 	   {
	 	   sk75=1.0*more75/bw75_50;// if proportion of less 75% AF is high (more 1?), more likely to be a mixture
	        }
	  return sk75;
}


// make decision if mix25 noise or mixture: update it to zero if needed
//if mix25 is more than n25 (here =25), it is likely to be diploid allele noise  :  17 Nov+ amount AF<25%  skewed (too high)
//int MakeDecisionMix(int mix25, float sk25, float thr25)
int MakeDecisionMix(int mix25, float sk25, float thr25,float sk75,float thr75)// update percentage mix

{
       int mix25_update;
       mix25_update=0.0;
       if ((sk25 >= thr25) && (sk75 >= thr75))// thr of mixture=mix25 ;thr25 was 0.95 for path, 0.62 for human
       {
		   mix25_update=mix25;
		   //printf("percentage mixture= %d\n", mix25);
	   }


        return mix25_update;
}



//==========================================fill the thr bins for skewness
void get_bins_ploidy	(int bins_ploid[], int plo)//
{
	if (plo>=2){

			  bins_ploid[0]=24;
			  bins_ploid[1]=45;
			  bins_ploid[2]=55;
			  bins_ploid[3]=77;

			  // mode 1: mix<25% diploid
			              //n25=26;
			             // n49=46;
			             // n51=60;
			             // n75=83;


	  }

	  if (plo<2){

	  			  bins_ploid[0]=26;
	  			  bins_ploid[1]=50;
	  			  bins_ploid[2]=51;
	  			  bins_ploid[3]=75;

	  			  // mode 1: mix<25% haploid
	  			              //n25=26;
	  			             // n49=50;
	  			             // n51=51;
	  			             // n75=79;


	  }
}
//--------------------skewness dep on ploidy
void get_skew_ploidy	(float thr_sk[], int plo)// // how to consider skew25, skew75 depending on ploidy
{
        if (plo<2)
        {

	   	  			  thr_sk[0]=0.8;
	   	  			  thr_sk[1]=0.6;
	    }

  // give value for thr25=    0.55 for human  diploid; 1.0, 0.8 for pathogens-haploid
              //thr25=0.25;
              //thr75=0.5;//complementary skewness threshold


	   	if (plo>=2)
	   	{

	   	  			  thr_sk[0]=0.5;
	   	  			  thr_sk[1]=0.4;
	    }


}
//====================================
 int GetPeakComplement25(int data[], int n75)
 {

  // input data=AAF
	 int Imax;
     int i,mixc25,k,n;
     int mixcc25;//100-mix
     int peak[25], perc_peak[25];

     // initialise mixc25, in case no peaks exist

     mixc25=0.0;
     mixcc25=0.0;

// find max peak among peaks, first 25 percent only, data is maf(AAF) vector
     n=6;// not look at last n bins- noise for 100% SNPs
//what if no peaks?
     peak[0]=0;

     // find peaks
     k=0;
     for (i=n75; i<100-n; i++)
     {
		 if ( ((data[i]-data[i-1]) >0 ) && ((data[i+1]-data[i]) < 0))// && ( Imax < data[i]))//here peak in [i]
         {

             peak[k] = data[i];
             perc_peak[k]=i;
             k=k+1;

         }
	 }
// find max peak/s if any exist

  if (peak[0]>0)
  {
	 Imax = peak[1];

	 	      for (i=1; i<k; i++)
	 	      {
	 	          if( Imax < peak[i])
	 	              Imax = peak[i];
	 	              //mixc25=i;
	           }

       printf("Mode 1 com: max peak complement25 = %d\n", Imax);

   // find percent giving Imax, first one if several are max
       for (i=0; i<k; i++)
  	   {
  	          if( peak[i]==Imax)
  	            { mixc25=perc_peak[i];

  	              break;}
       }
   }// if (peak[1]>0)-exist any peak <25%

   mixcc25=100-mixc25;
      printf("Mode 1 com:percentage giving this complement max = %d\n", mixcc25);
       return mixcc25;

}
//===========================
//==================================================================
//----------------------max peak from the 50- 75
 int GetPeakComplement25_50(int data[], int n51, int n75)
 {

  // input data=AAF
	 int Imax;
     int i,mixc25_50,k,n;
     int mixcc25_50;// 100-mixc
     int peak[25], perc_peak[25];

     // initialise mixc25_50, in case no peaks exist

     mixc25_50=0.0;
     mixcc25_50=0.0;

// find max peak among peaks, first 25 percent only, data is maf(AAF) vector

//what if no peaks?
     peak[0]=0;

     // find peaks
     k=0;
     for (i=n51; i<n75; i++)
     {
		 if ( ((data[i]-data[i-1]) >0 ) && ((data[i+1]-data[i]) < 0))// && ( Imax < data[i]))//here peak in [i]
         {

             peak[k] = data[i];
             perc_peak[k]=i;
             k=k+1;

         }
	 }
// find max peak/s if any exist

  if (peak[0]>0)
  {
	 Imax = peak[1];

	 	      for (i=1; i<k; i++)
	 	      {
	 	          if( Imax < peak[i])
	 	              Imax = peak[i];
	 	              //mixc25_50=i;
	           }

       printf("Mode 2 com: max peak complement 25_50 = %d\n", Imax);

   // find percent giving Imax, first one if several are max
       for (i=0; i<k; i++)
  	   {
  	          if( peak[i]==Imax)
  	            { mixc25_50=perc_peak[i];

  	              break;}
       }
   }// if (peak[1]>0)-exist any peak <25%
      mixcc25_50=100-mixc25_50;
      printf("Mode 2 com: percentage giving this complement max = %d\n", mixcc25_50);
       return mixcc25_50;

}
// ============================================MODE1
// confidence of mixture <25%-mode1
float Confidence25(float sk25,float sk75,int mix25,int mixc25, float thr25, float thr75, float avDP4)
{
//likelihoods


    float conf,lik1,likc1;
    float E1,EC1;
    int peak_dist;

    peak_dist=8;

    //--------------------
 printf("Mode1:\n ");
    conf=1.0;//initial

    E1=0.0;
    if (mix25 >0 ) E1=1.0;// peak <25 does Exist
    printf("exist mix25 = %.2f\n", E1);

    EC1=0.0;
    if (abs(mixc25-mix25) < peak_dist )    // close enough
    EC1=1.0;
    //if (mixc25 >0 ) EC1=1.0;
    printf(" exist complem mix25 = %.2f\n", EC1);

    lik1=sk25*sk25;
    if (sk25 >= thr25)
    {lik1=sk25;}


    printf(" likely mix25 = %.2f\n", lik1);

    likc1=sk75*sk75;
    if (sk75 >= thr75)
    {likc1=sk75;}
   printf(" likely complem mix25 = %.2f\n", likc1);

      // larger any how!
          conf=(1-1/avDP4)*(E1+EC1)*0.5*lik1*likc1;

          // if (abs(mixc25-mix25) < peak_dist )// close enough
          //     conf=conf*1.1;

            if (conf>1)
            {conf=1.0;}

     return conf;
}



//======================
// confidence of mixture 25%-50%  mode2
float Confidence25_50(float sk25,float sk75,int mix25_50,int mixc25_50, float thr25, float thr75, float avDP4)
{
//likelihoods


    float conf,lik1,likc1;
    float E1,EC1;
    int peak_dist;

    peak_dist=8;

    //--------------------
 printf("Mode2:\n ");
    conf=1.0;//initial

    E1=0.0;
    if (mix25_50 >0 ) E1=1.0;// peak <25 does Exist
    printf("exist mix25_50 = %.2f\n", E1);

    EC1=0.0;
   // if (mixc25_50 >0 ) EC1=1.0;

    if (abs(mixc25_50-mix25_50) < peak_dist )    // close enough
	    EC1=1.0;

    printf(" exist complem mix25_50 = %.2f\n", EC1);

    lik1=sk25*sk25;//assume sk25<1
    if (sk25 < thr25)// humph 25-20 is larger than humph25=sk25!
    {lik1=sk25;}
    printf(" likely mix25_50 = %.2f\n", lik1);

    likc1=sk75*sk75;
    if (sk75 < thr75)
    {likc1=sk75;}
   printf(" likely complem mix25_50 = %.2f\n", likc1);

      // larger any how!
          conf=(1-1/avDP4)*(E1+EC1)*0.5*lik1*likc1;

            if (conf>1)
            {conf=1.0;}

     return conf;
}
//-------------------------
int let2int(char letter)
{

	int result;

	//initiate
	result=0;

	if (letter=='a' || letter=='A')
	    result=1;
	if (letter=='c' || letter=='C')
	    result=2;
	if (letter=='g' || letter=='G')
	    result=3;
	if (letter=='t' || letter=='T')
	    result=4;

	    return result;
	}
