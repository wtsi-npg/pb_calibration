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

computes AF (alternative to Ref allele frequency) for each error base call from extracted vcf

applies minimum depth threshold, both for Ref and Alternative alleles fwd rev (total min depth will be
4*MIN_DEPTH

input1 thrfile  (thRR thAA  mu  std  ploidy )
                 3     2    82  110     2

input2: extracted vcf: pos, D , DP4 (DP4 is four fields)


output1: mix=
mixture mode1 likelihood1
12               0.87
mixture mode2 likelihood2
12               0.87
mixture mode3 likelihood3
12               0.87

output2: .distrAF=histogram of AF with significance for each bin
% count
8 0
9 0
10 385
11 780


usage:./get_mixture name.thr extracted.vcf name.mix name.distrAF

*/

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

// ******** predefined user constants,

#define NRID 6      // Number of variant ids=fields: pos, chromosome, Depth,DP4,(no PV4+maf) in esigtracted_vcf files
#define NPAR 5      // masig Number of PARameters and arguments
#define Nthr 5     //  Number of thresholds and muD1 stdD1 12 June + ploidy=1 (haploid) (or 2-diploid)

// ******** declarations of  functions
int     GetAf (int []);// computes AF percentage for a variant position
float   GetMu( int data[]);
float   GetStd( int data[]);

void skewnesses (int data[], float sk[], float skc[], int st[], int stc[]);
void MixMode(int mix[], int st[], int data[], int dist);
void MixModeComplement(int mixc[], int stc[], int data[], int dist);
int GetMaxPeakInterval(int data[], int n1, int n2);//------max peak from the interval
void Likely (float Lik[], float sk[],float skc[], int mixc[],float avDP4);

//=====================================================MAIN

int main    (int argc, char *argcv[]) {

    // ---------------------flags
    int firstSixAreOK  = 1;
    int canWriteMix        = 1;
    int canWriteDistrib      = 1;

    FILE *thrFile, *extract_vcfFile, *mixFile, *distribFile ; //file handles

        int n,i;

      // ---------- params from input1 rams(now four+1 values thrR thrA muD stdD ploidy)
		 int muD,stdD;
		 int thR,thA, ploid,plo;


    //---------------params from input2 vcfq
	    int DP4[4];//
	    int D,pos;// Depth, genome pos of error

	    //----------stats to compute
	    int avDP4;// of actual sDP4 depth after RA & bad regions filtering
	    int sumDP4=0;//across all pos for mu std
	    int sumDbefore=0;

		float sDP4;
		float filtered_left;

       // ----------------settings for skewnesses and their likelihoods: hardcoded within functions

       int dist;//distance , bin width for hist to compute skewness

      //---------------results to compute:
       int count_afterF=0;// after filtering with min depths
	   int count_beforeF=0;// before filtering with min depths

       int AF;

       //---------------------------arrays
       int histAF[100];// to store mafs percentages
       int st[4], stc [4];// starts of summing intervals for histo
       float sk[3],skc[3];// three skewnesses for mode0,1,2
  	   int mix[3],mixc[3]; //mixtures three modes and their complements
	   float Lik[3];// likelihood or confidences three modes

	//-----------------------------------------open files for read/write

    if(argc < NPAR)//four input_output files are submitted
	{
	        fprintf(stderr, "not enough of parms |input_output files\n");
	        fprintf(stderr, "usage:./get_mixture input1.thr input2.vcfq output1.mix output2.distr\n");
	        exit(EXIT_FAILURE);
    }

    thrFile=fopen(argcv[1],"r");
			    if (thrFile == NULL) {
			    fprintf(stderr, "cannot open first input _threshold file %s: %s\n", argcv[1], strerror(errno));
			    exit(EXIT_FAILURE);
    }


    extract_vcfFile=fopen(argcv[2],"r");
	          if (extract_vcfFile == NULL) {
	          fprintf(stderr, "cannot open first input _.vcfq file %s: %s\n", argcv[2], strerror(errno));
	          exit(EXIT_FAILURE);
    }

    mixFile=fopen(argcv[3],"w");
	        if (mixFile == NULL) {
	        fprintf(stderr, "cannot open first output mix_conf file %s for writing: %s\n", argcv[3], strerror(errno));
	        exit(EXIT_FAILURE);
    }

    distribFile=fopen(argcv[4],"w");
            if (distribFile == NULL) {
            fprintf(stderr, "cannot open second output distAF file %s for writing: %s\n", argcv[4], strerror(errno));
            exit(EXIT_FAILURE);
    }


         fprintf(stderr, "get_mixture: parameters\n");


    // ===SETTINGS     parameters for skewness measurements and peak consideretion: if Diploid then more noise around 50%

             plo=0;// initial ploidy
             dist=11;
             avDP4 = 0;//initial average depth

             // -----------initiate zero vector for histograme and its peaks
			 for(i=0;i<100;i++)//
	         {
	              histAF[i]=0;
	         }

		      //--------------starts for skews and complement skews bins
			  for(i=0;i<4;i++)//
			  {
			       st[i]=0;
			 	   stc[i]=0;
			  }
		      // ------------initial skwnesses, mix, conf for modes 0,1,2
              for(i=0;i<3;i++)//
			  {
			 	  sk[i]=0.0;
			 	  skc[i]=0.0;
			 	  mix[i]=0;
			 	  mixc[i]=0;
			 	  Lik[i]=0.0;
		      }

   //1. scan-read thr file	and assign current precomputed thresholds mu sd ploidy for filtering
   while( (n = fscanf(thrFile,"%d %d %d %d %d", &thR, &thA,&muD, &stdD,&ploid)) >= 0)
   {
  		      if (n != Nthr) // incorrect format input thr file
  		      {
              fprintf (stderr, "corrupted input thrFile format\n");
  		      exit(EXIT_FAILURE);
  		      }

        plo=ploid;
   }
        fprintf(stderr, "ploidy= %d\n", plo);


   //2. scan-read main vcf file:6 fields of input2 file	// no PV4: do it to compute AAF=AF

    while( (n = fscanf(extract_vcfFile,"%d,%d,%d,%d,%d,%d", &pos, &D, &DP4[0], &DP4[1], &DP4[2], &DP4[3])) >= 0 && firstSixAreOK == 1 && canWriteMix == 1)
    {
         if( n != NRID )     // incorrect format vcfq file
		 {
			  fprintf(stderr,"reading strange line= %d\n", n);
			  continue;//skips this line if not end
         }
              count_beforeF++;
              sumDbefore=sumDbefore+D;

         // filter for DP4 separately for ref and alternative alleles+ filter abnormal large Depth
           if ( DP4[0] >= thR && DP4[1]>= thR && DP4[2]>= thA && DP4[3]>= thA && D < (muD+2*stdD))// 0 ref alt are ok
           {
                AF = GetAf(DP4);
                histAF[AF-1]++;

		        // for EACH POSITION, compute depths after thrRA and stupid Depth filtering:
                       // average run depth after Q30, sum(DP4) across pos,=  sumDP4
                 count_afterF++;
				 sDP4=DP4[0]+DP4[1]+DP4[2]+DP4[3];
				 sumDP4=sumDP4+sDP4;
			}// end filters and AAF computing for one line DP4

			//? -------------output stupidly deep positions if (D >= (msD[0]+msD[1])) { fill the file

        // removing space symbols before carriage return
        do  {
            n = fgetc (extract_vcfFile);
        }
        while ((char)n != '\n');

   } //END of  while loop across vcfq file   AF is computed!


             fprintf(stderr,"count after filtering = %d\n",count_afterF);

  if ( count_afterF >200 )// non empty, large enough after filtering vcfq
  {
  //=======================================analysis of AF distribution to find mixture sample

         // stats after filtering :compute average depth :default and after filtering

          avDP4=ceil(sumDP4/count_afterF);//actual average cov  after Filt
	      filtered_left=(100.0*count_afterF/count_beforeF);// percentage of filteerd out regions ( ? do we need?)
          fprintf(stderr,"percentage left after filtering = %.2f\n",filtered_left);
          fprintf(stderr,"coverage after filtering = %d\n",avDP4);

          //skewnesses (histAF, sk, skc);
          skewnesses (histAF, sk, skc, st,stc);
          for (i=0;i<3;i++)
          fprintf(stderr,"sk = %.2f\n",sk[i]);
          for (i=0;i<3;i++)
          fprintf(stderr,"skc = %.2f\n",skc[i]);

          // mixtures per mode; confidences per mixture

           fprintf(stderr,"Results:\n");
           MixMode(mix, st, histAF,dist);
           MixModeComplement(mixc, stc, histAF,dist);
           Likely (Lik, sk,skc, mixc,avDP4);

	}// if not empty data

//==========================================================OUTPUTS
     //====================== output2    write AF
		      for(i=0;i<101;i++)
		      {
		   		  if( fprintf(distribFile,"%d %d\n",i, histAF[i]) <= 0 ) {
		             canWriteDistrib = 0;
		   	      }

              }//  for i=100 cycle

           //=================MAIN output1  mix file

	    if (fprintf(mixFile,"mix low freq=%d\nconfidence low freq= %.4f\nmix middle freq= %d\nconfidence middle freq=%.4f\nmix high freq= %d\nconfidence high freq=%.4f\nAvActDepth=%d\nmin_depthR=%d\nmin_depthA=%d\n", mix[0],Lik[0],mix[1],Lik[1],mix[2],Lik[2],avDP4,thR,thA)<= 0 )
           {
	  	  	  	                canWriteMix = 0;

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
        exit(EXIT_FAILURE);
    }

    fprintf(stderr," get mixture is performed.\n");
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

//1============================new skewnesses for three modes
void skewnesses (int data[], float sk[], float skc[], int st[], int stc[])
{

     int i,j,d,n1,n2;
     float s[4];//sums
     float sc[4];

     d=11;//

     //1.  starts, starts complementary
     st[0]=5;
     stc[0]=51;
     for (i=1;i<4;i++)
     {
		 st[i]=st[i-1]+d;
		 stc[i]=stc[i-1]+d;
		 //printf(" starts=%d\n",st[i]);
         //printf(" starts com=%d\n",stc[i]);
   }

     //2. sums between starts: 3 sums and 3 sums complementary

      for(i=0;i < 4;i++)
      {
         n1=st[i];
         n2=stc[i];
         s[i]=0;
         sc[i]=0;
         for (j=n1;j<n1+d;j++)
         {
	   		s[i]=s[i]+data[j];
	   		 //printf(" sums=%.2f\n",s[i]);
	   	 }
	   	 for (j=n2;j<n2+d;j++)
		 {
		 	sc[i]=sc[i]+data[j];
	   	 }
	   	 //printf(" sums=%.2f\n",s[i]);
	   	 //printf(" sums com=%.2f\n",sc[i]);
	  }

     //3.     // compute skewnesses

     for (i=0;i<3;i++)
     {
          if (s[i+1]>0) sk[i]=s[i]/s[i+1];
          if (sc[i]>0) skc[i]=sc[i+1]/sc[i];
     }


}

//2============================================

void MixMode(int mix[], int st[], int data[], int dist)

{
	 int i;
	 int n1,n2;
	 int m;

	 for (i=0;i< 3;i++)
	 {
		 n1=st[i];
		 n2=n1+dist;
	     m=GetMaxPeakInterval(data, n1, n2);//2.1 function
         mix[i]=m;
         printf(" mixes in modes=%d\n",mix[i]);
	 }

}
//2.1.1============================================

void MixModeComplement(int mixc[], int stc[], int data[], int dist)

{
	 int i;
	 int n1,n2;
	 int m;

	 for (i=0;i< 3;i++)
	 {
		 n1=stc[i+1];
		 n2=n1+dist;
	     m=GetMaxPeakInterval(data, n1, n2);//2.1 function
         mixc[i]=m;
         printf(" mixes in modes=%d\n",mixc[i]);
	 }

}
//2.1==================================================================
//----------------------max peak from the interval
 int GetMaxPeakInterval(int data[], int n1, int n2)
 {
 // input data=AAF
	 int max_peak;
     int i,mixx,k,j;
     int peak[100], perc_peak[100];
 // initialise mix, in case no peaks exist
     mixx=0.0;
// find max peak among peaks in the interval only, data is maf(AAF) vector

//what if no peaks?
     peak[0]=0;
     perc_peak[0]=0;

     // find peaks
     k=0;
     for (i=n1; i<n2; i++)
     {
		 if ( ((data[i]-data[i-1]) >=0 ) && ((data[i+1]-data[i]) < 0))// && ( max_peak < data[i]))//here peak in [i]
         {

             peak[k] = data[i];
             perc_peak[k]=i;
             //printf("peaks= %d\n", peak[k]);
             //printf("peaks percent= %d\n", perc_peak[k]);
             k=k+1;
         }
	 }
// found max peak/s if any exist

//if one peak
  if (k==1)
  {
	max_peak = peak[0];
	 mixx=perc_peak[0];
  // printf("Mode i:max peak  = %d\n", max_peak);
  // printf("Mode i: percentage giving this max = %d\n", mixx);
  }

  if (k>1)//more than one peak
  {
	 max_peak = peak[0];

	 	      for (j=1; j<k; j++)
	 	      {
	 	          if( max_peak < peak[j])
	 	              max_peak = peak[j];
	 	              //mix=i;
	           }
      // printf("Mode i:max peak  = %d\n", max_peak);

   // find percentage giving max_peak, first one if several are max
       for (j=0; j<k; j++)
  	   {
  	          if( peak[j]==max_peak)
  	            { mixx=perc_peak[j];

  	             break;
			    }
       }
        //printf("Mode i: percentage giving this max = %d\n", mixx);
   }// if (peak[1]>0)-exist any peak <25%

       return mixx;

}//==================================================================

//--confidences modes0,1,2
void Likely (float Lik[], float sk[],float skc[], int mixc[],float avDP4)
{
  int i;
  float l,lc,lp;

  for (i=0;i<3;i++)
  {
	  lp=0.0;
	  if (mixc[i]>0)
	  lp=skc[2-i];
	  l=sk[i];
	  lc=skc[2-i];
	  Lik[i]=(1-1/avDP4)*(l+lc+lp)/3;
	   printf("Mode i: likelihood of mix = %.2f\n", Lik[i]);

}
}






