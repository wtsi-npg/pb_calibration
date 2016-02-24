/*  File: get_thr_ploidy.c // filters 1 1 vcfq file and computes av std Depth  depending on Depth of vcfq file

 * Authors: designed by Irina Abnizova (ia1) and edited by Steve Leonard(srl)
 *
  Last edited:
  18 Jan 2016 - real std
  11 June: one thr file. output is e.g.  2 2 41 113
  9 June 2014 :n1 from pipe of three:filters 1 1 vcfq file and computes av std Depth  depending on Depth of vcfq file
  28 May : not output filtered vcfq
  22 May  get_filters 1 1 .vcfq  RA_ms.filters
  21 May  filters bad varians: abnormal Depth
  30 April: adds Filteref 1,1 vcfq + distr_1_1
  29 April
  // computes suitable thresholds bases on avDP4=actual (not cov before as earlier!!!

  3 Feb-put thrR thrA into file thrD ;29 Jan 2014  - get autom thr, starting with default ones 1,5
  still old fields of vcf input

  *-------------------------------------------------------------------
 * Description: computes computes threshold file depending on Depth of vcfq file

 * Exported functions:
 * HISTORY:

computes minimum depth threshold, both for Ref and Alternative alleles fwd rev (total min depth will be
4*MIN_DEPTH


input: MIN_DEPTH_R(integer)=1[default]  MIN_DEPTH_A(integer)=1[default]   extracted=vcfq,

output1:
thrR  thrA
3      1
output2:  .distAF_1_1
output3:  .vcfqF_1_1  filtered 1,1
usage:./get_thrRA MIN_DEPTH_R(integer) MIN_DEPTH_A(integer) name.vcf  ra_i_j.thr(Depth)  name.distr name.vcfqF_1_1

*/


#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

// ******** predefined user constants,

#define NRID 6      // Nunber of variant ids=fields: pos, Depth,DP4,(no PV4+maf) in extracted_vcf files
#define NPAR 7      // Number of PARameters and arguments


// *********    Global variables min_depths for Ref and Alternative alleles
int thR = 0;
int thA = 0;

// ******** declarations of  functions

int   GetAf (int []);// computes AF percent for a variant position
int   GetMax (int []);// percent giving max AAF
int ploidy (int ind);// defines ploidy from results of GetMax ploid=2 if ~50%
int GetMax50	(int data[]);// local max around 50%

int main    (int argc, char *argcv[]) {

    // flags
    int firstSixAreOK  = 1;
    int canWriteDistrib      = 1;
    int canWriteFilteredVCF      = 1;

    FILE *extract_vcfFile, *thrFile, *distribFile, *filtered_vcfFile; //file handles

    int n,i,ind,ploid, ind50;//
    int count_afterF=0;// after filtering with min depths
    int count_beforeF=0;// before filtering with min depths
    int DP4[4];//
    int D,pos;// Depth, genome pos of error
    int AF;

     //arrays
	int histAF[100];// to store mafs percentages


    //stats
    int sumDP4=0,sumDP4_2=0;//across all pos for mu std
    int sumDbefore=0;
    float avDbefore;
	float Q30frac;
	float sDP4;
	float perc_left;

    //to compute
    int avDP4,stdDP4,stdDP4_1,avDP41;// of actual sDP4 depth after 1 1 filtering
    int thRR, thAA;// recommended based on actual depth avDP4 and stdDP4
    float ratio;

      if(argc < NPAR)//five input_output files are submitted
    {
        printf("not enough of parms |input_output files\n");
        printf("usage:./get_thr thrR thrA input2 output1 output2 output3\n");
        return -1;
    }

    if (sscanf(argcv[1],"%d",&thR) == EOF) {
                printf("Failed to convert the first argument %s to integer\n", argcv[1]);
            }


     if (sscanf(argcv[2],"%d",&thA) == EOF) {
                printf("Failed to convert the second argument %s to integer\n", argcv[2]);
            }


    extract_vcfFile=fopen(argcv[3],"r");
    if (extract_vcfFile == NULL) {
      printf("cannot open first input _vcf.txt file %s\n", argcv[3]);
      return -1;
    }
 // two outputs
     thrFile=fopen(argcv[4],"w");
            if (thrFile == NULL) {
            printf("cannot open first output file thrD %s\n", argcv[4]);
            return -1;
    }

     distribFile=fopen(argcv[5],"w");
           if (distribFile == NULL) {
           printf("cannot open second output file distAF %s\n", argcv[5]);
           return -1;
	 }

     filtered_vcfFile=fopen(argcv[6],"w");
		    if (filtered_vcfFile == NULL) {
		    printf("cannot open second output file filteredVCF %s\n", argcv[6]);
           return -1;
     }

        printf("get_thr_ploidy\n");

         // initiate zero vector for histograme
		   			 for(i=0;i<100;i++)//
		   	         {
		   	              histAF[i]=0;
				     }


    //6 fields of input file,CSV reading
    while( (n = fscanf(extract_vcfFile,"%d,%d,%d,%d,%d,%d", &pos, &D, &DP4[0], &DP4[1], &DP4[2], &DP4[3])) >= 0 && firstSixAreOK == 1)// && canWriteAF == 1)
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

           if ( DP4[0] >= thR && DP4[1]>= thR && DP4[2]>= thA && DP4[3]>= thA)// && D < 10)// 0 ref are ok
           {

			     // compute AAF histo
				AF = GetAf(DP4);
                histAF[AF-1]++;

			     // for EACH POSITION, compute actual depths sDP4 after thrRA filtering:

                                // average run depth after Q25, sum(DP4) across pos,=  sumDP4
                                        count_afterF++;
                                        sDP4=DP4[0]+DP4[1]+DP4[2]+DP4[3];// for one position
                                        sumDP4=sumDP4+sDP4;//for av
                                        sumDP4_2=sumDP4_2+sDP4*sDP4;//for std

                                        Q30frac=sDP4/D;// how actual  quality Depth differs from default

             if( fprintf(filtered_vcfFile,"%d,%d,%d,%d,%d,%d\n",pos, D, DP4[0], DP4[1],DP4[2], DP4[3]) <= 0 ) {
			          canWriteDistrib = 0;
			          //"%d,%d,%d,%d,%d,%d", pos, D, DP4[0], DP4[1],DP4[2], DP4[3]))
	         }

            }// end filter DP4

        // removing space symbols before carriage return
        do  {
            n = fgetc (extract_vcfFile);
        }
        while ((char)n != '\n');

    } //END of outer while loop for all vcfq file, distrAF is computed

    // check ploidy
    ind=GetMax(histAF);
    printf("percent of max AAF= %d\n", ind);

    ind50=GetMax50(histAF);
    printf("percent of max AAF around 50= %d\n", ind50);

     ploid=ploidy(ind);
     printf("ploidy= %d\n", ploid);
 //====================== output2    write AF
      for(i=0;i<101;i++)
      {
		  if( fprintf(distribFile,"%d %d\n",i, histAF[i]) <= 0 ) {
          canWriteDistrib = 0;
	      }

       }//  for i=100 cycle
    // compute average depth : default and actual after 1,1 and q25 filtering

      // stats after 1,1 filtering

       avDP4=ceil(sumDP4/count_afterF);//actual average cov  after Filt
       //stdDP4=ceil(sqrt(sumDP4_2/count_afterF));// why did I approximate this like that?...
       stdDP4=ceil(sqrt(sumDP4_2/count_afterF -avDP4*avDP4));//
       avDbefore=sumDbefore/count_beforeF;
       perc_left=(100.0*count_afterF/count_beforeF);
     printf("average Depth after filtering= %d\n", avDP4);
     printf("std Depth appr= %d\n", stdDP4);
     printf("percent data left after q25 & min_depth filtering thrR thrA = %.2f\n", perc_left);

 // adjust actual depth to Bad (deep covered) regions and recommend thr RA

        avDP41=avDP4;
        ratio=(sqrt(sumDP4_2/count_afterF))/(sumDP4/count_afterF);//stdDP4/avDP4;
        printf("cv actual depth = %.2f\n", ratio);
        if ( ratio > 1.5) //bad regions
           {avDP41=ceil(0.8*avDP4);
           printf("moderately deep covered regions here\n");
           printf("adjusted(avDP4) = %d\n", avDP41);
           }
         if (ratio > 3) //very bad regions
		            {avDP41=ceil(0.7*avDP4);
		            printf("extremely deep covered regions here\n");
		            printf("adjusted(avDP4) = %d\n", avDP41);
           }

           // recommendation for thrRR AA
           if ( avDP41 > 0 && avDP41 < 4 )
		     	     { thRR = 100;thAA=100;//stupidly high?
		     	     printf("better not to bother to look for mixture: too low actual coverage=%d\n",avDP4);
		     	     }
		     	     if ( avDP41 >= 4 && avDP41 < 10  )
		     		      { thRR = 1;thAA=1;
		     	     }
		     	     if ( avDP41 >= 10 && avDP41 < 70  )
		     		 { thRR = 2;thAA=2;
		     	     }
		     	     if ( avDP41 >= 70 && avDP41 < 90  )
		     		 { thRR = 3;thAA=2;
		     	     }
		     	     if ( avDP41 >= 90 )
		     		 { thRR = 3;thAA=3;
		  	     }

 // fill in the output files

      //====================== output1   write thrs, avD, stdD and ploidy
        fprintf(thrFile,"%d %d %d %d %d\n",thRR,thAA,avDP4,stdDP4,ploid);

    fclose(extract_vcfFile);
    fclose(thrFile);
    fclose(distribFile);
    fclose(filtered_vcfFile);

       // checking write/read
        if( firstSixAreOK  == 0)// || canWriteFilteredVCF == 0)
        {
		        printf ("Error during execution. Details: \n");
		        printf ("\tfirstSixAreOK %d\n",  firstSixAreOK);
		        printf ("\tcanWriteDistrib %d\n",      canWriteDistrib);
		        printf ("\tcanWriteFilteredVCF %d\n",   canWriteFilteredVCF);
		        printf ("Execution aborted\n");
		        return -1;
        }

    printf(" done get_thr_ploidy.\n");
    return 0;
}//main

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
    //(long) (sig+0.5)
    return AF;
}
//--------------fins max value and its bin/index; returns bin; data is histo

int GetMax	(int data[])
{
	int Imax;
	int i,ind;

	Imax = data[0];
	for (i=1; i<80; i++)// for histo: not latest bins (for bad reference ~100)
	{
		if( Imax < data[i])
		{
			Imax = data[i];
			ind=i;// bin/index for max value
		}
	}
	return ind;
}

//-----------------------identify ploidy
int ploidy (int ind)
{

	int ploid;

	ploid=1;

	if ( ind >40 && ind < 60)
	ploid=2;

	return ploid;
}

//-----------------------double check local max around 50+-10

int GetMax50	(int data[])
{
	int Imax;
	int i,ind;

	Imax = data[0];
	for (i=40; i<60; i++)// for histo: not latest bins (for bad reference ~100)
	{
		if( Imax < data[i])
		{
			Imax = data[i];
			ind=i;// bin/index for max value
		}
	}
	return ind;
}
