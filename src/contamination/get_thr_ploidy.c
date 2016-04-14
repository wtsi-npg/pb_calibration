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
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>


// ******** predefined user constants,

#define NRID 6      // Number of variant ids=fields: pos,Depth,DP4,(no PV4+maf) in extracted_vcf files
#define NPAR 7      // Number of PARameters and arguments


// *********    Global variables min_depths for Ref and Alternative alleles
int thR = 0;
int thA = 0;

// ******** declarations of  functions

int   GetAf (int []);// computes AF percent for a variant position
int   GetMax (int []);// percent giving max AAF
int ploidy (int ind);// defines ploidy from results of GetMax ploid=2 if ~50%
int GetMax50	(int data[]);// local max around 50%
int let2int(char letter);


int main    (int argc, char *argcv[]) {

    // -----error flags
    int firstSixAreOK  = 1;
    int canWriteDistrib      = 1;
    int canWriteFilteredVCF      = 1;

    FILE *extract_vcfFile, *thrFile, *distribFile, *filtered_vcfFile; //file handles

    // ---------------to read vcfa/vcfq
    static const int line_size = 8192;//stupidly large
    char line[line_size];

    int k,i,ind,ploid, ind50;//
    int count_afterF=0;// after filtering with min depths
    int count_beforeF=0;// before filtering with min depths

    //----------from vcfa/vcfq input
    int DP4[4];//
    int pos,D;// Depth, genome pos of error

     //--------------arrays
	int histAF[100];// to store mafs percentages


    //-----------------stats
    int sumDP4=0,sumDP4_2=0;//across all pos for mu std
    int sumDbefore=0;
    float sDP4;
	float perc_left;

    //-----------to compute fot outputs

    int AF;
    int avDP4,stdDP4,avDP41;// of actual sDP4 depth after 1 1 filtering
    int thRR, thAA;// recommended based on actual depth avDP4 and stdDP4
    float ratio;


//--------------------------opening
    if(argc < NPAR)//five input_output files are submitted
    {
        fprintf(stderr, "not enough of parms |input_output files\n");
        fprintf(stderr, "usage:./get_thr thrR thrA input2 output1 output2 output3\n");
        exit(EXIT_FAILURE);
    }

    if (sscanf(argcv[1],"%d",&thR) == EOF) {
	     fprintf(stderr, "Failed to convert the first argument %s to integer\n", argcv[1]);
	}

    if (sscanf(argcv[2],"%d",&thA) == EOF) {
	      fprintf(stderr, "Failed to convert the second argument %s to integer\n", argcv[2]);
	}

     extract_vcfFile=fopen(argcv[3],"r");
	 if (extract_vcfFile == NULL) {
	     fprintf(stderr, "cannot open first input _vcf.txt file %s: %s\n", argcv[3], strerror(errno));
	     exit(EXIT_FAILURE);
    }

     // -----------------------three outputs

     thrFile=fopen(argcv[4],"w");
	 if (thrFile == NULL) {
	      fprintf(stderr, "cannot open first output file thrD %s: %s\n", argcv[4], strerror(errno));
	      exit(EXIT_FAILURE);
     }

     distribFile=fopen(argcv[5],"w");
	 if (distribFile == NULL) {
	       fprintf(stderr, "cannot open second output distAF file %s: %s\n", argcv[5], strerror(errno));
	       exit(EXIT_FAILURE);
	 }

     filtered_vcfFile=fopen(argcv[6],"w");
	 if (filtered_vcfFile == NULL) {
	 	   fprintf(stderr, "cannot open third output filteredVCF file %s: %s\n", argcv[6], strerror(errno));
	       exit(EXIT_FAILURE);
     }

            fprintf(stderr,"get_thr_ploidy\n");

       // initiate zero vector for histograme
		   			 for(i=0;i<100;i++)//
		   	         {
		   	              histAF[i]=0;
				     }
       // initialise values
                   thRR=0;
                   thAA=0;
                   avDP41=0;
                   avDP4=0;
                   stdDP4=0;
                   ploid=1;

     //read main vcfq/a file: 6 fields of input file,CSV reading
     while (fgets(line, line_size, extract_vcfFile))
     {
	// read what is in this string=line
      k = sscanf(line, "%d,%d,%d,%d,%d,%d", &pos,&D, &DP4[0], &DP4[1], &DP4[2], &DP4[3]); //  && firstSixAreOK == 1)

              if( k != NRID )     // incorrect format
	          {
	            printf("reading strange line= %d\n", k);
	            continue;//skips this line if not end
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

           //==============  output3 filtered 1 1 vcfq file
             if( fprintf(filtered_vcfFile,"%d,%d,%d,%d,%d,%d\n",pos,D, DP4[0], DP4[1],DP4[2], DP4[3]) <= 0 )
             {
			        canWriteDistrib = 0;
	         }

         }// end filter f1 DP4

    } //END of outer while loop for all vcfq file, distrAF is computed

         fprintf(stderr,"count of good lines total before filter = %d\n", count_beforeF);
         fprintf(stderr,"count of good lines total after filter = %d\n", count_afterF);

     //====================== output2    write AF
	    for(i=0;i<101;i++)
	    {
			  if( fprintf(distribFile,"%d %d\n",i, histAF[i]) <= 0 ) {
	          canWriteDistrib = 0;
		      }
        }//  for i=100 cycle

  // -----------------------if large enough vcfq file
 if ( count_afterF >200 )// non empty, large enough after filtering vcfq
 {


       // -----------  check ploidy detection
       ind=GetMax(histAF);// percent giving max peak in hist <80%
       fprintf(stderr,"percent of max AAF= %d\n", ind);

       ind50=GetMax50(histAF);//percent giving max peak in hist within (40,60)%
       fprintf(stderr,"percent of max AAF around 50= %d\n", ind50);//

       ploid=ploidy(ind);// automated detection of ploidy: 1 or 2
       fprintf(stderr,"ploidy= %d\n", ploid);



     // compute stats: average depth : default and actual after 1,1 and q25 filtering

       avDP4=ceil(sumDP4/count_afterF);//actual average cov  after Filt
       stdDP4=ceil(sqrt(sumDP4_2/count_afterF -avDP4*avDP4));//
       perc_left=(100.0*count_afterF/count_beforeF);

     fprintf(stderr,"average Depth after filtering= %d\n", avDP4);
     fprintf(stderr,"std Depth appr= %d\n", stdDP4);
     fprintf(stderr,"percent data left after q25 & min_depth filtering thrR=1 thrA=1 = %.2f\n", perc_left);

    // adjust actual depth to Bad (deep covered) regions and recommend thr RA

        avDP41=avDP4;
        ratio=(sqrt(sumDP4_2/count_afterF-avDP4*avDP4))/(sumDP4/count_afterF);//stdDP4/avDP4;
        fprintf(stderr,"cv actual depth = %.2f\n", ratio);

          if ( ratio > 1.5) //bad regions
          {
			 avDP41=ceil(0.8*avDP4);
             fprintf(stderr,"moderately deep covered regions here\n");
             fprintf(stderr,"adjusted(avDP4) = %d\n", avDP41);
          }

          if (ratio > 3) //very bad regions
		  {
			  avDP41=ceil(0.7*avDP4);
		      fprintf(stderr,"extremely deep covered regions here\n");
		      fprintf(stderr,"adjusted(avDP4) = %d\n", avDP41);
           }

           // recommendation for thrRR AA, depending on adjasted Depth
           if ( avDP41 > 0 && avDP41 < 4 )
		     	     {
					 thRR = 100;thAA=100;//stupidly high?
		     	     printf("better not to bother to look for mixture: too low actual coverage=%d\n",avDP4);
		     	     }
		   if ( avDP41 >= 4 && avDP41 < 10  )
		     		 {
				     thRR = 1;thAA=1;
		     	     }
		   if ( avDP41 >= 10 && avDP41 < 70  )
		     		 {
					 thRR = 2;thAA=2;
		     	     }
		   if ( avDP41 >= 70 && avDP41 < 90  )
		     		 {
				     thRR = 3;thAA=2;
		     	     }
		   if ( avDP41 >= 90 )
		     		 {
					  thRR = 3;thAA=3;
		  	         }

	}// if large enough vcfq file

      //====================== output1   write thrs, avD, stdD and ploidy
        fprintf(thrFile,"%d %d %d %d %d\n",thRR,thAA,avDP41,stdDP4,ploid);

    fclose(extract_vcfFile);
    fclose(thrFile);
    fclose(distribFile);
    fclose(filtered_vcfFile);

       // checking write/read
        if( firstSixAreOK  == 0)// || canWriteFilteredVCF == 0)
        {
		        fprintf(stderr,"Error during execution. Details: \n");
		        fprintf(stderr,"\tfirstSixAreOK %d\n",  firstSixAreOK);
		        fprintf(stderr,"\tcanWriteDistrib %d\n",      canWriteDistrib);
		        fprintf(stderr,"\tcanWriteFilteredVCF %d\n",   canWriteFilteredVCF);
		        fprintf(stderr,"Execution aborted\n");
		        exit(EXIT_FAILURE);
        }

    fprintf(stderr," done get_thr_ploidy.\n");
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

