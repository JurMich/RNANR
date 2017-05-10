/****************************************************************************************/
/*     Reads an RNA Sequence and outputs all locally optimal secondary structures       */ 
/****************************************************************************************/
 
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <mpfr.h>
#include "rna.h"
#include "output_file.h"
#include "base_pairs.h"
#include "flat_structures.h"
#include "locally_optimal_structures.h"
#include "energies.h"
#include "counting_double.h"
#include "counting_mpfr.h"

int FILEBP; 
int FILEOUT;
int VERBOSE;
int COUNTING;
int PARTITION;
int FLATS; 
int SAMPLING;
int NONREDUN;
int GRAMMAR;
int TIMEINFO;
int USEMPFR;
float ZSAMPLING; 
char * nameFileOut; 
char * nameFileInFasta;
char * nameFileGrammar;
char * nameFileBasePair;

/* Print usage */
void usage(){
  printf("\nRNANR -i sequence.fa [-o outputfile] [-a int] [-c] [-u] [-f] [-g outputfile] [-d file.bra] [-e file.bra] [-k int] [-l int] [-m int] [-n int] [-p int]  [-q int] [-t int] [-x] [-h] [-v]\n\n");
  printf("   -a <int> \n"); 
  printf("        minimum helix length. Default is 3.\n"); 
  printf("   -c count\n"); 
  printf("   -u calculates the number of flat structures;\n");
  printf("   -f calculates Boltzmann's partition function.\n");
  printf("   -w uses MPFR to augment precision of compututation, allowing deeper samplig (slower)."); 
  printf("   -g outputfile \n");
  printf("       prints all flat structures in adapted grammar.\n");
  printf("   -s <int> \n"); 
  printf("       non-redundant stochastic backtrack, giving entered number of sampled structures.\n");
  printf("   -z <int> \n"); 
  printf("       non-redundant stochastic backtrack, giving a number of structures that constitute entered percentage of Boltzmann's partition function.\n");
  printf("   -r available only with option, enables redundancy for stochastic backtrack. Does not work for -z, as to prevent possible very long execution times\n"); 
  printf("   -d file.bra\n");
  printf("       set of base pairs that will be used for the foldings.\n");
  printf("   -e file.bra\n");
  printf("       set of base pairs that will be used for the foldings in addition to the A-U, G-C, G-U pairings\n");  
  printf("   -i sequence.fa \n"); 
  printf("        input sequence in Fasta format \n");
  printf("   -k <int> \n");
  printf("       minimum hairpin loop length. Default is 3.\n");
  printf("   -l <int> \n"); 
  printf("        maximum hairpin loop length. Default is 12.\n"); 
  printf("   -m <int> \n");
  printf("        maximum loop length size. Default is 7.\n"); 
  printf("   -n <int> \n"); 
  printf("        maximum span of a base pair. Default is the length of the input RNA sequence.\n"); 
  printf("   -o outputfile \n"); 
  printf("        name of the output file (text file). This file contains the list of all locally \n"); 
  printf("        optimal secondary structures found in bracket-dot format with the number of base \n");
  printf("        pairs present in the structure. Default is stdout. (does not work with -c)\n");
  printf("   -p <int> \n"); 
  printf("        minimum percentage of paired bases (does not work with -c).\n");       
  printf("        Default is 0.\n");
  printf("   -q <int> \n");
  printf("        maximum degree of a multiloop. Default is: no limit.\n");
  printf("   -t <int>\n");
  printf("        temperature (in °C), default value is 37°C.\n");
  printf("   -b <int>\n");
  printf("        temperature scaling factor that does not influence loop energies. Default value is 1.0.\n");
  printf("   -x CPU execution time. Shows a CPU execution time for different parts of code. \n");
  printf("   -h print this help and exit\n"); 
  printf("   -v verbose. In this mode, the set of base pairs and the set of flat structures is output.\n");
  
  printf("\n\nSee the README file for more information about these options. \n\n"); 
  
}

/* Read parameters */
void parse_params(int argc, char **argv){
  
  char char_read;
  
  while((char_read = getopt(argc, argv, "a:b:cufwd:g:e:i:k:l:m:n:b:vo:p:q:s:t:z:xr")) != EOF){
    switch(char_read){ 
    case 'm' : /* bulge length */
      MAX_BULGE_SIZE= atoi(optarg);
      if (MAX_BULGE_SIZE<0){
	printf("\nBad parameter value (-m): This value should be positive.\n\n"); 
	usage(); 
	exit(EXIT_SUCCESS);
      }
      break; 
    case 'n': /* max helix span */
      MAX_HELIX_SCOPE=atoi(optarg);
      if (MAX_HELIX_SCOPE<2){
	printf("\nBad parameter value (-n): This value should larger than 1.\n\n"); 
	usage(); 
	exit(EXIT_SUCCESS);
      }
      break;
    case 'c' : /* counting */
      COUNTING=1; 
      break; 
    case 'u': /* get numbers of flat structures */
      FLATS=1;
      break;
    case 'f' : /* calculating Boltzmann's partition function */
      PARTITION=1;
      break;
    case 'w' : /* uses MPFR instead of standard doubles */
	  USEMPFR=1;
	  break;
    case 'g' : /* printing in adapted grammar */
      nameFileGrammar = (char*) malloc (strlen(optarg)+1 * sizeof(char));
      strcpy (nameFileGrammar, optarg);
      GRAMMAR=1;
      break;    
    case 's' : /* stochastic backtrack */
      SAMPLING= atoi(optarg);
      if(SAMPLING<0){
	  	printf("Bad parameter value (-s): The value entered is too high or negative. \n\n"); 
		usage(); 
	    exit(EXIT_SUCCESS);  
	  }
      break;
    case 'r' : /* redundant stochastic backtrack */
      NONREDUN=1;
      break; 
    case 'i' : /* input sequence in Fasta format*/
      nameFileInFasta=(char*) malloc ((strlen(optarg)+1) * sizeof(char));
      strcpy (nameFileInFasta, optarg);
      break; 
    case 'o' : /* output file */
      nameFileOut = (char*) malloc (strlen(optarg)+1 * sizeof(char));
      strcpy (nameFileOut, optarg);
      FILEOUT=1;
      break;
    case 'd' : /* base pair file */
      nameFileBasePair = (char*) malloc (strlen(optarg)+1 * sizeof(char));
      strcpy (nameFileBasePair, optarg);
      FILEBP=1;
      break;
    case 'e' : /* base pair file */
      nameFileBasePair = (char*) malloc (strlen(optarg)+1 * sizeof(char));
      strcpy (nameFileBasePair, optarg);
      FILEBP=2;
      break;
    case 'v': /*verbose */
      VERBOSE=1; 
      break;
    case 'k' :/* minimum size of the terminal loop */
      MIN_LOOP_SIZE= atoi(optarg); 
      if (MIN_LOOP_SIZE<=0){
	printf("Bad parameter value (-k): This value should be positive.\n\n"); 
	usage(); 
	exit(EXIT_SUCCESS);
      }
      break; 
    case 'a': /* minimum length of a helix */
      MIN_HELIX_LENGTH=atoi(optarg);
      if (MIN_HELIX_LENGTH<=0){
	printf("Bad parameter value (-a): This value should be positive.\n\n"); 
	usage(); 
	exit(EXIT_SUCCESS);
      }
      break; 
    case 'l': /* maximum size of the terminal loop */
      MAX_LOOP_SIZE= atoi(optarg); 
      if (MAX_LOOP_SIZE<0){
	printf("Bad parameter value (-l): This value should be positive or 0.\n\n"); 
	usage(); 
	exit(EXIT_SUCCESS);
      }
      break; 
    case 'p' : /* minimal number of base pairs */
      MIN_PERCENT_PAIRING=atoi(optarg); 
      if ((MIN_PERCENT_PAIRING<0) ||(MIN_PERCENT_PAIRING>100)) {
	printf("Bad parameter value (-p): This value should range between 0 and 100.\n\n"); 
	usage(); 
	exit(EXIT_SUCCESS);
      }
      break;
       case 'q': /* maximum degree of a multiloop */
      MAX_DEGREE=atoi(optarg);
      if (MAX_DEGREE<=0){
	printf("Bad parameter value (-q): This value should be positive.\n\n"); 
	usage(); 
	exit(EXIT_SUCCESS);
      }
      break;
    case 't':
      TEMP =  atof(optarg);
      if (TEMP<-273){ 
	printf("Bad parameter value (-t): This value should be higher than -273°C.\n\n"); 
	usage(); 
	exit(EXIT_SUCCESS);	 
      }
      break;
    case 'b':
      TEMPSCALE =  atof(optarg);
      if (TEMPSCALE<0){ 
	printf("Bad parameter value (-b): This value should be positive.\n\n"); 
	usage(); 
	exit(EXIT_SUCCESS);	 
      }
      break; 
    case 'x':
	  TIMEINFO=1;
	  break; 
	case 'z':
	  ZSAMPLING=atof(optarg);
	  if ((ZSAMPLING<=0.0) || (ZSAMPLING>1.0)){ 
		usage("Bad parameter value (-z): This value should be between 0 and 1."); 
	    exit(EXIT_SUCCESS);	
	  }
	  break;  
    case '?':
      if(optopt == 'o')
	printf("Option -%c requires an argument \n\n", optopt);
      usage();
      exit(EXIT_SUCCESS);
    } /* end switch */
  } /* end while */

 if (nameFileInFasta[0]==' ') {
    printf("\n! ! MISSING INPUT FILE (-i) ! !  \n \n"); 
    usage(); 
    exit(EXIT_SUCCESS);
  }

}

void check_and_update_params(int length){  
 
  if (MAX_HELIX_SCOPE==0)
    MAX_HELIX_SCOPE=length; 
  if (MAX_DEGREE==0)
    MAX_DEGREE=length;
  
  printf("\n"); 
  printf("Minimum length of a helix (-a): %d\n", MIN_HELIX_LENGTH);
  printf("Minimum size of a terminal hairpin loop (-k): %d\n", MIN_LOOP_SIZE); 
  printf("Maximum size of a terminal hairpin loop (-l): %d\n", MAX_LOOP_SIZE); 
  printf("Maximum size of an internal loop (-m): %d\n", MAX_BULGE_SIZE);
  printf("Maximum span of a base pair (-n): %d\n", MAX_HELIX_SCOPE);
  printf("Maximum degree od a multiloop (-q): %d\n", MAX_DEGREE); 

  if (MAX_HELIX_SCOPE < MIN_LOOP_SIZE + 2*MIN_HELIX_LENGTH){
    printf("\n! ! Bad parameter value (-n): This value is too small.\n\n"); 
    usage(); 
    exit(EXIT_SUCCESS); 
  }
 if (MAX_LOOP_SIZE < MIN_LOOP_SIZE){
    printf("\n! ! Bad parameter value (-l): This value should be higher than the minimum size of the hairpin loop (-k).\n\n"); 
    usage(); 
    exit(EXIT_SUCCESS);
  }
  if (COUNTING){ 
    printf("\n Counting locally optimal structures only.\n");
  }
  else if(PARTITION){
	printf("\n Counting Boltzmann's partition function only.\n");
  }		
  else {
    printf("Minimum percentage of positions in a base pairing (-p): %d\n", MIN_PERCENT_PAIRING); 
    if (FILEOUT)
      printf("Output file: %s\n", nameFileOut);
  }
  printf("\n"); 
}

/*builds up a general*/
void print_time_general(plain_sequence * rna_seq, double BP_t, double flat_t){
  printf("\nExecution time details for sequence of length %d nt: \n", rna_seq->size);
  printf(" - base pair determination: %f \n", BP_t);
  printf(" - flat structure computation: %f \n", flat_t);	
}


int main(int argc, char **argv){
	
  char RNAname[200];  
  FILE * outfile; 
   
  plain_sequence  * rna_seq;
  
  //srand(time(NULL)); // moved to counting_double.c since that is only one using this kind of generator
  
  opterr = 0;
  MIN_LOOP_SIZE=3; 
  MAX_LOOP_SIZE=12;
  MAX_BULGE_SIZE=7; 
  MIN_PERCENT_PAIRING=0;
  MIN_HELIX_LENGTH=3;   
  MAX_DEGREE=0; 
  FILEOUT=0;
  VERBOSE=0;
  COUNTING=0;
  PARTITION=0;
  FLATS=0; 
  SAMPLING=0;
  GRAMMAR=0;
  USEMPFR=0;
  TEMP=37.;
  TEMPSCALE=1.; 
  nameFileOut=""; 
  nameFileInFasta=" ";
  nameFileGrammar="";
  nameFileBasePair=""; 

  if (argc<2) {
    usage(); 
    exit(EXIT_SUCCESS);
  }
  
  /* for testing purpose only */
  struct timespec start_time, end_time;
  double BP_t;
  double flat_t;
  double DP_t;
  
  parse_params(argc, argv);
  RT = TEMPSCALE*0.0019872370936902486 * (273.15 + TEMP) * 100;
  rna_seq= (plain_sequence *) get_plain_sequence(nameFileInFasta, RNAname);	 
  check_and_update_params(rna_seq->size);
  if(!USEMPFR){
	  if(NONREDUN && ZSAMPLING){
		printf("Warning : Non-redundancy option '-r' cannot be used with an option '-z' : '-r' ignored.\n");
		NONREDUN = 0;  
	  }
	  
	  /* generate all base pairs */
	  if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &start_time);
	  compute_all_base_pairs(rna_seq, FILEBP, nameFileBasePair);
	  if (TIMEINFO){
		clock_gettime(CLOCK_MONOTONIC, &end_time);
		BP_t = (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/1.0e9 ; 
	  }
	  if (VERBOSE) 
		display_base_pairs(rna_seq);
	  
	  /* generate all flat structures */
	  if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &start_time);
	  build_all_flat_structures(rna_seq); 
	  if (TIMEINFO){
		clock_gettime(CLOCK_MONOTONIC, &end_time);
		flat_t = (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/1.0e9 ;
	  }
	  if (VERBOSE)
		display_all_flat_structures(rna_seq); 
	  
	  if (COUNTING){
		/* count all loc opt structures */ 
		if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &start_time);  
		count_all_locally_optimal_structures(rna_seq);
		if(TIMEINFO){
		  clock_gettime(CLOCK_MONOTONIC, &end_time);
		  DP_t = (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/1.0e9 ;
		  print_time_general(rna_seq, BP_t, flat_t);
		  printf(" - dynamic programming - counting: %f \n", DP_t);
		} 
	  }
	  else if(PARTITION){
		if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &start_time);
		double partition_fun;
		get_partition_function_double(&partition_fun, rna_seq);
		printf("Boltzmann's partition function: %.9f \n", partition_fun);
		if(TIMEINFO){
		  clock_gettime(CLOCK_MONOTONIC, &end_time);
		  DP_t = (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/1.0e9 ;
		  print_time_general(rna_seq, BP_t, flat_t);
		  printf(" - dynamic programming matrix calculation: %f \n", DP_t);
		} 
	  }
	  else if(FLATS){
		printf("Total number of non-empty flat structures: %ld \n", count_all_flat_structures(rna_seq));  
	  }
	  else if(GRAMMAR){
		printf("All flat structures are exported to file in adapted grammar. Terminal output also available below: \n\n");
		print_all_flat_structures_pile(rna_seq, nameFileGrammar); 
	  }	  
	  else{ 
		if (SAMPLING>0){
		  /* stochastic backtrack*/ 
		  int i;
		  folding* folded_rna;
		  int struc_count = 0;
		  if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &start_time);  
		  folded_rna = stochastic_backtrack_locally_optimal_structures_double(SAMPLING, rna_seq, NONREDUN, TIMEINFO, 0, &struc_count, &DP_t);
		  if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &end_time);
		  display_plain_sequence(rna_seq); 
		  for (i=1;i<=SAMPLING;i++){
			print_structure(folded_rna->structures[i], folded_rna->part_fcis[i], folded_rna->energy_ref[i], rna_seq);
		  }
		  if(TIMEINFO){
			print_time_general(rna_seq, BP_t, flat_t);
			printf(" - sampling time (%d samples): %f \n", SAMPLING, DP_t);
		  } 
		}
		else if(ZSAMPLING>0){
		  /* stochastic backtrack*/ 
		  int i;
		  folding* folded_rna;
		  int struc_count = 0;
		  if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &start_time);  
		  folded_rna = stochastic_backtrack_locally_optimal_structures_double(SAMPLING, rna_seq, NONREDUN, TIMEINFO, ZSAMPLING, &struc_count, &DP_t);
		  if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &end_time);
		  display_plain_sequence(rna_seq);
		  for (i=1;i<=struc_count;i++){
			print_structure(folded_rna->structures[i], folded_rna->part_fcis[i], folded_rna->energy_ref[i], rna_seq);
		  }
		  if(TIMEINFO){
			print_time_general(rna_seq, BP_t, flat_t);
			printf(" - sampling time (%d samples): %f \n", SAMPLING, DP_t);
		  } 
		}
		else
		{
		  /* print all loc opt structures */
		  printf("-- Locally optimal secondary structures\n\n"); 
		  if (!FILEOUT) 
			display_plain_sequence(rna_seq);
		  outfile=open_outputfile(FILEOUT, nameFileOut); 
		  if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &start_time);  
		  generate_all_locally_optimal_structures(rna_seq,outfile);
		  if(TIMEINFO){
			clock_gettime(CLOCK_MONOTONIC, &end_time);
			print_time_general(rna_seq, BP_t, flat_t);
			printf(" - determining of structures: %f \n", DP_t);
		  } 
		  close_outputfile(FILEOUT, outfile); 
		  if (FILEOUT) 
			printf("See %s", nameFileOut); 
		}
	  }
	  
	  free_base_pairs(rna_seq); 
	  free_plain_sequence(rna_seq);
	  printf("\nBye bye\n");
  }else{
	  if(NONREDUN && ZSAMPLING){
		printf("Warning : Non-redundancy option '-r' cannot be used with an option '-z' : '-r' ignored.\n");
		NONREDUN = 0;  
	  }
	  
	  /* generate all base pairs */
	  if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &start_time);
	  compute_all_base_pairs(rna_seq, FILEBP, nameFileBasePair);
	  if (TIMEINFO){
		clock_gettime(CLOCK_MONOTONIC, &end_time);
		BP_t = (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/1.0e9 ; 
	  }
	  if (VERBOSE) 
		display_base_pairs(rna_seq);
	  
	  /* generate all flat structures */
	  if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &start_time);
	  build_all_flat_structures(rna_seq); 
	  if (TIMEINFO){
		clock_gettime(CLOCK_MONOTONIC, &end_time);
		flat_t = (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/1.0e9 ;
	  }
	  if (VERBOSE)
		display_all_flat_structures(rna_seq); 
	  
	  if (COUNTING){
		/* count all loc opt structures */ 
		if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &start_time);  
		count_all_locally_optimal_structures(rna_seq);
		if(TIMEINFO){
		  clock_gettime(CLOCK_MONOTONIC, &end_time);
		  DP_t = (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/1.0e9 ;
		  print_time_general(rna_seq, BP_t, flat_t);
		  printf(" - dynamic programming - counting: %f \n", DP_t);
		} 
	  }
	  else if(PARTITION){
		if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &start_time);
		mpfr_t partition_fun;
		mpfr_init2(partition_fun, 100);
		get_partition_function_mpfr(&partition_fun, rna_seq);
		printf("Boltzmann's partition function: %.9f \n", mpfr_get_d(partition_fun, MPFR_RNDN));
		if(TIMEINFO){
		  clock_gettime(CLOCK_MONOTONIC, &end_time);
		  DP_t = (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/1.0e9 ;
		  print_time_general(rna_seq, BP_t, flat_t);
		  printf(" - dynamic programming matrix calculation: %f \n", DP_t);
		} 
	  }
	  else if(FLATS){
		printf("Total number of non-empty flat structures: %ld \n", count_all_flat_structures(rna_seq));  
	  }
	  else if(GRAMMAR){
		printf("All flat structures are exported to file in adapted grammar. Terminal output also available below: \n\n");
		print_all_flat_structures_pile(rna_seq, nameFileGrammar); 
	  }	  
	  else{ 
		if (SAMPLING>0){
		  /* stochastic backtrack*/ 
		  int i;
		  folding* folded_rna;
		  int struc_count = 0;
		  if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &start_time);  
		  folded_rna = stochastic_backtrack_locally_optimal_structures_mpfr(SAMPLING, rna_seq, NONREDUN, TIMEINFO, 0, &struc_count, &DP_t);
		  if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &end_time);
		  display_plain_sequence(rna_seq); 
		  for (i=1;i<=SAMPLING;i++){
			print_structure(folded_rna->structures[i], folded_rna->part_fcis[i], folded_rna->energy_ref[i], rna_seq);
		  }
		  if(TIMEINFO){
			print_time_general(rna_seq, BP_t, flat_t);
			printf(" - sampling time (%d samples): %f \n", SAMPLING, DP_t);
		  } 
		}
		else if(ZSAMPLING>0){
		  /* stochastic backtrack*/ 
		  int i;
		  folding* folded_rna;
		  int struc_count = 0;
		  if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &start_time);  
		  folded_rna = stochastic_backtrack_locally_optimal_structures_mpfr(SAMPLING, rna_seq, NONREDUN, TIMEINFO, ZSAMPLING, &struc_count, &DP_t);
		  if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &end_time);
		  display_plain_sequence(rna_seq);
		  for (i=1;i<=struc_count;i++){
			print_structure(folded_rna->structures[i], folded_rna->part_fcis[i], folded_rna->energy_ref[i], rna_seq);
		  }
		  if(TIMEINFO){
			print_time_general(rna_seq, BP_t, flat_t);
			printf(" - sampling time (%d samples): %f \n", SAMPLING, DP_t);
		  } 
		}
		else
		{
		  /* print all loc opt structures */
		  printf("-- Locally optimal secondary structures\n\n"); 
		  if (!FILEOUT) 
			display_plain_sequence(rna_seq);
		  outfile=open_outputfile(FILEOUT, nameFileOut); 
		  if (TIMEINFO) clock_gettime(CLOCK_MONOTONIC, &start_time);  
		  generate_all_locally_optimal_structures(rna_seq,outfile);
		  if(TIMEINFO){
			clock_gettime(CLOCK_MONOTONIC, &end_time);
			print_time_general(rna_seq, BP_t, flat_t);
			printf(" - determining of structures: %f \n", DP_t);
		  } 
		  close_outputfile(FILEOUT, outfile); 
		  if (FILEOUT) 
			printf("See %s", nameFileOut); 
		}
	  }
	  
	  free_base_pairs(rna_seq); 
	  free_plain_sequence(rna_seq);
	  printf("\nBye bye\n");		
  }	  
  return 0;
}

 
