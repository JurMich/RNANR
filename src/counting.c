/**************************************************/
/*                   counting                     */
/**************************************************/

#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/loop_energies.h>
#include "tree.h"
#include "rna.h"
#include "base_pairs.h"
#include "flat_structures.h"
#include "sort.h"
#include "stack.h"
#include "counting.h"

TYPE ** number_of_locopt_structures;
TYPE ** partition_function_table; 

typedef struct {int flat_start; int flat_end;} flat_pair;
 
/* parameters related with Vienna package*/
vrna_fold_compound_t *E_fold_cp;


/********************************************/
/*    operators for the locopt structures   */
/********************************************/

TYPE sum(TYPE a, TYPE b){
  return a+b;
}

TYPE product(TYPE a, TYPE b){
  return a*b; 
}

TYPE empty_set(){ /* identity element for PLUS */
  return 0;
}

TYPE empty_structure(){ /* identity element for PRODUCT */
  return 1;
}

TYPE apply_energy_term(TYPE dG, TYPE val){ 
  /* modifying term for energy (Boltzmann factor)*/
  return exp(-dG/RT)*val;   /*Boltzmann's cte already contains temperature */
  return dG*val;
}

/**********    end operators    ***********/


/********************************************/
/*       structure and misc. utilities      */
/********************************************/
void print_structure(char * structure, TYPE energy, TYPE energy_ref, plain_sequence * rna){
  int i;
  for (i=1; i<=rna->size; i++)
    printf("%c", structure[i]);
  printf(" - Free Energy : %.3f \n", (-RT*log(energy))/100);
}

void add_base_pair(int i, int j, char * structure){
  structure[i] = '(';
  structure[j] = ')';
}

void add_helix(int i, int j, int thickness, char * structure){
  int k;
  for (k=0;k<thickness;k++){
    add_base_pair(i+k, j-k, structure);
  }
}

/* this converts base into a number compatible with ViennaRNA's
 * vrna_md_t->pair/ptype definition*/
int base2int(char base){
  switch(base){
	case 'a': return 1;
	case 'A': return 1;
	case 'c': return 2;
	case 'C': return 2;
	case 'g': return 3;
	case 'G': return 3;
	case 'u': return 4;
	case 'U': return 4;
	default : return 0;  
  }	
}

/* obtains a type as defined by ViennaRNA of a pair entered in parameters.*/
int get_type(char base_5, char base_3){ /* base_5 - from 5' end, base_3 - from 3' end */
  int b_5, b_3, type;
  b_5 = base2int(base_5);
  b_3 = base2int(base_3);
  type = E_fold_cp->params->model_details.pair[b_5][b_3];
  return type;
}

TYPE random_number(TYPE upper_bound){
  return (double)rand()*upper_bound/((double)RAND_MAX + 1.);
}

/* prints a character to terminal repeatedly*/
void repeat_print(char*to_print, int n_repeats){
  int i;
  for(i=0; i<n_repeats; i++){
	printf("%s ", to_print);
  }	
}

/* checks if flat pair is already in table */
int is_in_table(int flat_start, int flat_end, flat_pair * pile, int pile_size){
   int i;
   for(i = 0; i<pile_size; i++){
      if((pile[i].flat_start == flat_start)&&(pile[i].flat_end == flat_end)){
		 return 1;  
	  } 
   } 
   return 0;	
}

/************ Free energy utilities ************/

/* Energy of (i+1,j-1) stacking over (i,j), i+1<j-1 */
TYPE stacking_energy(int i, int j){
  TYPE stack_energy;
  stack_energy = vrna_E_stack(E_fold_cp, i, j);  /*ViennaRNA*/
  return stack_energy;
}

/* Energy of helix spanning from (i,j) to (i+length,j-length), i+length<j-length */
TYPE helix_energy(int i, int j, int length){
  TYPE dG;
  int k;
  
  dG = 0.;
  for(k=0;k<length-1;k++){
	dG += stacking_energy(i+k, j-k);
  }
  return dG;
}

/* Energy of multiloop with alpha helixes and beta unpaired bases 
 * - si - a table of a 5' nucleotide coordinates next to multiloop
 * - sj - a table of a 3' nucleotide coordinates next to multiloop  */
TYPE multiloop_energy(int alpha, int beta, int *ptypes, int *si, int *sj){
  	int i, closing_pty, unpaired_pty; /*to extract from Vienna package*/
  	TYPE mtloop_E;
  	closing_pty = E_fold_cp->params->MLclosing; /* ML closing penalty */
  	unpaired_pty = E_fold_cp->params->MLbase; /* unpaired base penalty */
  	mtloop_E = closing_pty + unpaired_pty*beta;
  	for(i = 0; i <= alpha; i++){
	   mtloop_E += E_MLstem(ptypes[i], E_fold_cp->sequence_encoding[si[i]],\
	    E_fold_cp->sequence_encoding[sj[i]], E_fold_cp->params);
	}
  	return mtloop_E;
}

/* Energy of internal/bulge loop defined by two base pairs (i,l) and (j,k), i<k<l<j */
TYPE internal_loop_energy(int i, int k, int l, int j, plain_sequence *rna){
  TYPE it_energy;
  int size1, size2, ptype1, ptype2, base_i, base_j, base_k, base_l;
  
  size1 = k-i-1;
  size2 = j-l-1;
  ptype1 = get_type(rna->label[i], rna->label[j]);
  ptype2 = get_type(rna->label[l], rna->label[k]);
  base_i = E_fold_cp->sequence_encoding[i+1];
  base_j = E_fold_cp->sequence_encoding[j-1];
  base_k = E_fold_cp->sequence_encoding[k-1];
  base_l = E_fold_cp->sequence_encoding[l+1];
  it_energy = E_IntLoop(size1, size2, ptype1, ptype2, base_i, base_j, base_k, base_l, E_fold_cp->params);
  return it_energy;
}

/* computes energy of a hairpin */
TYPE hairpin_energy(int i, int j, plain_sequence * rna){
  TYPE hairpin_energy;
  int hsize, ptype;
  char *hseq; /* sequence of unpaired bases in hairpin */
  hsize = j - i - 1;
  ptype = get_type(rna->label[i], rna->label[j]);
  hseq = (char*) malloc ((hsize+1)*sizeof(char));
  strncpy(hseq, rna->label+i+1, hsize);
  hseq[hsize] = '\0';
  hairpin_energy = E_Hairpin(hsize, ptype, E_fold_cp->sequence_encoding[i+1], E_fold_cp->sequence_encoding[j-1], hseq, E_fold_cp->params);
  //hairpin_energy = vrna_E_hp_loop(E_fold_cp, i, j);
  //printf("ptype: %d, energy: %f, hseq: %s, hsize: %d\n", ptype, hairpin_energy, hseq, hsize);
  free(hseq);
  return hairpin_energy;
}

/* computes energy of an exterior loop */
TYPE ext_loop_energy(int i, int j, plain_sequence *rna){
  TYPE ext_loop_energy;
  int ptype, seq_i, seq_j;
  ptype = get_type(rna->label[i], rna->label[j]);
  seq_i = (i > 1) ? E_fold_cp->sequence_encoding[i-1] : -1;
  seq_j = (j < rna->size) ? E_fold_cp->sequence_encoding[j+1] : -1;
  ext_loop_energy = E_ExtLoop(ptype, seq_i, seq_j, E_fold_cp->params); 
  return ext_loop_energy;	
}

int is_entirely_contained(int i, int j, plain_sequence * rna){
  return ((i==1) && (j==rna->size));
}

/***********   end utilities     ************/

void init_locopt_table(int n){
  int x, y; 
  number_of_locopt_structures=(TYPE **) malloc ((n+1)*sizeof(TYPE *));
  for (x=0; x<=n; x++){
    number_of_locopt_structures[x]=(TYPE *) malloc ((n+1)*sizeof(TYPE));
    for (y=0; y<=n; y++){
      number_of_locopt_structures[x][y]= empty_set();/* yann */
    }
  } 
} 

// changes pointers within table (modifies directly flat_table via pointers)
void rewire(int x, int y, linktab * linktable){
	if(linktable->length > 0){
		flat_table[x][y] = linktable->links[0].fcell;
		for(int	i = 0; i<(linktable->length-1); i++){
			linktable->links[i].fcell->next = linktable->links[i+1].fcell;
		}
		linktable->links[linktable->length-1].fcell->next = NULL;
	}
}


TYPE count_all_locally_optimal_structures(plain_sequence * rna){
  TYPE n;
  int x, y, i, j;
  flat_cell * current_flat_structure;
  int current_base_pair;
  int thickness;
  init_locopt_table(rna->size);
  for (x=rna->size; x>=1; x--){
    for (y=x+1; y<=rna->size; y++){
      if (flat_table[x][y]!=NULL){
	current_flat_structure=flat_table[x][y];
	do {
	  current_base_pair=current_flat_structure->current;
	  i=get_flat_structure_start(current_base_pair);
	  j=get_flat_structure_end(current_base_pair); 
      /* TODO: The stacking case -> fetch energy bonus from Vienna package */
      if ((i==x) && (j==y) && !is_entirely_contained(i,j,rna)){ /* exception to cases where i=1 OR j=rna->size */
	    /* thickness= 1 */
        n = number_of_locopt_structures[x+1][y-1];  /* yann  */
	  }
	  else{
        /* TODO: Energy contribution ( distinguish Multiloop + internal loops) */
        /* TODO: Auxiliary function to compute overall contribution of flat structure*/
        n=empty_structure(); /* yann */
	    while(current_base_pair != 0){
	      i=get_flat_structure_start(current_base_pair);
	      j=get_flat_structure_end(current_base_pair);
	      get_flat_structure_suffix(current_base_pair);
	      thickness=get_BP(i,j); 
	      if (thickness>MIN_HELIX_LENGTH) 
		thickness=MIN_HELIX_LENGTH;
	      n = product(n, number_of_locopt_structures[i+thickness][j-thickness]); /* yann */
	      current_base_pair=get_flat_structure_suffix(current_base_pair);
	    } /* end while */
	  }/* endif */
	  number_of_locopt_structures[x][y] = sum(n,number_of_locopt_structures[x][y]); /* yann */
	  current_flat_structure=current_flat_structure->next;
	}while (current_flat_structure != NULL);
      }/*endif else */
      else{
	number_of_locopt_structures[x][y]=empty_structure(); 
      }
    }/* end for y */
  }/* end for x */
  printf("\nNumber of structures: %ld\n", (long int)number_of_locopt_structures[1][rna->size]); 
  return number_of_locopt_structures[1][rna->size];
}

/* stacks every flat structure into a pile and counts them */
int count_all_flat_structures(plain_sequence * rna){
	int x, y, i, j;
	int n_flat = 0;
	flat_cell * current_flat_structure;
	int current_base_pair;
	int thickness;
	int pile_size;
	flat_pair *flat_pile = (flat_pair*) malloc(sizeof(flat_pair));
	flat_pair tmp_pair;
	tmp_pair.flat_start = 1;
	tmp_pair.flat_end = rna->size;
	flat_pile[0] = tmp_pair;
	pile_size = 1;
	for(int k=0; k<pile_size; k++){
		x = flat_pile[k].flat_start;
		y = flat_pile[k].flat_end;
		if (flat_table[x][y]!=NULL){
	current_flat_structure=flat_table[x][y];
	do {
	  n_flat++;		
	  current_base_pair=current_flat_structure->current;
	  i=get_flat_structure_start(current_base_pair);
	  j=get_flat_structure_end(current_base_pair); 
      if ((i==x) && (j==y) && !is_entirely_contained(i,j,rna)){ /* exception to cases where i=1 OR j=rna->size */
	    /* thickness= 1 */
	    if(!is_in_table(x+1, y-1, flat_pile, pile_size)){
			tmp_pair.flat_start = x+1;
			tmp_pair.flat_end = y-1;
			pile_size++;
			flat_pile = (flat_pair*) realloc(flat_pile, sizeof(flat_pair)*pile_size);
			flat_pile[pile_size-1] = tmp_pair;
		}
	  }
	  else{
        /* TODO: Energy contribution ( distinguish Multiloop + internal loops) */
        /* TODO: Auxiliary function to compute overall contribution of flat structure*/
	    while(current_base_pair != 0){			
	      i=get_flat_structure_start(current_base_pair);
	      j=get_flat_structure_end(current_base_pair);
	      thickness=get_BP(i,j); 
	      if (thickness>MIN_HELIX_LENGTH) 
		    thickness=MIN_HELIX_LENGTH;		  
	      if(!is_in_table(i+thickness, j-thickness, flat_pile, pile_size)){
			tmp_pair.flat_start = i+thickness;
			tmp_pair.flat_end = j-thickness;
			pile_size++;
			flat_pile = (flat_pair*) realloc(flat_pile, sizeof(flat_pair)*pile_size);
			flat_pile[pile_size-1] = tmp_pair;
	      }
	      current_base_pair=get_flat_structure_suffix(current_base_pair);
	    } /* end while */
	  }/* endif */
	  current_flat_structure=current_flat_structure->next;
	}while (current_flat_structure != NULL);
      }/*endif else */
	}
	return n_flat;
}



/* stacks every flat structure into a pile, and in parallel prints them one by one */
void print_all_flat_structures_pile(plain_sequence * rna, char * file_path){
	int x, y, i, j, i0, l;
	int len_unpaired;
	flat_cell * current_flat_structure;
	int current_base_pair;
	int thickness;
	int pile_size;
	FILE *flat_base;
    flat_base = fopen(file_path, "w");
	flat_pair *flat_pile = (flat_pair*) malloc(sizeof(flat_pair));
	flat_pair tmp_pair;
	tmp_pair.flat_start = 1;
	tmp_pair.flat_end = rna->size;
	flat_pile[0] = tmp_pair;
	pile_size = 1;
	for(int k=0; k<pile_size; k++){
		x = flat_pile[k].flat_start;
		y = flat_pile[k].flat_end;
		if (flat_table[x][y]!=NULL){
	current_flat_structure=flat_table[x][y];
	do {
	  printf("F%dF%d -> fXX <<< ", x, y);	
	  fprintf(flat_base, "Flat_strcture|%d|%d|::", x, y);	
	  current_base_pair=current_flat_structure->current;
	  i=get_flat_structure_start(current_base_pair);
	  j=get_flat_structure_end(current_base_pair); 
      if ((i==x) && (j==y) && !is_entirely_contained(i,j,rna)){ /* exception to cases where i=1 OR j=rna->size */
	    /* thickness= 1 */
	    if(!is_in_table(x+1, y-1, flat_pile, pile_size)){
			tmp_pair.flat_start = x+1;
			tmp_pair.flat_end = y-1;
			pile_size++;
			flat_pile = (flat_pair*) realloc(flat_pile, sizeof(flat_pair)*pile_size);
			flat_pile[pile_size-1] = tmp_pair;
		}
        printf("n F%dF%d n ", x+1, y-1);
        fprintf(flat_base, "(%d,%d);[%d,%d];", x, y, x+1, y-1);
	  }
	  else{
        /* TODO: Energy contribution ( distinguish Multiloop + internal loops) */
        /* TODO: Auxiliary function to compute overall contribution of flat structure*/
        i0 = x - 1;
	    while(current_base_pair != 0){			
	      i=get_flat_structure_start(current_base_pair);
	      j=get_flat_structure_end(current_base_pair);
	      len_unpaired = i - i0 - 1;
	      i0 = j;
	      thickness=get_BP(i,j); 
	      if (thickness>MIN_HELIX_LENGTH) 
		thickness=MIN_HELIX_LENGTH;
		  if(len_unpaired>0)
		    printf( "sx%d ", len_unpaired);
		    repeat_print("n", thickness);
		  for(l = 0; l<thickness; l++)
		    fprintf(flat_base, "(%d,%d);", i+l, j-l);			  
	      if(!is_in_table(i+thickness, j-thickness, flat_pile, pile_size)){
			tmp_pair.flat_start = i+thickness;
			tmp_pair.flat_end = j-thickness;
			pile_size++;
			flat_pile = (flat_pair*) realloc(flat_pile, sizeof(flat_pair)*pile_size);
			flat_pile[pile_size-1] = tmp_pair;
	      }
		  if((j-i-2*thickness+1)>1){
		    printf("F%dF%d ", i+thickness, j-thickness);
		    fprintf(flat_base,"[%d,%d];", i+l, j-l);
		  }
		  else
		    printf("sx1 ");  
		    repeat_print("n", thickness);
	        current_base_pair=get_flat_structure_suffix(current_base_pair);
	    } /* end while */
	    len_unpaired = y - i0;
	    if(len_unpaired>0)
	      printf("sx%d ", len_unpaired);
	  }/* endif */
	  printf("\n");
	  fprintf(flat_base, "\n");
	  current_flat_structure=current_flat_structure->next;
	}while (current_flat_structure != NULL);
      }/*endif else */
      else{
	    printf("F%dF%d -> fC <<< sx%d \n",x,y,y-x+1);	 
	    fprintf(flat_base, "Flat_strcture|%d|%d|::();\n", x, y); 
      }
	}	
	fclose(flat_base);
}


/* Creates a table of Boltzmann partition*/
void init_partition_table(int n){
  int x, y; 
  partition_function_table=(TYPE **) malloc ((n+1)*sizeof(TYPE *));
  for (x=0; x<=n; x++){
    partition_function_table[x]=(TYPE *) malloc ((n+1)*sizeof(TYPE));
    for (y=0; y<=n; y++){
      partition_function_table[x][y]= empty_set();/* yann */
    }
  } 
}


TYPE get_partition_function(plain_sequence * rna){
  TYPE partition_fci; /* terms of partition function */
  int x, y, i, j, i0; 
  flat_cell * current_flat_structure;
  int current_base_pair;
  int thickness; 
  int alpha; /* number of helixes in loop */
  int beta; /* number of unpaired bases */
  char *rna_seq = (char*) malloc(sizeof(char)*(rna->size+1));
  strncpy(rna_seq, rna->label+1, rna->size);
  rna_seq[rna->size] = '\0';
  
  init_partition_table(rna->size);
  vrna_md_t md;
  set_model_details(&md);
  md.temperature = TEMP;  
  E_fold_cp = vrna_fold_compound(rna_seq, &md, VRNA_OPTION_DEFAULT);
  
  for (x=rna->size; x>=1; x--){
    for (y=x+1; y<=rna->size; y++){
	   linktab * linktable = start_linktab();
      if (flat_table[x][y]!=NULL){
	current_flat_structure=flat_table[x][y];
	do {
	  current_base_pair=current_flat_structure->current;
	  i=get_flat_structure_start(current_base_pair);
	  j=get_flat_structure_end(current_base_pair); 
      if ((i==x) && (j==y) && !is_entirely_contained(i,j,rna)){ /* exception to cases where i=1 OR j=rna->size */
	    /* thickness= 1 */
        partition_fci = apply_energy_term(stacking_energy(x-1,y+1), partition_function_table[x+1][y-1]);  /* yann  */ 
        add_link_element(linktable, partition_fci, current_flat_structure);
	  }
	  else{
        partition_fci = empty_structure();
        alpha = 0;
        beta = 0;
        i0 = x - 1;
		int *ptypes =  (int *) malloc (sizeof(int)); /* table of pair types for multibranch loops */
		int *si = (int *) malloc (sizeof(int));
		int *sj = (int *) malloc (sizeof(int));
		ptypes[0] =  get_type(rna -> label[y+1], rna -> label[x-1]); /* closing pair type */ 
		si[0] = y;
		sj[0] = x;     
	    while(current_base_pair != 0){
		  alpha ++;
	      i=get_flat_structure_start(current_base_pair);
	      j=get_flat_structure_end(current_base_pair);
	      /* external loop contribution */
		  if(is_entirely_contained(x,y,rna)){
	     	 partition_fci = apply_energy_term(ext_loop_energy(i, j, rna), partition_fci);
		  }
	      get_flat_structure_suffix(current_base_pair);
	      ptypes =  (int *) realloc (ptypes, (alpha+1)*sizeof(int));
	      si = (int *) realloc (si, (alpha+1)*sizeof(int));
	      sj = (int *) realloc (sj, (alpha+1)*sizeof(int));	      
	      ptypes[alpha] = get_type(rna -> label[i], rna -> label[j]); /* enclosed pair type(s) */
	      si[alpha] = i-1;
	      sj[alpha] = j+1;
		  beta += (i - i0 - 1);  
	      i0 = j;
	      thickness=get_BP(i,j);
	      if (thickness>MIN_HELIX_LENGTH) 
		thickness=MIN_HELIX_LENGTH;   
	      partition_fci = product(partition_fci, apply_energy_term(helix_energy(i,j,thickness),\
	       partition_function_table[i+thickness][j-thickness])); /* yann */
	      current_base_pair = get_flat_structure_suffix(current_base_pair);
	    } /* end while */
	    if((x != 1) || (y != rna->size)){  
	      if(alpha == 1){ /* internal loop or bulge since stack was taken care of before*/
            partition_fci = apply_energy_term(internal_loop_energy(x-1, i, j, y+1, rna), partition_fci); 
		  }else{ /*multiloop*/
			beta += (y - i0);
		    partition_fci = apply_energy_term(multiloop_energy(alpha, beta, ptypes, si, sj), partition_fci);
		  }	
		}
		free(ptypes);
		free(si);
		free(sj);
	  }/* endif */
	  /* next line is sum because we are adding energies for different secondary structures*/ 
	  partition_function_table[x][y] = sum(partition_fci, partition_function_table[x][y]); /* yann */
	  add_link_element(linktable, partition_fci, current_flat_structure);
	  current_flat_structure=current_flat_structure->next;
	}while (current_flat_structure != NULL);
      }/*endif else */
      else if ((x>1) && (y<rna->size)){ /* exception to cases where i=1 OR j=rna->size */
	    partition_function_table[x][y] = apply_energy_term(hairpin_energy(x-1,y+1, rna), empty_structure());
	  }
	  linktable = sort(linktable);
	  rewire(x, y, linktable);
	  free_linktab(linktable);
    }/* end for y */
  }/* end for x */
  
  free(rna_seq); 
  return partition_function_table[1][rna->size];
}

/* probabilistic backtracking of random sampling (redundant and non-redundant) */

/* calculates energy of a structure (for comparison) */
TYPE get_reference_energy(vrna_fold_compound_t *E_fold_cp, char * structure, plain_sequence * rna){
	char *fold = (char*) malloc(sizeof(char)*(rna->size+1));
    strncpy(fold, structure+1, rna->size);
    fold[rna->size] = '\0';
	float energy = vrna_eval_structure (E_fold_cp, fold);
	free(fold);
	return energy;
}


void stochastic_backtrack_locally_optimal_structure_rec(int x, int y, plain_sequence * rna, char * structure, 
 double * energy, tr_node ** actual_node, TYPE * part_fci);

/* DP_t and struc_count are references to execution time of DP and number of structures respectively */
folding* stochastic_backtrack_locally_optimal_structures(int number_of_structures, plain_sequence * rna, 
 int is_non_redun, int use_timer, float zqlimit, int * struc_count, double * DP_t){
  clock_t begin_time = clock();
  clock_t current_time;
  double time_spent;
  /* Precompute/cache #locOpts */
  int i,j, continuing_on;
  TYPE max_struct, total_Boltzmann;
  /* timer-related variables */
  struct timespec start_time, end_time;
  TYPE cumulative_en;
  /* struct holding RNA structure and energy*/
  folding *folding_energy = (folding*)malloc (sizeof(folding));
  max_struct = (TYPE)count_all_locally_optimal_structures(rna);
  if((max_struct < number_of_structures) && !(is_non_redun)){
	 printf("Error: maximum number of structures is %ld. Use number lower or equal to this with option -s.\n\
	 Alternatively, you can also use option -z for non-redundant sampling. \n\n", (long int)max_struct); 
	 exit(EXIT_SUCCESS);
  } 
  
  /* creates tree to ensure non-redundant sampling */
  tr_node * non_red_tree = create_root();
  tr_node * actual_node = non_red_tree; // pointer to actual node
    
  total_Boltzmann = get_partition_function(rna);
  folding_energy->structures = (char**) malloc (sizeof(char*));
  folding_energy->energies = (TYPE*) malloc (sizeof(TYPE));
  folding_energy->energy_ref = (TYPE*) malloc (sizeof(TYPE));
  if(use_timer) clock_gettime(CLOCK_MONOTONIC, &start_time);
  /* this is basically fancy 'for' loop, since it allows two different types of samplings without much redundancy */
  i = 1;
  continuing_on = 1;
  while(continuing_on){
	//printf("\nGenerating new structure: %d\n", i);
	folding_energy->structures = (char**) realloc (folding_energy->structures, (i+1)*sizeof(char*));
    folding_energy->energies = (TYPE*) realloc (folding_energy->energies, (i+1)*sizeof(TYPE));
    folding_energy->energy_ref = (TYPE*) realloc (folding_energy->energy_ref, (i+1)*sizeof(TYPE));
	folding_energy->energies[i] = empty_structure();
	
    char* structure;
    TYPE* energy  = (TYPE*) malloc (sizeof(TYPE));
    TYPE* part_fci  = (TYPE*) malloc (sizeof(TYPE));
    structure = (char*) malloc ((rna->size+1)*sizeof(char));
    for (j=1;j<=rna->size;j++){
      structure[j] = '.';
    }
    *energy = 1.;
    *part_fci = partition_function_table[1][rna->size]; // = Z
    stochastic_backtrack_locally_optimal_structure_rec(1, rna->size, rna, structure, energy, &actual_node, part_fci);   
    folding_energy->structures[i] = structure;
    folding_energy->energies[i] = *energy;
    cumulative_en += (*energy);
    //printf("Total sum of chosen partition part: %0.3f \n", cumulative_en);
    folding_energy->energy_ref[i] = get_reference_energy(E_fold_cp, structure, rna);
    if(!is_non_redun) actual_node = traceback_to_root(actual_node, *energy);
    current_time = clock();
    time_spent = (double)(current_time-begin_time) / CLOCKS_PER_SEC;
    printf("clock: - %f - ", time_spent);
    print_structure(structure, *energy, get_reference_energy(E_fold_cp, structure, rna), rna);
    free(energy);
    i+=1;
    if(zqlimit != 0){
		if(zqlimit <= (float)cumulative_en/(float)total_Boltzmann){
			continuing_on = 0;
			*struc_count = i-1;
		}
	}else{
		if(i > number_of_structures){
			continuing_on = 0;
		}	
	}
  } 
  if(use_timer){
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    *DP_t = (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/1.0e9 ;
  } 
  return folding_energy;
}

void stochastic_backtrack_flat_structure_rec(flat_cell * current_flat_structure, plain_sequence * rna, 
 char * structure, double * energy, tr_node ** actual_node, TYPE * part_fci){
  int current_base_pair;
  int i, j; 
  int thickness; 
						
  current_base_pair = current_flat_structure->current;
  while(current_base_pair != 0){
    i=get_flat_structure_start(current_base_pair);
    j=get_flat_structure_end(current_base_pair);
	thickness=get_BP(i,j); 
    if (thickness>MIN_HELIX_LENGTH) 
      thickness=MIN_HELIX_LENGTH;
    add_helix(i, j, thickness, structure);
    stochastic_backtrack_locally_optimal_structure_rec(i+thickness, j-thickness, rna, structure, energy, actual_node, part_fci);
    current_base_pair=get_flat_structure_suffix(current_base_pair);
  } /* end while */  
}

void stochastic_backtrack_locally_optimal_structure_rec(int x, int y, plain_sequence * rna, char * structure, 
 double * energy, tr_node ** actual_node, TYPE * part_fci){
  TYPE r;
  TYPE partition_fci;
  TYPE e_contribution, local_energy; /* e_contribution - contribution of rna element (ML etc.) to energy */
  int i, j, i0;
  flat_cell * current_flat_structure;
  int current_base_pair;
  int thickness; 
  int alpha; /*number of helixes in pair*/
  int beta; /*number of unpaired bases*/
  TYPE intermediate1;
  TYPE intermediate2;
  
  intermediate1 = partition_function_table[x][y]/(*part_fci);
  intermediate2 = total_weight_par(*actual_node)*intermediate1;

  r = random_number(partition_function_table[x][y] - intermediate2);
  //printf("Upper_lim: %.30f\n", partition_function_table[x][y] - (total_weight_par(*actual_node)*partition_function_table[x][y])/(*part_fci));
  //printf("random r: %.30f, Partial fci val : %f, Energy ; %f table %f; weight %f \n", r, *part_fci, *energy, partition_function_table[x][y], total_weight_par(*actual_node));
  if (flat_table[x][y]!=NULL){
  	current_flat_structure=flat_table[x][y];
  	do {
  	  current_base_pair=current_flat_structure->current;
  	  i=get_flat_structure_start(current_base_pair);
  	  j=get_flat_structure_end(current_base_pair); 
      if ((i==x) && (j==y) && !( (i==1) && (j==rna->size) )){
		  e_contribution = stacking_energy(x-1,y+1);
          r -= (apply_energy_term(e_contribution, partition_function_table[x+1][y-1])
           - (partition_function_table[x][y]*tr_node_weight(current_flat_structure->current, *actual_node))/(*part_fci));
           //printf("substracting %f, total %.30f, r %.30f\n", tr_node_weight(current_flat_structure->current, *actual_node),  partition_fci - (partition_function_table[x][y]*tr_node_weight(current_flat_structure->current, *actual_node))/(*part_fci), r);
          if (r<0){
			add_base_pair(x,y,structure);
			*actual_node = add_if_nexists(current_flat_structure->current, *actual_node);
			*energy = apply_energy_term(e_contribution, *energy);
			*part_fci *= (apply_energy_term(e_contribution, partition_function_table[x+1][y-1])/partition_function_table[x][y]);
			//printf("Partial fci val : %f, Energy ; %f\n", *part_fci, *energy);
			stochastic_backtrack_locally_optimal_structure_rec(x+1, y-1, rna, structure, energy, actual_node, part_fci);
			return ;
          }
  	  }
  	  else{
		local_energy = 1.;
		partition_fci=1.;
		alpha = 0;
		beta = 0;
		i0 = x - 1;
		int *ptypes = (int *) malloc (sizeof(int)); /* table of pair types for multibranch loops */
		int *si = (int *) malloc (sizeof(int));
		int *sj = (int *) malloc (sizeof(int));
		ptypes[0] =  get_type(rna -> label[y+1], rna -> label[x-1]); /* closing pair type */ 
		si[0] = y;
		sj[0] = x;
  	    while(current_base_pair != 0){
		  alpha ++;
  	      i=get_flat_structure_start(current_base_pair);
  	      j=get_flat_structure_end(current_base_pair);
  	      ptypes =  (int *) realloc (ptypes, (alpha+1)*sizeof(int));
  	      si = (int *) realloc (si, (alpha+1)*sizeof(int));
	      sj = (int *) realloc (sj, (alpha+1)*sizeof(int));	      
	      ptypes[alpha] = get_type(rna -> label[i], rna -> label[j]); /* enclosed pair type(s) */
	      si[alpha] = i-1;
	      sj[alpha] = j+1;
		  beta += (i - i0 - 1);  
          i0 = j;
  	      thickness=get_BP(i,j); 
  	      if (thickness>MIN_HELIX_LENGTH) 
  		      thickness=MIN_HELIX_LENGTH;
  		  e_contribution = helix_energy(i,j,thickness);
  	      partition_fci = product(partition_fci, apply_energy_term(e_contribution, partition_function_table[i+thickness][j-thickness])); /* yann */
  	      local_energy = apply_energy_term(e_contribution, local_energy);
  	      // apply bonuses/penalties for external loop
  	      if(is_entirely_contained(x, y, rna)){
			partition_fci = apply_energy_term(ext_loop_energy(i, j, rna), partition_fci);
            local_energy = apply_energy_term(ext_loop_energy(i, j, rna), local_energy);
          }
  	      current_base_pair=get_flat_structure_suffix(current_base_pair);
  	    } /* end while */
  	    if((x != 1) || (y != rna->size)){  
	      if(alpha == 1){ /* internal loop or bulge since stack was taken care of before*/		  
			e_contribution =  internal_loop_energy(x-1, i, j, y+1, rna);
		  }else{ /*multiloop*/
			beta += (y - i0);
			e_contribution =  multiloop_energy(alpha, beta, ptypes, si, sj);
		  }	
          partition_fci = apply_energy_term(e_contribution, partition_fci);
          local_energy = apply_energy_term(e_contribution, local_energy);
		}
  	    /* Now we have the overall "weight" of flat structure, should we choose it? */
  	    r -= (partition_fci - (partition_function_table[x][y]*tr_node_weight(current_flat_structure->current, *actual_node))/(*part_fci));
  	    //printf("substracting %f, total %.30f, r %.30f\n", tr_node_weight(current_flat_structure->current, *actual_node), partition_fci - (partition_function_table[x][y]*tr_node_weight(current_flat_structure->current, *actual_node))/(*part_fci), r);
          if (r<0){
			 *energy = product(local_energy, *energy);
			 *part_fci *= (partition_fci/partition_function_table[x][y]);
			 //printf("Partial fci val : %f, Energy ; %f\n", *part_fci, *energy);
			 *actual_node = add_if_nexists(current_flat_structure->current, *actual_node);
             stochastic_backtrack_flat_structure_rec(current_flat_structure, rna, structure, energy, actual_node, part_fci);
             /* Warning: Potential memory leak below */
             return;
          }
          free(ptypes);
          free(si);
          free(sj); 	    
  	  }/* endif */
  	  current_flat_structure=current_flat_structure->next;
  	} while (current_flat_structure != NULL);
  }/*endif else */
  else{
	if((x != 1) || (y != rna->size))
		*energy = apply_energy_term(hairpin_energy(x-1, y+1, rna), *energy);
    return;
  }
}

/* Uses the stack. No memoization */
/* This function is not for you, Yann */
void count_all_locally_optimal_structures2(plain_sequence * rna){
  TYPE r=1; 
  init_flat_stack(rna);  
  while (succ_flat_stack(rna->size)==1){
    if ((stack_nb_of_bp*200)>=rna->size*MIN_PERCENT_PAIRING){
      r++ ; 
    }
  }
  printf("\n(debug): %.3f\n",r); 
  free_flat_table(rna); 
  free_flat_list(); 
}
