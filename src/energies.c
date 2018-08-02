/*
* This file constitutes a part of RNANR program.
*
* RNANR is free software to which the terms of the GNU Affero General 
* Public License apply as published by the Free Software Foundation, 
* either of version 3, or (at your option) any later version.
*
* This program is distributed without ANY WARRANTY or RESPONSIBILITY.
* You can distribute and change its contents, youŕe not allowed to SELL it.
* For more precisions, see the GNU Affero General Public License.
*
* You should have received a copy of the GNU Affero General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

/**************************************************/
/*                    energies                    */
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
#include <ViennaRNA/loops/all.h>
#include "rna.h"
#include "base_pairs.h"
#include "flat_structures.h"
#include "sort.h"
#include "stack.h"
#include "energies.h" 

/* Contains functions used in mostly in counting.h to compute energies of
 * different structural elements among others.
 */

/**********    end operators    ***********/


/********************************************/
/*       structure and misc. utilities     */
/******************************************/
void print_structure(char * structure, double part_fci, float energy_ref, plain_sequence * rna, FILE * outfile){
  int i;
  if(!outfile)
	outfile = stdout;
  for (i=1; i<=rna->size; i++)
    fprintf(outfile, "%c", structure[i]);
  fprintf(outfile, " - Free Energy : %.3f \n", (-RT*log(part_fci))/100);
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

/* checks if interval [i,j] covers the entirety of sequence */
int is_entirely_contained(int i, int j, plain_sequence * rna){
  return ((i==1) && (j==rna->size));
}


/************ Free energy utilities ************/

/* Energy of (i+1,j-1) stacking over (i,j), i+1<j-1 */
double stacking_energy(int i, int j){
  double stack_energy;
  stack_energy = vrna_E_stack(E_fold_cp, i, j);  /*ViennaRNA*/
  return stack_energy;
}

/* Energy of helix spanning from (i,j) to (i+length,j-length), i+length<j-length */
double helix_energy(int i, int j, int length){
  double dG;
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
double multiloop_energy(int alpha, int beta, int *ptypes, int *si, int *sj){
  	int i, closing_pty, unpaired_pty; /*to extract from Vienna package*/
  	double mtloop_E;
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
double internal_loop_energy(int i, int k, int l, int j, plain_sequence *rna){
  double it_energy;
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
double hairpin_energy(int i, int j, plain_sequence * rna){
  double hairpin_energy;
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
double ext_loop_energy(int i, int j, plain_sequence *rna){
  double ext_loop_energy;
  int ptype, seq_i, seq_j;
  ptype = get_type(rna->label[i], rna->label[j]);
  seq_i = (i > 1) ? E_fold_cp->sequence_encoding[i-1] : -1;
  seq_j = (j < rna->size) ? E_fold_cp->sequence_encoding[j+1] : -1;
  ext_loop_energy = E_ExtLoop(ptype, seq_i, seq_j, E_fold_cp->params); 
  return ext_loop_energy;	
}

/***********   end utilities     ************/

void init_locopt_table(int n){
  int x, y; 
  number_of_locopt_structures=(long int **) malloc ((n+1)*sizeof(long int *));
  for (x=0; x<=n; x++){
    number_of_locopt_structures[x]=(long int *) malloc ((n+1)*sizeof(long int));
    for (y=0; y<=n; y++){
      number_of_locopt_structures[x][y]= 0;/* yann */
    }
  } 
} 


long int count_all_locally_optimal_structures(plain_sequence * rna){
  int n, x, y, i, j;
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
        n=1; /* yann */
	    while(current_base_pair != 0){
	      i=get_flat_structure_start(current_base_pair);
	      j=get_flat_structure_end(current_base_pair);
	      get_flat_structure_suffix(current_base_pair);
	      thickness=get_BP(i,j); 
	      if (thickness>MIN_HELIX_LENGTH) 
		thickness=MIN_HELIX_LENGTH;
	      n *= number_of_locopt_structures[i+thickness][j-thickness]; /* yann */
	      current_base_pair=get_flat_structure_suffix(current_base_pair);
	    } /* end while */
	  }/* endif */
	  number_of_locopt_structures[x][y] += n; /* yann */
	  current_flat_structure=current_flat_structure->next;
	}while (current_flat_structure != NULL);
      }/*endif else */
      else{
	number_of_locopt_structures[x][y]=1; 
      }
    }/* end for y */
  }/* end for x */
  printf("\nNumber of structures: %ld\n", number_of_locopt_structures[1][rna->size]); 
  return number_of_locopt_structures[1][rna->size];
}


/* stacks every flat structure into a pile and counts them */
long int count_all_flat_structures(plain_sequence * rna){
	long int x, y, i, j;
	long int n_flat = 0;
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

/* calculates energy of a structure using RNAfold as control
 * to be compared with energies returned by RNANR
 * (used for testing/debuffing purposes).
 */
 
float get_reference_energy(vrna_fold_compound_t *E_fold_cp, char * structure, plain_sequence * rna){
	char *fold = (char*) malloc(sizeof(char)*(rna->size+1));
    strncpy(fold, structure+1, rna->size);
    fold[rna->size] = '\0';
	float energy = vrna_eval_structure (E_fold_cp, fold);
	free(fold);
	return energy;
}

/***** light version of counting (slower) ******/

/* Uses the stack. No memoization */
/* This function is not for you, Yann */
void count_all_locally_optimal_structures2(plain_sequence * rna){
  long int r=1; 
  init_flat_stack(rna);  
  while (succ_flat_stack(rna->size)==1){
    if ((stack_nb_of_bp*200)>=rna->size*MIN_PERCENT_PAIRING){
      r++ ; 
    }
  }
  printf("\n(debug): %ld\n",r); 
  free_flat_table(rna); 
  free_flat_list(); 
}

