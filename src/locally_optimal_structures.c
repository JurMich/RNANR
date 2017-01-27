/**************************************************/
/*    locally optimal secondary structures        */
/**************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*#include <ViennaRNA/eval.h>*/
#include "rna.h"
#include "base_pairs.h"
#include "flat_structures.h"
#include "stack.h"

char *s; 

void init_pretty_string(plain_sequence * rna){
  s=(char*) malloc ((rna->size+1)*sizeof(char)); 
}


void  pretty_print_locopt_structure(plain_sequence * rna, FILE * outfile){
  int i; 
  int k, k_start, k_end ;
  int current_flat_structure;
  int x,y;
  if ((stack_nb_of_bp*200)>=rna->size*MIN_PERCENT_PAIRING){
    for (i=1; i<=rna->size; i++)
      s[i]='.';
    for (k=1;k<=stack_height; k++){
      current_flat_structure=(flat_stack[k]->order)->current; /* index of the current flat structure */
      x=flat_stack[k]->start;
      y=flat_stack[k]->end;
      do {
	k_start = get_flat_structure_start(current_flat_structure);
	k_end = get_flat_structure_end(current_flat_structure); 
	s[k_start]='(';
	s[k_end]=')';
	for (i=2; i<=flat_stack[k]->thickness; i++){
	  s[k_start+i-1]='[';
	  s[k_end-i+1]=']';
	}
	current_flat_structure=get_flat_structure_suffix(current_flat_structure); 
      } while (current_flat_structure>0) ;
    }/* end for */
    for (i=1; i<=rna->size; i++)
      fprintf(outfile, "%c", s[i]);
    /*char *tmp_rna = (char*) malloc(sizeof(char)*(rna->size+1));
    for(i=1;i<=rna->size;i++) 
      tmp_rna[i-1] = rna->label[i];
    char *tmp_s = (char*) malloc(sizeof(char)*(rna->size+1));
    for(i=1;i<=rna->size;i++)
	  tmp_s[i-1] = s[i];
    tmp_rna[rna->size] = '\0';
    tmp_s[rna->size] = '\0';
    float freenergy = vrna_eval_structure_simple (tmp_rna, tmp_s);
    free(tmp_s);
    free(tmp_rna);*/
    fprintf (outfile, "  --  %d\n", stack_nb_of_bp);
    /*fprintf(outfile, "  --  Free Energy : %.3f \n", freenergy);*/
  }/* end if */
}


void compressed_print_locopt_structure(plain_sequence * rna){
  int i; 
  int k, k_start, k_end ;
  int current_flat_structure;
  int x,y;
  
  if ((stack_nb_of_bp*200)>=rna->size*MIN_PERCENT_PAIRING){
    for (i=1; i<=rna->size; i++)
      s[i]='.';
    for (k=1;k<=stack_height; k++){
      current_flat_structure=(flat_stack[k]->order)->current; 
      x=flat_stack[k]->start;
      y=flat_stack[k]->end;
      do {
	k_start = get_flat_structure_start(current_flat_structure);
	k_end = get_flat_structure_end(current_flat_structure); 
	s[k_start]='(';
	s[k_end]=')';
	for (i=2; i<=flat_stack[k]->thickness; i++){
	  s[k_start+i-1]='[';
	  s[k_end-i+1]=']';
	}
	current_flat_structure=get_flat_structure_suffix(current_flat_structure); 
      } while (current_flat_structure>0) ;
    }/* end for */
    
    for (i=1; i<=rna->size; i++)
      printf("%c", s[i]);
    printf ("  --  %d\n", stack_nb_of_bp);
  }/* end if */
}

void generate_all_locally_optimal_structures(plain_sequence * rna, FILE *outfile){
  init_flat_stack(rna);  
  init_pretty_string(rna);
  pretty_print_locopt_structure(rna, outfile);
  while (succ_flat_stack(rna->size)==1){
    pretty_print_locopt_structure(rna, outfile);
  } 
  free_flat_table(rna); 
  free_flat_list(); 
}
