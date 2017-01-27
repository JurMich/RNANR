/**************************************************/
/*                      stack                     */
/**************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rna.h"
#include "base_pairs.h"
#include "flat_structures.h"
#include "stack.h"


stack_cell ** tmp_stack;

int number_of_base_pairs(int h){
 return get_flat_structure_nb_of_bp(flat_stack[h]->order->current)*flat_stack[h]->thickness;
}

void init_flat_stack_rec(int x, int y, int rna_size);

/* k: index of the flat structure in flat_list */
void init_flat_stack_index(int k, char thickness, int rna_size){
  if (k>0){
    init_flat_stack_rec(get_flat_structure_start(k)+thickness, get_flat_structure_end(k)-thickness, rna_size);
    init_flat_stack_index(get_flat_structure_suffix(k), thickness, rna_size);
  }
}

/* x..y : range on the RNA sequence */
/* Push the first flat structure for x..y */
void init_flat_stack_rec(int x, int y, int rna_size){
  
  if ((x<=0)||(y<=0)){
    printf("erreurs de valeurs: x=%d, y=%d\n", x, y); 
    exit(EXIT_SUCCESS);
  }

  if (flat_table[x][y]!=NULL){
    stack_height ++;
    flat_stack[stack_height]->start=x;
    flat_stack[stack_height]->end=y;
    flat_stack[stack_height]->order=flat_table[x][y];
    int i = get_flat_structure_start(flat_table[x][y]->current);
    int j = get_flat_structure_end(flat_table[x][y]->current);
     
    if ( (i==x) && (j==y)
    && !((x==1) && (y==rna_size))){
      flat_stack[stack_height]->thickness= 1;
    }
    else{
      flat_stack[stack_height]->thickness= MIN_HELIX_LENGTH;
    }
    stack_nb_of_bp += number_of_base_pairs(stack_height);
    init_flat_stack_index(flat_table[x][y]->current, flat_stack[stack_height]->thickness, rna_size);
  }
}


void init_flat_stack(plain_sequence * rna){
  int i; 
  tmp_stack=(stack_cell **) malloc((rna->size/2)*sizeof(stack_cell *));
  flat_stack=(stack_cell **) malloc((rna->size/2)*sizeof(stack_cell *));
  for (i=1; i<rna->size/2; i++){
    tmp_stack[i]= (stack_cell*) malloc(sizeof(stack_cell));
    flat_stack[i]= (stack_cell*) malloc(sizeof(stack_cell));
  }
  stack_height=0;
  stack_nb_of_bp=0;
  init_flat_stack_rec(1, rna->size, rna->size);
}


void display_stack(stack_cell ** stack, int s_height){
  int i; 

  for (i=s_height;i>=1; i--){
    printf("%d ", (stack[i]->order)->current);
    display_one_flat_structure( (stack[i]->order)->current);
    printf("\n");
  }
}


int succ_flat_stack(int rna_size){
  int tmp_stack_height=0;
  int i; 
  while( (stack_height>0) && ((flat_stack[stack_height]->order)->next==NULL)){
    while ((tmp_stack_height>0) && (flat_stack[stack_height]->end > tmp_stack[tmp_stack_height]->end)) {
      tmp_stack_height--;
    }
    tmp_stack_height++;       
    tmp_stack[tmp_stack_height]->start=flat_stack[stack_height]->start;
    tmp_stack[tmp_stack_height]->end=flat_stack[stack_height]->end;
    tmp_stack[tmp_stack_height]->order=flat_stack[stack_height]->order;
    tmp_stack[tmp_stack_height]->thickness=flat_stack[stack_height]->thickness;
    stack_nb_of_bp -= number_of_base_pairs(stack_height); 
    stack_height--;   
  }
  if (stack_height>0){
    /* on remplace l'element de tete de la pile flat_stack */
    stack_nb_of_bp-= number_of_base_pairs(stack_height); 
    flat_stack[stack_height]->order=(flat_stack[stack_height]->order)->next;
    int i = get_flat_structure_start((flat_stack[stack_height]->order)->current);
    int j = get_flat_structure_end((flat_stack[stack_height]->order)->current);
    
    if (  
	 (i==flat_stack[stack_height]->start ) 
	&& 
	 (j==flat_stack[stack_height]->end ) 
	&& !((flat_stack[stack_height]->start==1) && (flat_stack[stack_height]->end==rna_size)) ){
      flat_stack[stack_height]->thickness= 1;
    }
    else{
      flat_stack[stack_height]->thickness= MIN_HELIX_LENGTH;
    }

    stack_nb_of_bp+= number_of_base_pairs(stack_height); 
    while ((tmp_stack_height>0) && (flat_stack[stack_height]->end > tmp_stack[tmp_stack_height]->end)){
      tmp_stack_height--;
    }   
    init_flat_stack_index( flat_stack[stack_height]->order->current,  flat_stack[stack_height]->thickness, rna_size);
    for (i=tmp_stack_height; i>0; i--)
      init_flat_stack_rec(tmp_stack[i]->start, tmp_stack[i]->end, rna_size); 
    return 1; 
  }
  return 0; 
}




