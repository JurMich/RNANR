#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rna.h"
#include "base_pairs.h"
 

/*****************************************************************************/
/*                                  display                                  */
/*****************************************************************************/


 /* BP_table[x][y]->thickness = 1 : x-y and x can be paired together  */
 /* x: final position of the pairing */
 /* y: length of the pairing */
typedef struct {char thickness; int index;} bp_cell; 
bp_cell *** BP_table; 

void display_BP_list(){
  int k; 
  printf("\n--- List of all possible base pairs t----------\n"); 
  for (k=1; k<=BP_list[0].final;k++){
    printf("%d -> (%d-%d) [%d nt]\n",k, BP_list[k].final- BP_list[k].span, BP_list[k].final,  BP_table[ BP_list[k].final][BP_list[k].span]->thickness );
  }
  printf("\n-------\n%d base pairs found  \n--------------\n", BP_list[0].final); 
  
}

void display_BP_table(plain_sequence * rna){
  int x,y; 
  for (x=1; x<=rna->size; x++){
    for (y=x+1; y<=rna->size; y++){
      if (BP_table[y][y-x]->thickness>0)
	printf("(%d,%d): thickness %d, index %d\n",x,y, BP_table[y][y-x]->thickness, BP_table[y][y-x]->index );
    }
  }
}

void display_base_pairs(plain_sequence *rna){
  int x,y; 
  printf("-- Set of all possible base pairs \n\n"); 
  for (x=1; x<=rna->size; x++)
    for (y=x+1; y<=rna->size; y++)
      if (BP_table[y][y-x]->thickness>0)
	printf("(%d,%d) ", x,y);
  printf("\n"); 
}

void init_BP_external(plain_sequence *rna, char * BP_file_name){
  int sequence_length=rna->size; 
  FILE * BP_file; 
  int * positions; 
  int position_top; 
  char *s; 
  int i; 
  
  positions=(int *) malloc (sizeof(int)*sequence_length/2); 
  s = (char*) malloc (sizeof(char)*sequence_length); 
  BP_file=fopen(BP_file_name, "rt");
  if(BP_file==NULL){
    printf("error load : %s\n",BP_file_name);
    exit(EXIT_FAILURE);
  }
  
  while (!feof(BP_file)){
    fgets(s,sequence_length,BP_file);
    position_top=0; 
    for (i=0; i<sequence_length; i++){
      switch (s[i]) {
      case '(': 
	position_top++; 
	positions[position_top]=i; 
	break; 
      case  ')': 
	BP_table[i][i-positions[position_top]]->thickness =1; 
	position_top--;
	break; 
      }
    }
  }
  fclose(BP_file);
  free(positions);
}

void init_BP_ACGU( plain_sequence * rna){
  
  int x, y; 

 for (x=1; x<=rna->size; x++)
   for (y=MIN_LOOP_SIZE+1; y<MIN(x,MAX_HELIX_SCOPE+1); y++)
     if (are_complementary(rna->label[x], rna->label[x-y]) )
       BP_table[x][y]->thickness=1;        
}


 void compute_all_base_pairs(plain_sequence *rna, int BPFILE, char * BP_file_name){
   
   int x,y;
   int nb_BP = 0;
   BP_table=(bp_cell ***)malloc((rna->size+2)*sizeof(bp_cell **)); 
   for   (x=0; x<=rna->size+1; x++){
     BP_table[x]= (bp_cell **) malloc((rna->size+2)*sizeof(bp_cell*));
     for (y=0; y<=rna->size+1; y++){
       BP_table[x][y]=(bp_cell *) malloc(sizeof(bp_cell));
       BP_table[x][y]->thickness=0;
       BP_table[x][y]->index=0; 
     }
   }

   /* mark all possible base pairs with 1 */
   switch (BPFILE){
   case  0:   
     init_BP_ACGU (rna);  
     break; 
   case 1: init_BP_external(rna, BP_file_name); 
     break; 
   case 2:   init_BP_ACGU (rna);
     init_BP_external(rna, BP_file_name); 
     break; 
   } 

   /* compute the thickness of each base pair */
   for (x=1; x<=rna->size; x++)
     for (y=MIN_LOOP_SIZE+1; y<MIN(x,MAX_HELIX_SCOPE+1); y++)
       if (BP_table[x][y]->thickness>0)
	 BP_table[x][y]->thickness=BP_table[x-1][y-2]->thickness+1;
   
   
   /* remove all base pairs that do not satisfy the MIN_HELIX_LENGTH threshold*/
   for (x=1; x<=rna->size; x++){
     for (y=MIN_LOOP_SIZE+1; y<MIN(x, MAX_HELIX_SCOPE+1); y++){
       if  (BP_table[x][y]->thickness<MIN_HELIX_LENGTH){
		   
	 if  ((x+MIN_HELIX_LENGTH>rna->size) || (y+MIN_HELIX_LENGTH>=x)){
	   BP_table[x][y]->thickness=0; 
	 }
	 else if (y+2*MIN_HELIX_LENGTH < MAX_HELIX_SCOPE + 1){
       if (BP_table[x+MIN_HELIX_LENGTH][y+2*MIN_HELIX_LENGTH]->thickness<MIN_HELIX_LENGTH){
	     BP_table[x][y]->thickness=0; 
       }
	 }
	 
       }
     }
   }
   
   /* remove all base pairs that do not satisfy the MIN_HELIX_LENGTH threshold 
   for (x=1; x<=rna->size; x++){
     for (y=MIN_LOOP_SIZE+1; y<MIN(x, MAX_HELIX_SCOPE+1); y++){
       if  (BP_table[x][y]->thickness<MIN_HELIX_LENGTH){
		 BP_table[x][y]->thickness=0; 
       }
     }
   } */
   
   /* compute the index of each base pair and store it in BP_list */
   for (x=1; x<=rna->size; x++){
     for (y=MIN_LOOP_SIZE+1; y<MIN(x, MAX_HELIX_SCOPE+1); y++){
       if  (BP_table[x][y]->thickness>0)
	 nb_BP++; 
     }
   }
   BP_list= (bp *) malloc ((nb_BP+1)*sizeof(bp)); 
   BP_list[0].final=nb_BP;
   nb_BP=1;  
   for (x=1; x<=rna->size; x++){
     for (y=MIN_LOOP_SIZE+1; y<MIN(x, MAX_HELIX_SCOPE+1); y++){
       if  (BP_table[x][y]->thickness>0){
	 BP_list[nb_BP].final=x;
	 BP_list[nb_BP].span=y;
	 BP_table[x][y]->index=nb_BP; 
	 nb_BP++; 
       }
     }  
   }
 /*display_BP_list(); */
 /*display_BP_table(rna);*/
 }


void free_base_pairs(plain_sequence * rna){

  int x; 
  for   (x=1; x<=rna->size; x++){
    free(BP_table[x]);
  }
  free(BP_table); 
  free(BP_list);
}

/***********************************/
/*               get               */
/***********************************/

/* retourne la longueur de la tige issue des positions i et j (i<j)
   0 s'il n'y a pas d'appariements entre i et j*/ 
int get_BP(int i, int j){
  return BP_table[j][j-i]->thickness;
}


int get_BP_index(int i, int j){
 return BP_table[j][j-i]->index;
}

