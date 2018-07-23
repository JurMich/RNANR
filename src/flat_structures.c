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

/***********************************************/
/*                Flat structures              */
/***********************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rna.h"
#include "base_pairs.h"
#include "flat_structures.h"

long int MAX_MEMORY; 

char ** Useful_flat_structures; 
int index_max; 

typedef struct {int first_bp_index; int suffix; char nb_of_bp;} flat_struct;
flat_struct * flat_list; /* flat_list: list of all possible  flat structures */

flat_cell *** flat_table_aux; /****/

char *s; 


/**********************************/
/*        initialisation          */
/**********************************/ 

void init_flat_list(plain_sequence * rna){
  MAX_MEMORY= 1000000;
  flat_list= (flat_struct *) malloc (MAX_MEMORY* sizeof(flat_struct)); 
  index_max=0;
}

void free_flat_list(){
  free(flat_list);  
}

void init_flat_table(plain_sequence * rna){
  int x, y; 
  
  flat_table=(flat_cell ***) malloc ((rna->size+2)*sizeof(flat_cell**));
  for   (x=0; x<=rna->size+1; x++){
    flat_table[x]=(flat_cell **) malloc ( (rna->size+2)*sizeof(flat_cell *));
    for  (y=0; y<=rna->size+1; y++){
      flat_table[x][y]= NULL; 
    }
  }   
  flat_table_aux=(flat_cell ***) malloc ((rna->size+2)*sizeof(flat_cell**));
  for   (x=0; x<=rna->size+1; x++){
    flat_table_aux[x]=(flat_cell **) malloc ((rna->size+2)*sizeof(flat_cell *));
    for  (y=0; y<=rna->size+1; y++){
      flat_table_aux[x][y]= NULL; 
    }
  }  
}


void free_flat_structure(int x, int y){
 flat_cell * ptr1, * ptr2; 
       ptr1=flat_table[x][y]; 
      while (ptr1!=NULL){
 	ptr2=ptr1->next;
	free(ptr1);
	ptr1=ptr2;
      }
      flat_table[x][y]=NULL; 
}

void free_flat_table(plain_sequence * rna){
  int x,y; 
  flat_cell * ptr1, * ptr2; 
  for (x=1; x<=rna->size+1; x++){
    for (y=1; y<=rna->size+1; y++){
      ptr1=flat_table[x][y]; 
      while (ptr1!=NULL){
 	ptr2=ptr1->next;
	free(ptr1);
	ptr1=ptr2;
      }/* end while */
    }/*end for y*/
    free(flat_table[x]);
  }/*end for x */
  free(flat_table); 
}


void init_useful_flat_structures(plain_sequence * rna){
  int x, y, k;
  
   Useful_flat_structures=(char **)malloc((rna->size+2)*sizeof(char *));
  
  for (x=1; x<= rna->size; x++){
    Useful_flat_structures[x]= (char *) malloc((rna->size+2)*sizeof(char));
    for (y=1; y<=rna->size; y++)
      Useful_flat_structures[x][y]=0; 
  }

  /* step 1: flat structures that are covered by a helix  */
  for (x=1; x<=MAXI(1,rna->size-MIN_HELIX_LENGTH); x++)
    for (y=x+MIN_HELIX_LENGTH; y<=rna->size; y++)
      if (get_BP(x,y)>=MIN_HELIX_LENGTH){
	for (k=MAXI(1,y-MIN_HELIX_LENGTH - MAX_BULGE_SIZE); k<=MAXI(1,y-MIN_HELIX_LENGTH); k++)
	  if ((k>0) && (k<=rna->size)) /* a verifier */
	    Useful_flat_structures[x+MIN_HELIX_LENGTH][k]=2;
      }
  
  for (k=MAXI(1,rna->size - MAX_BULGE_SIZE); k<=rna->size ; k++)
    Useful_flat_structures[1][k]=2;
  
  /* step 2 : all suffixes */
  for (x=1; x<rna->size; x++)
    for (y=x; y<=rna->size; y++)
      if (Useful_flat_structures[x][y]==2)
	for (k=x+1; k<y; k++)
	  if  (Useful_flat_structures[k][y]==0)
	    Useful_flat_structures[k][y]=1;
}

flat_struct get_flat_structure(int k){
  return flat_list[k]; 
}

int get_flat_structure_start(int k){
  flat_struct f=get_flat_structure(k); 
  int x= BP_list[f.first_bp_index].final; 
  int p= BP_list[f.first_bp_index].span;
  return x-p; 
}

/* k is an index in  flat_list */
int get_flat_structure_end(int k){
  flat_struct f=get_flat_structure(k); 
  return  BP_list[f.first_bp_index].final; 
}

int get_flat_structure_suffix(int k){
  flat_struct f=get_flat_structure(k); 
  return f.suffix; 
}

int  get_flat_structure_nb_of_bp(int k){
  flat_struct f=get_flat_structure(k); 
  return f.nb_of_bp; 
}


/*  Add the new flat structure (x,y,k) to flat_list, and output   */
/*  its index. This new element is added at the end of the    */
/*  list.                                                     */
int get_new_flat_index(int x, int y, int k){
  flat_struct current_flat_structure;
  int current_bp_index= get_BP_index(x,y);  /* index of the base pair (x,y) */ 
  /* creation of a new flat structure  */
  current_flat_structure.first_bp_index=current_bp_index; 
  current_flat_structure.suffix=k;
  current_flat_structure.nb_of_bp= get_flat_structure(k).nb_of_bp+1;
  index_max++;
  if (index_max==MAX_MEMORY){
    MAX_MEMORY = 2*MAX_MEMORY;
    flat_list=realloc(flat_list, MAX_MEMORY*sizeof(flat_struct));
    if (flat_list==NULL)
       fprintf(stderr, " Out of memory - too many flat structures. No more memory available\n");
  }
  flat_list[index_max]=current_flat_structure;
  return index_max;
}

/******************************/
/*         display            */
/******************************/

void display_one_flat_structure(int k){

  if (k>0){
    printf("(%d,%d)", get_flat_structure_start(k), get_flat_structure_end(k)); 
    display_one_flat_structure(flat_list[k].suffix);
  }
}

void display_all_flat_structures(){

  int k; 
  printf("\n%d flat structures found\n\n", index_max);
  printf("first bp index : %d\n", get_flat_structure(0).first_bp_index);
  for (k=1; k<=index_max; k++){
    printf("%d -> ", k); 
    display_one_flat_structure(k);
    printf("\n");
  }
}

void display_flat_range(int x, int y){
  
  flat_cell * ptr=flat_table[x][y];

  if (ptr!=NULL){
    printf("[%d..%d] -> ", x, y);   
    while (ptr!= NULL) {
      display_one_flat_structure(ptr->current);    
      printf("|");
      ptr=ptr->next;
    }
    printf("\n"); 
  }
}

void display_all_flat_structures_range(plain_sequence * rna){
  int x,y; 

  printf("-- Set of all possibles maximal flat structures\n\n");
 for (x=1; x<=MAXI(1,rna->size-2*MIN_HELIX_LENGTH-MIN_LOOP_SIZE+1); x++)
    for (y=x+2*MIN_HELIX_LENGTH+MIN_LOOP_SIZE-1; y<=rna->size; y++)
      display_flat_range(x,y);
}

/* TRUE, if the flat structure with index k is conflicting with all base pairs of the form $(x,i)$ 
 with $x<i<=y$. Additional hypotheses: k does not equal 0, and y is the last position of k and is paired in k.  */
int left_extension(int x, int y, int k, plain_sequence* rna){
 
  int i; 
  int start_pos=get_flat_structure_start(k);
  int end_pos=get_flat_structure_end(k) ; 
 
  /* Case 1: one can create a base pair starting at position x that extends */
  /* the first helix of the flat structure k : [(((..)))].(...)(...) */
  if (
      (start_pos==x+1) 
      && (end_pos<y) 
      && (get_BP(x,end_pos+1)>0) 
      && (get_flat_structure_start(get_flat_structure_suffix(k))>end_pos+1)
      )
    return 0;
  /* End case 1 */
  
  /* Case 2: it is possible to create a brand new helix [[[..]]] starting at */
  /* position x and compatible with k: [[[.]]].(((..))).(((..))) or          */
  /* [[[.(((..))).(((..)))]]]..(((...)))                                     */
  if (start_pos-x >=MIN_HELIX_LENGTH) {
    start_pos=x+2*MIN_HELIX_LENGTH+MIN_LOOP_SIZE-1;  
    do{
      end_pos=get_flat_structure_start(k)-1;  
      for (i=start_pos; i<=end_pos;  i++){
	if (get_BP(x,i)>=MIN_HELIX_LENGTH) 
	  return 0; 
      } 
      start_pos=get_flat_structure_end(k)+MIN_HELIX_LENGTH;
      k=get_flat_structure_suffix(k);
    }while (k!=0); 
  }
  /* End case 2 */

  return 1; 
}


/* TRUE, if the flat structure with index k is in conflict with all base pairs of the form       */
/* (i,y), x<=i<y, with the additional information that the last paired position in k is q (q<y) */
/* More precisely  :                                                                            */
/* 0: the flat structure k is not selected                                                      */
/* 1: the flat structure k is added to flat_table                                               */
/* 2: the flat structure k is added to flat_table_ext                                           */
/*                                                                                              */
/* Prerequisite: k != 0                                                                         */
int right_extension(int x, int y, int k, int q, plain_sequence * rna){
  int i; 
  int k_old=0; 
  int k_old_old=0;
  int start_pos, end_pos; 
  int k0=get_flat_structure_start(k);  

  /* Case 1: a brand new helix can be created from position y (y-q>=MIN_HELIX_LENGTH), */
  /* We try each possible base pair between i and y, to search for non-conflicting i's */ 

  /* Case 1.a: the helix starts after position q : .(((..))).(((..))).[[[.]]]*/
  for (i=q+1; i<=y-2*MIN_HELIX_LENGTH-MIN_LOOP_SIZE+1; i++)
    if (get_BP(i,y)>=MIN_HELIX_LENGTH) 
      return 0;
  
  /* Case 1.b: the helix starts before position q*/
  start_pos=x;
  if (y-q>=MIN_HELIX_LENGTH) {
    do{
      end_pos=get_flat_structure_start(k)-MIN_HELIX_LENGTH; 
      for (i=start_pos; i<=end_pos;  i++)
	if (get_BP(i,y)>=MIN_HELIX_LENGTH) 
	  return 0; 
      start_pos=get_flat_structure_end(k)+1;
      k_old_old=k_old; 
      k_old=k; 
      k=get_flat_structure_suffix(k);
    }while (k!=0); 

    /* we know for sure that start_pos is the position after the last base pair. */
    /* We investigate [start_pos.. XXX] */  
    for (i=start_pos; i<=y-2*MIN_HELIX_LENGTH-MIN_LOOP_SIZE; i++)
      if (get_BP(i,y)>=MIN_HELIX_LENGTH) 
	return 0; 
  }/* end if (y-q) */
  
  /* End of case 1 */

  /*  Case 2: position y can be paired to form a base pair that extends the last */
  /* helix of the flat  structure */ 
  /* We first reach the last base pair of the flat structure. */
  /* This loop is used  only if  !(y-q>=MIN_HELIX_LENGTH) */
  if (q==y-1){
    while (k!=0){
      k_old_old=k_old;
      k_old=k;
      k=get_flat_structure_suffix(k);
    }/* end while */
    
    start_pos=get_flat_structure_start(k_old);
    
    if (k_old_old==0) {
      /* the flat structure contains a single helix that starts at position x*/ 
      if (start_pos==x) 
	return 1; 
      else 
	if (get_BP(start_pos-1,y)>0) 
	  return 0; 
    }
    
    /* from this point, the flat structure contains at least two helices */
    if (k_old_old>0){
      if (start_pos>get_flat_structure_end(k_old_old)+1) {
	if (get_BP(start_pos-1,y)>0) 
	  return 0; 
      }
    }
  }
  /* End of case 2 */

  /* Case 3: x and y can be paired together */
  /* to extend an upper helix  */
  if (get_BP(x,y)==0) 
    return 1;   
  /* we know for sure that q<y  */   
  if ((k0>x) && (y<rna->size))
    return 2; /* the structure is added to flat_table_ext */
  else 
    return 1; /* the structure is added to flat_table */
  /* End of case 3 */
}


/* Insert the flat structure of index k before the cell designated by ptr. */
/* This cell is either in flat_table or in flat_table_ext                  */
flat_cell * push(flat_cell * ptr, int k){
  flat_cell * new_ptr;
  new_ptr= (flat_cell *) malloc( sizeof(flat_cell));
  if  (new_ptr ==NULL){
    fprintf(stderr, "\n Out of memory \n");
    exit(EXIT_FAILURE);
  } 
  new_ptr->current = k; 
  new_ptr->next = ptr; 
  return new_ptr; 
}


int is_useful(int x, int y, int thickness){
  int i;
  /* we are below the MAX_LOOP_SIZE threshold */
  if (y-x < MAX_LOOP_SIZE + (2*thickness)) 
    return 1;
  if (thickness==1)
    return ((flat_table[x+1][y-1]!=NULL) || (flat_table_aux[x+1][y-1]!=NULL));
  for (i=MIN_HELIX_LENGTH; i<=thickness; i++){
     if (flat_table[x+i][y-i]!=NULL) return 1;
     if (flat_table_aux[x+i][y-i]!=NULL) return 1;
  }
  return 0;
}

/* (x,p) k */
/* y is the last position in k, and is paired */
void create_new_flat_structures(int x,int p,int y,int k,plain_sequence * rna){
  int current_flat_structure; 
  int status_left, status_right;
  int z, t;
  int maxt, maxtt, minz;

  current_flat_structure= get_new_flat_index(x,p,k);
  flat_table[x][y]=push(flat_table[x][y],current_flat_structure);
  
  /* right extensions with no left extension */
  t=y+1;
  status_right=1;
  maxt=MIN(rna->size,y+MAX_BULGE_SIZE);
  maxtt=maxt;
  while ((maxtt>y) && (Useful_flat_structures[x][maxtt]!=2))
     maxtt--; 
  while ((t<=maxtt) && (status_right!=0)) {
    status_right=right_extension(x,t,current_flat_structure, y, rna); 
    if (status_right>0)
	flat_table_aux[x][t]=push(flat_table_aux[x][t], current_flat_structure);
    t++;
    } 
  if (status_right==0)
    maxt=t--;
    
  minz=MAXI(1,x-MAX_BULGE_SIZE);
  z=x-1;
  status_left=1;
  while (z>=minz && status_left!=0){
    status_left=left_extension(z,y,current_flat_structure, rna);
    if (status_left>0){
      flat_table[z][y]=push(flat_table[z][y],current_flat_structure);
      status_right=1;
      t=y+1;
      maxtt=maxt;
      while ((maxtt>y) && (Useful_flat_structures[z][maxtt]!=2))
	maxtt--; 
      while (t<=maxtt && status_right!=0) {
	status_right=right_extension(z,t,current_flat_structure, y, rna); 
	if (status_right==1)
	  flat_table_aux[z][t]=push(flat_table_aux[z][t], current_flat_structure);
	  t++;
	  } 
      if (status_right==0)
	  maxt=t--;
    }
    z--;
  } /* end while z */
	  
}      


/* calcul du tableau flat_table */
void build_all_flat_structures( plain_sequence * rna){
  
  int x,y,p ;  
  flat_cell * ptr;
  
  /* flat_table[i][j]: list of all flat structures for the range [i..j] */
  /* Each flat structure is characterized by an index from flat_list */
  
  init_flat_list(rna); 
  init_flat_table(rna);
  init_useful_flat_structures(rna); 
  
  for (x=MAXI(1,rna->size-MIN_HELIX_LENGTH+1); x>0; x--){ 
    for (p=x+MIN_LOOP_SIZE; p<=rna->size; p++){
      /* helix internal extension */
      if ((get_BP(x,p)<MIN_HELIX_LENGTH) && (get_BP(x,p)>0) && is_useful(x,p,1)){
	  flat_table_aux[x][p]= push(flat_table_aux[x][p], get_new_flat_index(x,p,0));
      }/* end if */
      /* creation of all flat structures starting with the base pair (x,p) */
      if ((get_BP(x,p)>=MIN_HELIX_LENGTH) && is_useful(x,p,get_BP(x,p))) {
	/* we create the flat structure {(x,p)} and all derived flat structures*/
	create_new_flat_structures(x,p,p,0,rna); 
	/* we create all other flat structures starting with (x,p) */  
	for (y=p+2*MIN_HELIX_LENGTH+MIN_LOOP_SIZE; y<=rna->size; y++){
	  if (Useful_flat_structures[x][y]>0){
	    ptr=flat_table[p+1][y];
	    if (ptr!=NULL) { 
	      while (ptr != NULL){ /*  /!\  */
		if  (get_flat_structure_nb_of_bp(ptr->current)<MAX_DEGREE){
		  create_new_flat_structures(x,p,y,ptr->current,rna); 
		}
		ptr=ptr->next;
	    }/* end while */
	  } /* end if */
	}/* end if useful */ 
      }/* end for y */
      }/* end if p */
    }/* end for p */
  } /* end for x */
  
  
  for (x=rna->size; x>0; x--){ 
    for (y=x+1; y<=rna->size; y++){
      ptr=flat_table[x][y];
      if (ptr==NULL){
	flat_table[x][y]=flat_table_aux[x][y];
      }
      else{
	while (ptr->next!=NULL)
	  ptr=ptr->next;
	ptr->next=flat_table_aux[x][y];
      }
    }
  }
} 
  
