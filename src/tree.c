/**************************************************/
/*                     tree                       */
/**************************************************/

#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "tree.h"

// creates new tree node with undefined weight with parent node
tr_node * create_tr_node(int id, tr_node * parent){
  tr_node * new_tr_node = (tr_node *) malloc(sizeof(tr_node));
  new_tr_node->id = id;
  new_tr_node->parent = parent;
  new_tr_node->ll_start = NULL;
  new_tr_node->weight = 0;
  return new_tr_node;
}

// creates new node for linked list
ll_node * create_ll_node(int id){
  ll_node * new_ll_node = (ll_node *) malloc(sizeof(ll_node));
  new_ll_node->id = id;
  new_ll_node->ll_next = NULL;
  new_ll_node->next_node = NULL;
  return new_ll_node;	
}

// creates root (start of a tree)
tr_node * create_root(){
  tr_node * root = create_tr_node(-1, NULL);
  return root;  	
}

// inserts a new tr_node (along with creating ll_node) into tree
tr_node * insert_tr_node(int id, ll_node * ptr_prev, ll_node * ptr, tr_node * last_node){
  // create new ll_node	 
  ll_node * new_ll_node = create_ll_node(id);
  // relinking of previous and next node if existing(reordering)	
  if(ptr_prev) 
     ptr_prev->ll_next = new_ll_node; 
  else
     last_node->ll_start = new_ll_node;
  
  if(ptr){ 
    new_ll_node->ll_next = ptr;
  }
  
  // now create joined tr_node	  
  tr_node * new_tr_node = create_tr_node(id, last_node);
  new_ll_node->next_node = new_tr_node;
  return new_tr_node;
}

// creates node with id if such node doesn't already exist for last_node
tr_node * add_if_nexists(int id, tr_node*last_node){
  ll_node *ptr = last_node->ll_start;
  ll_node *ptr_prev = NULL;
  while(ptr){
	int id_tmp = ptr->id; 
	if(id_tmp == id){
	  return ptr->next_node;
	}
	  
	if(id_tmp>id){
	  tr_node * new_tr_node = insert_tr_node(id, ptr_prev, ptr, last_node);
	  return new_tr_node;    
	}
	ptr_prev = ptr;
	ptr = ptr->ll_next;
  }
  // if here then reached end of ll -> add node to its end
  tr_node * new_tr_node = insert_tr_node(id, ptr_prev, NULL, last_node);
  return new_tr_node;  	
}

// returns weight of an element if it exists; otherwise returns 0
TYPE tr_node_weight(int id, tr_node * last_node){
  ll_node *ptr = last_node->ll_start;
  while(ptr){
    if(ptr->id == id)
      return ptr->next_node->weight;
    ptr = ptr->ll_next;
  }
  return 0.;
}

// returns sum of all weights for given parent
TYPE total_weight_par(tr_node * last_node){
  ll_node *ptr = last_node->ll_start;
  TYPE weight_sum = 0.;
  while(ptr){
    weight_sum += ptr->next_node->weight;
    ptr = ptr->ll_next;
  }
  return weight_sum;
}

// updates weight of a given node
void update_weight(tr_node * node, TYPE weight){
  node->weight += weight;
}

// tracebacks to root while updating values for each node passed
tr_node * traceback_to_root(tr_node * leaf, TYPE struct_weight){
  update_weight(leaf, struct_weight);
  while(leaf->parent){
	update_weight(leaf->parent, struct_weight);
    leaf = leaf->parent;
  }
  return leaf; 	
}