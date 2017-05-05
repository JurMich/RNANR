/**************************************************/
/*                     sort                       */
/**************************************************/


#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "rna.h"
#include "flat_structures.h"
#include "sort.h"

// creates an empty linktab
linktab * start_linktab(){
	linktab * linktable = (linktab *) malloc(sizeof(linktab));
	linktable->links = malloc(0);
	linktable->length = 0;
	return linktable;	
}

// adds one element to linktab
void add_link_element(linktab * linktable, double Boltz_partition, flat_cell * fcell){
	flat_sort new_sort;
	new_sort.Boltz_partition = Boltz_partition;
	new_sort.fcell = fcell;
	linktable->links = (flat_sort *) realloc(linktable->links, (linktable->length+1)*sizeof(flat_sort));
	linktable->links[linktable->length] = new_sort;
	linktable->length++;
}

/* creates a new tab from part of old tab, needs previous two functions.
 * No externally visible, just makes sorting function easier
 */
linktab * copytab_part(linktab * linktable, int start, int end){
	linktab * copytab;
	copytab = start_linktab();
	if ((linktable->length<start)||(linktable->length<end)){
		return copytab;  // return empty table
	}
	for(int i = start; i<end; i++){
		add_link_element(copytab, linktable->links[i].Boltz_partition, linktable->links[i].fcell);	
	}
	return copytab;
}

// destructor
void free_linktab(linktab * linktable){
	free(linktable->links);
	free(linktable);
}

/*#####################
 *##### Mergesort #####
 *#####################
 */
 
// merging part of sorting function
linktab * merge_sort(linktab * first, linktab * second){
	linktab * linktable;
	int curr_p1, curr_p2;
	linktable = start_linktab();
	curr_p1 = 0;
	curr_p2 = 0;
	// compare values within tables until we ran out of tables
	while((curr_p1<first->length)&&(curr_p2<second->length)){
		if(first->links[curr_p1].Boltz_partition>second->links[curr_p2].Boltz_partition){
			add_link_element(linktable, first->links[curr_p1].Boltz_partition, first->links[curr_p1].fcell);
			curr_p1++;	
		}else{
			add_link_element(linktable, second->links[curr_p2].Boltz_partition, second->links[curr_p2].fcell);
			curr_p2++;
		}
			
	}
	// now add what stayed out of tables
	if(curr_p1<first->length){
		for(int a = curr_p1; a<first->length; a++){
			add_link_element(linktable, first->links[a].Boltz_partition, first->links[a].fcell);
		}
	}	
	
	if(curr_p2<second->length){
		for(int a = curr_p2; a<second->length; a++){
			add_link_element(linktable, second->links[a].Boltz_partition, second->links[a].fcell);
		}
	}
	
	return linktable;
}

// sorting function
linktab * sort(linktab * linktable){
	int halfway;
	linktab * first_half;
	linktab * second_half;
	
	// list of length <= 1 is always sorted, and we want to shut process if table is empty
	if(linktable->length <= 1){
		return linktable;
	}
	
	// split table in half
	halfway = (int)linktable->length/2;
	first_half = copytab_part(linktable, 0, halfway);
	second_half = copytab_part(linktable, halfway, linktable->length);
	
	// do recursively the same to two halves
	first_half = sort(first_half);
	second_half = sort(second_half);
	
	// sort and merge two halves
	linktable = merge_sort(first_half, second_half);
	
	free_linktab(first_half);
	free_linktab(second_half);
	
	return linktable;	
}