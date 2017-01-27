#include "rna.h"

#ifndef BASE_PAIRS_H_
#define BASE_PAIRS_H_

int MIN_LOOP_SIZE; 
int MAX_LOOP_SIZE;
int MAX_BULGE_SIZE; 
int MIN_PERCENT_PAIRING;
int MAX_HELIX_SCOPE; /* default= size of the input sequence */
int MIN_HELIX_LENGTH; 
int MAX_DEGREE; /* default= size of the input sequence */

typedef struct {int final; int span; char thickness;} bp; 


bp * BP_list; /* list of all possible base pairs  */

void compute_all_base_pairs(plain_sequence *rna, int BPFILE, char * BP_file_name);
void free_base_pairs(plain_sequence * rna); 
void display_base_pairs(plain_sequence *rna);

int get_BP(int i, int j); 
int get_BP_index(int i, int j); 

#endif /*BASE_PAIRS_H*/
