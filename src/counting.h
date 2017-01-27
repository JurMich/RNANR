#ifndef COUNTING_H_
#define COUNTING_H_

double TEMP; /*temperature */
double TEMPSCALE; /*temperature scaling factor*/
double RT; /*Boltzmann's cte*/

/* holds structures in bracket form and energies*/
typedef struct folding {char** structures; TYPE* energies; TYPE* energy_ref;} folding; 

TYPE random_number(TYPE upper_bound);
void print_structure(char * structure, TYPE energy, TYPE energy_ref, plain_sequence * rna);

int count_all_flat_structures(plain_sequence * rna);
TYPE count_all_locally_optimal_structures(plain_sequence * rna);

void print_all_flat_structures(plain_sequence * rna, char * file_path);
void print_all_flat_structures_pile(plain_sequence * rna, char * file_path);

folding* stochastic_backtrack_locally_optimal_structures(int number_of_structures, plain_sequence * rna,\
 int is_non_redun, int use_timer, float zqlimit, int * struc_count, double * DP_t);

TYPE get_partition_function(plain_sequence * rna);

void print_struct(plain_sequence * rna);

void count_all_locally_optimal_structures2(plain_sequence * rna);

#endif
