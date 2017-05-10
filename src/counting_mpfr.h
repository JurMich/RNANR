#ifndef COUNTING_MPFR_H_
#define COUNTING_MPFR_H_

#include <mpfr.h>

double TEMP; /*temperature */
double TEMPSCALE; /*temperature scaling factor*/
double RT; /*Boltzmann's cte*/

/* holds structures in bracket form and energies*/

void print_all_flat_structures_mpfr(plain_sequence * rna, char * file_path);
void print_all_flat_structures_pile_mpfr(plain_sequence * rna, char * file_path);
 
folding* stochastic_backtrack_locally_optimal_structures_mpfr(int number_of_structures, plain_sequence * rna,\
 int is_non_redun, int use_timer, float zqlimit, int * struc_count, double * DP_t);

void get_partition_function_mpfr(mpfr_t * partition_function, plain_sequence * rna);

#endif