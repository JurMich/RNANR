#ifndef COUNTING_DOUBLE_H_
#define COUNTING_DOUBLE_H_

double TEMP; /*temperature */
double TEMPSCALE; /*temperature scaling factor*/
double RT; /*Boltzmann's cte*/

/* holds structures in bracket form and energies*/
//typedef struct folding {char** structures; TYPE* part_fcis; float* energy_ref;} folding as general form
typedef struct folding_d {char** structures; double* part_fcis; float* energy_ref;} folding_d; 

# define print_structure_double(structure, energy, energy_ref, rna) print_structure(structure, energy, energy_ref, rna)
void print_structure_double(char * structure, double energy, double energy_ref, plain_sequence * rna);

long int count_all_flat_structures(plain_sequence * rna);
long int count_all_locally_optimal_structures(plain_sequence * rna);

# define print_all_flat_structures_double(rna, file_path) print_all_flat_structures(rna, file_path)
void print_all_flat_structures_double(plain_sequence * rna, char * file_path);
# define print_all_flat_structures_pile_double(rna, file_path) print_all_flat_structures_pile(rna, file_path)
void print_all_flat_structures_pile_double(plain_sequence * rna, char * file_path);

# define stochastic_backtrack_locally_optimal_structures_double(number_of_structures, rna, is_non_redun, use_timer, zqlimit, struc_count, DP_t) stochastic_backtrack_locally_optimal_structures(number_of_structures, rna, is_non_redun, use_timer, zqlimit, struc_count, DP_t) 
folding_d* stochastic_backtrack_locally_optimal_structures_double(int number_of_structures, plain_sequence * rna,\
 int is_non_redun, int use_timer, float zqlimit, int * struc_count, double * DP_t);

# define get_partition_function_double(partition_function, rna) get_partition_function(partition_function, rna)
void get_partition_function_double(double * partition_function, plain_sequence * rna);

# define print_struct_double(rna) print_struct(rna)
void print_struct_double(plain_sequence * rna);

# define count_all_locally_optimal_structures2_double(rna) count_all_locally_optimal_structures2(rna)
void count_all_locally_optimal_structures2_double(plain_sequence * rna);

#endif