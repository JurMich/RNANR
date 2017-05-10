#ifndef ENERGIES_H_
#define ENERGIES_H_

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>

double TEMP; /*temperature */
double TEMPSCALE; /*temperature scaling factor*/
double RT; /*Boltzmann's cte*/

long int ** number_of_locopt_structures;

// parameters related with Vienna package
vrna_fold_compound_t *E_fold_cp;

typedef struct {int flat_start; int flat_end;} flat_pair;
typedef struct folding {char** structures; double* part_fcis; float* energy_ref;} folding; 

// utility functions
void print_structure(char * structure, double part_fci, float energy_ref, plain_sequence * rna);
void add_base_pair(int i, int j, char * structure);
void add_helix(int i, int j, int thickness, char * structure);
int base2int(char base);
int get_type(char base_5, char base_3);
int is_entirely_contained(int i, int j, plain_sequence * rna);

// energy computation functions
double stacking_energy(int i, int j);
double helix_energy(int i, int j, int length);
double multiloop_energy(int alpha, int beta, int *ptypes, int *si, int *sj);
double internal_loop_energy(int i, int k, int l, int j, plain_sequence *rna);
double hairpin_energy(int i, int j, plain_sequence * rna);
double ext_loop_energy(int i, int j, plain_sequence *rna);

// other functions used in counting.h
long int count_all_locally_optimal_structures(plain_sequence * rna);
long int count_all_flat_structures(plain_sequence * rna);
void print_all_flat_structures_pile(plain_sequence * rna, char * file_path);
float get_reference_energy(vrna_fold_compound_t *E_fold_cp, char * structure, plain_sequence * rna);
void count_all_locally_optimal_structures2(plain_sequence * rna);

#endif