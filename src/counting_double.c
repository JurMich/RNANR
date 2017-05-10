#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "rna.h"
#include "base_pairs.h"
#include "flat_structures.h"
#include "sort.h"
#include "stack.h"
#include "energies.h"
#include "counting_double.h"

// basic definitions necessary for "conting.c" to function

/******************/
/* Define Type    */
/******************/
#define TYPE double

TYPE ** partition_function_table; 


/**************************************************/
/*       operators for the locopt structures      */
/**************************************************/

#define INIT(a)         // associates memory, here does nothing
#define SET_TYPE_VAL(a, b) a=b   // sets value of a to b (which is TYPE)
#define SET_TYPE_VAL_FROM_DOUBLE(a, b) a=b //  sets value of a to b (which is double)
#define ADD(a,b) a=a+b       // equals to a += b
#define SUBSTRACT(a,b) a=a-b
#define MULTIPLY(a,b) a=a*b
#define DIVIDE(a,b) a=a/b
#define EXPONENT(a) a=exp(a)
#define IS_LOWER(a, b) (((a) < (b)) ? (1) : (0))
#define CLEAR(a)
#define TYPE_PRINT(a) printf("%.3f", a) // double only!  
#define DOUBLE_CAST(a) a // casts TYPE to double. Here does nothing (becomes important with MPFR, see relevant .h file)
#define INIT_RNG() srand(time(NULL)); // initiates RNG.
#define RANDOM_NUMBER(r, upper_bound) r = (double)rand()*upper_bound/((double)RAND_MAX + 1.)

// redefinition of fcis used in main

# define print_all_flat_structures(rna, file_path) print_all_flat_structures_double(rna, file_path)
# define print_all_flat_structures_pile(rna, file_path) print_all_flat_structures_pile_double(rna, file_path)
# define stochastic_backtrack_locally_optimal_structures(number_of_structures, rna, is_non_redun, use_timer, zqlimit, struc_count, DP_t) stochastic_backtrack_locally_optimal_structures_double(number_of_structures, rna, is_non_redun, use_timer, zqlimit, struc_count, DP_t)
# define get_partition_function(partition_function, rna) get_partition_function_double(partition_function, rna)
# define print_struct(rna) print_struct_double(rna)

// renaming internal functions to avoid linker clash: tree
# define create_tr_node(id, parent) create_tr_node_d(id, parent)
# define create_ll_node(id) create_ll_node_d(id)
# define create_root() create_root_d()
# define insert_tr_node(id, ptr_prev, ptr, last_node) insert_tr_node_d(id, ptr_prev, ptr, last_node)
# define add_if_nexists(id, par_node) add_if_nexists_d(id, par_node)
# define tr_node_weight(returned_val, id, par_node) tr_node_weight_d(returned_val, id, par_node)
# define total_weight_par(returned_val, par_node) total_weight_par_d(returned_val, par_node)
# define update_weight(node, weight) update_weight_d(node, weight)
# define traceback_to_root(lead, struct_weight) traceback_to_root_d(lead, struct_weight)

// renaming other internal functions 
# define apply_energy_term(returned_value, dG, val) apply_energy_term_d(returned_value, dG, val)
# define rewire(x, y, linktable) rewire_d(x, y, linktable)
# define init_partition_table(n) init_partition_table_d(n) 
# define stochastic_backtrack_flat_structure_rec(current_flat_structure, rna, structure, partit_fci, actual_node, denomin) stochastic_backtrack_flat_structure_rec_d(current_flat_structure, rna, structure, partit_fci, actual_node, denomin) 
# define stochastic_backtrack_locally_optimal_structure_rec(x,y, rna, structure, partit_fci, actual_node, denomin) stochastic_backtrack_locally_optimal_structure_rec_d(x,y, rna, structure, partit_fci, actual_node, denomin) 
  
// inclusion of other files   

/**************************************************/
/*                     tree                       */
/**************************************************/

// data structure definition
#include "tree.c"

/**************************************************/
/*                   counting                     */
/**************************************************/

// functions to 
#include "counting.c"