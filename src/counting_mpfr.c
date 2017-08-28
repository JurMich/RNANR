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

#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpfr.h> // MPFR
#include "rna.h"
#include "base_pairs.h"
#include "flat_structures.h"
#include "sort.h"
#include "stack.h"
#include "energies.h"
#include "counting_mpfr.h"

// basic definitions necessary for "conting.c" to function

/******************/
/* Define Type    */
/******************/
#define TYPE mpfr_t

TYPE ** partition_function_table; 

// parameters related with Vienna package
vrna_fold_compound_t *E_fold_cp;

/*** mpfr exclusive functions  and parameters ***/
gmp_randstate_t rng_state; //RNG

int precision(){ // returns precision
		return 200;
}


mpfr_rnd_t default_rnd(){  // returns default rounding mode (to nearest)
	return mpfr_get_default_rounding_mode(); 
}

/**************************************************/
/*       operators for the locopt structures      */
/**************************************************/

#define INIT(a) mpfr_init2(a, precision())   // associates memory, here does nothing
#define SET_TYPE_VAL(a, b) mpfr_set(a, b, default_rnd())   // sets value of a to b (which is TYPE)
#define SET_TYPE_VAL_FROM_DOUBLE(a, b) mpfr_set_d(a ,b, default_rnd()) //  sets value of a to b (which is double)
#define ADD(a,b) mpfr_add(a,a,b, default_rnd())       // equals to a += b
#define SUBSTRACT(a,b) mpfr_sub(a,a,b, default_rnd())
#define MULTIPLY(a,b) mpfr_mul(a,a,b, default_rnd())
#define DIVIDE(a,b) mpfr_div(a,a,b, default_rnd())
#define EXPONENT(a) mpfr_exp(a,a, default_rnd())
#define IS_LOWER(a, b) mpfr_less_p(a,b)
#define CLEAR(a) mpfr_clear(a)
#define TYPE_PRINT(a) mpfr_out_str(stdout, 10, 0, a, default_rnd()); // double only!  
#define DOUBLE_CAST(a) mpfr_get_d(a, default_rnd()) // casts TYPE to double. 
#define INIT_RNG() gmp_randinit_default (rng_state);  // initiates RNG
#define RANDOM_NUMBER(r, upper_bound) mpfr_urandomb(r, rng_state);\
                                      mpfr_mul(r,r, upper_bound, default_rnd())  

// redefinition of functions used outside included files

# define print_all_flat_structures(rna, file_path) print_all_flat_structures_mpfr(rna, file_path)
# define print_all_flat_structures_pile(rna, file_path) print_all_flat_structures_pile_mpfr(rna, file_path)
# define stochastic_backtrack_locally_optimal_structures(number_of_structures, rna, is_non_redun, use_timer, zqlimit, struc_count, DP_t) stochastic_backtrack_locally_optimal_structures_mpfr(number_of_structures, rna, is_non_redun, use_timer, zqlimit, struc_count, DP_t)
# define get_partition_function(partition_function, rna) get_partition_function_mpfr(partition_function, rna)
# define print_struct(rna) print_struct_mpfr(rna)

// renaming internal functions to avoid linker clash: tree
# define create_tr_node(id, parent) create_tr_node_m(id, parent)
# define create_ll_node(id) create_ll_node_m(id)
# define create_root() create_root_m()
# define insert_tr_node(id, ptr_prev, ptr, last_node) insert_tr_node_m(id, ptr_prev, ptr, last_node)
# define add_if_nexists(id, par_node) add_if_nexists_m(id, par_node)
# define tr_node_weight(returned_val, id, par_node) tr_node_weight_m(returned_val, id, par_node)
# define total_weight_par(returned_val, par_node) total_weight_par_m(returned_val, par_node)
# define update_weight(node, weight) update_weight_m(node, weight)
# define traceback_to_root(lead, struct_weight) traceback_to_root_m(lead, struct_weight)

// renaming other internal functions 
# define apply_energy_term(returned_value, dG, val) apply_energy_term_m(returned_value, dG, val)
# define rewire(x, y, linktable) rewire_m(x, y, linktable)
# define init_partition_table(n) init_partition_table_m(n) 
# define stochastic_backtrack_flat_structure_rec(current_flat_structure, rna, structure, partit_fci, actual_node, denomin) stochastic_backtrack_flat_structure_rec_m(current_flat_structure, rna, structure, partit_fci, actual_node, denomin) 
# define stochastic_backtrack_locally_optimal_structure_rec(x,y, rna, structure, partit_fci, actual_node, denomin) stochastic_backtrack_locally_optimal_structure_rec_m(x,y, rna, structure, partit_fci, actual_node, denomin) 


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