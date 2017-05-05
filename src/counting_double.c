#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/loop_energies.h>
#include "rna.h"
#include "base_pairs.h"
#include "flat_structures.h"
#include "sort.h"
#include "stack.h"
#include "counting_double.h"

// basic definitions necessary for "conting.c" to function

/******************/
/* Define Type    */
/******************/
#define TYPE double

long int ** number_of_locopt_structures;
TYPE ** partition_function_table; 

typedef struct {int flat_start; int flat_end;} flat_pair;

// parameters related with Vienna package
vrna_fold_compound_t *E_fold_cp;

/**************************************************/
/*       operators for the locopt structures      */
/**************************************************/

#define FOLDING folding_d
#define INIT(a)         // associates memory, here does nothing
#define SET_TYPE_VAL(a, b) a=b   // sets value of a to b (which is TYPE)
#define SET_TYPE_VAL_FROM_DOUBLE(a, b) a=b //  sets value of a to b (which is double)
#define ADD(a,b) a=a+b       // equals to a += b
#define SUBSTRACT(a,b) a=a-b
#define MULTIPLY(a,b) a=a*b
#define DIVIDE(a,b) a=a/b
#define EXPONENT(a) a=exp(a)
#define IS_NEGATIVE(a) (((a) < 0) ? (1) : (0))
#define CLEAR(a)
#define TYPE_PRINT(a) printf("%.3f", a) // double only!  
#define DOUBLE_CAST(a) a // casts TYPE to double. Here does nothing (becomes important with MPFR, see relevant .h file)
#define RANDOM_NUMBER(upper_bound) (double)rand()*upper_bound/((double)RAND_MAX + 1.)

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