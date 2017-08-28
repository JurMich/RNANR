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

#ifndef COUNTING_DOUBLE_H_
#define COUNTING_DOUBLE_H_

/* holds structures in bracket form and energies*/

void print_all_flat_structures_double(plain_sequence * rna, char * file_path);
void print_all_flat_structures_pile_double(plain_sequence * rna, char * file_path);
 
folding* stochastic_backtrack_locally_optimal_structures_double(int number_of_structures, plain_sequence * rna,\
 int is_non_redun, int use_timer, float zqlimit, int * struc_count, double * DP_t);

void get_partition_function_double(double * partition_function, plain_sequence * rna);

void print_struct_double(plain_sequence * rna);

void count_all_locally_optimal_structures2(plain_sequence * rna);

#endif