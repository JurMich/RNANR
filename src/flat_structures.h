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

#ifndef FLAT_STRUCTURES_H_
#define FLAT_STRUCTURES_H_

/* current is the index of the flat structure */
typedef struct flat_cell {int current; struct flat_cell * next; } flat_cell; 

/* Each flat structure is represented by a 3-uplet:  */
/* (start,end): positions of the leftmost base pair  */
/* suffix: index of the structure constituted by the other juxtaposed */
/* base pairs (this is a flat structure too) */

flat_cell *** flat_table; /* flat_table[x][y]: list of all flat structures on range x..y */  

/* mj_struct get_mj_list(int k); */

int get_flat_structure_start(int k);
int get_flat_structure_end(int k);
int get_flat_structure_suffix(int k);
int get_flat_structure_nb_of_bp(int k);

void display_one_flat_structure(int k); 
void display_all_flat_structures(); 

void build_all_flat_structures( plain_sequence * rna); 

void free_flat_list();
void free_flat_table(plain_sequence * rna);

#endif
