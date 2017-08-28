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

#ifndef STACK_H_
#define STACK_H_


typedef struct {
  int start ; 
  int end;
  flat_cell * order;
  char thickness;
  int nb_of_bp;
} stack_cell;

stack_cell ** flat_stack;

int stack_height;
int stack_nb_of_bp;

void init_flat_stack(plain_sequence * rna);
int succ_flat_stack(int rna_size);

#endif
