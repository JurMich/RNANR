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

#include "rna.h"

#ifndef BASE_PAIRS_H_
#define BASE_PAIRS_H_

int MIN_LOOP_SIZE; 
int MAX_LOOP_SIZE;
int MAX_BULGE_SIZE; 
int MIN_PERCENT_PAIRING;
int MAX_HELIX_SCOPE; /* default= size of the input sequence */
int MIN_HELIX_LENGTH; 
int MAX_DEGREE; /* default= size of the input sequence */

typedef struct {int final; int span; char thickness;} bp; 


bp * BP_list; /* list of all possible base pairs  */

void compute_all_base_pairs(plain_sequence *rna, int BPFILE, char * BP_file_name);
void free_base_pairs(plain_sequence * rna); 
void display_base_pairs(plain_sequence *rna);

int get_BP(int i, int j); 
int get_BP_index(int i, int j); 

#endif /*BASE_PAIRS_H*/
