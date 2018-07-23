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

#ifndef RNA_H_
#define RNA_H_


#define MIN(a,b)      ((a)>(b) ? (b):(a))
#define MAXI(a,b)      ((a)>(b) ? (a):(b))

typedef struct{int size;char * label;}plain_sequence; 

void display_plain_sequence(plain_sequence * rna, FILE * outfile);
plain_sequence * get_plain_sequence(char * inputFile,char * RNAname);
int are_complementary(char a, char b); 
int are_Watson_Crick(char a, char b);
void free_plain_sequence(plain_sequence * rna);
#endif /*RNA_H_*/
