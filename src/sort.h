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

#ifndef SORT_H_
#define SORT_H_

/* defines and creates a data structure allowing to sort all flat structures
 * between the same x and y by their Boltzmann factor
 */

typedef struct flat_sort {double Boltz_partition; flat_cell * fcell;} flat_sort;
typedef struct linktab {flat_sort * links; int length;} linktab;

linktab * sort(linktab * linktable);
linktab * start_linktab();
void add_link_element(linktab * linktable, double Boltz_partition, flat_cell * fcell);
void free_linktab(linktab * linktable);

#endif