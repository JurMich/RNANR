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

#include <stdio.h>
#include "output_file.h"

FILE *  open_outputfile(int FILEOUT, char * name){

  FILE * outfile; 

  if (FILEOUT==1)
    outfile=fopen(name, "w") ;
  else 
    outfile=stdout; 
  return outfile; 
}

void close_outputfile(int FILEOUT, FILE * outfile){
  if (FILEOUT==1)
    fclose(outfile); 
}


