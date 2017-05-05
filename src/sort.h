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