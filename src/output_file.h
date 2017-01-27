#include <stdio.h>

#ifndef OUPTUT_H_
#define OUTPUT_H_

FILE *  open_outputfile(int FILEOUT, char * name); 

void close_outputfile(int FILEOUT, FILE * outfile); 


#endif /*OUTPUT_H_*/
