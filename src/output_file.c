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


