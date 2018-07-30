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

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rna.h"


void wrong_format(char c){
      printf("WRONG FORMAT: Found %c.\n", c);
      exit(EXIT_SUCCESS);
}

char RNAcharacter(char c, char * FILENAME){ 
  c=toupper(c);
  if (c=='T') c='U';
  if (!((c=='A')||(c=='C')||(c=='G')||(c=='U')||(c=='N')))
    wrong_format(c);
  return c;
}

int MAX(int a, int b){
  if (a<b) 
    return b; 
  else 
    return a; 
}

int MINI(int a, int b){
  if (a<b) 
    return a; 
  else 
    return b; 
}


plain_sequence * get_plain_sequence(char * inputFile,char * RNAname){
  FILE *fp;
  plain_sequence * rna;
  int size=0;
  int i, skipped;
  char ch;
  char name[200]; 
  fp = fopen(inputFile, "rt");
  if(fp==NULL){
    printf("\nFile not found : %s\n",inputFile);
    exit(EXIT_FAILURE);
  }
  fgets(name,200,fp);
  i=1; 
  while ((name[i]!='\n')&&(name[i]!=NULL)) i++;
  size=MINI(199,i); 
  for (i=0; i<size-1; i++){
    RNAname[i]=name[i+1]; 
  }
  RNAname[size]='\0';
  size=0;
  while (fgetc(fp)!=EOF)  size++; 
  fclose(fp);
  if (size>100000) size=100000; 
  rna=(plain_sequence*) malloc(sizeof(plain_sequence)); 
  rna->label=(char*) malloc(sizeof(char)*(size+2));
  rna->size=size;
  fp = fopen(inputFile, "rt");
  fgets(name,200,fp);
  skipped = 0;
  for (i=1; i<=size; i++){
    fscanf(fp, "%c",&ch);
    if(ch != '\n'){
      rna->label[i-skipped]=RNAcharacter(ch,inputFile);  
    }else{
	  skipped += 1;	
	} 
  }
  fclose(fp);
  rna->size = size-skipped;
  rna->label[size+1-skipped]='\0';
  return rna;
}


/* display a plain sequence */
void display_plain_sequence(plain_sequence * rna, FILE * outfile){
 int i;
 if(!outfile){
	outfile = stdout;	
 }
 for(i=1;i<=rna->size;i++)
   fprintf(outfile, "%c",rna->label[i]);
 fprintf(outfile, "\n");
}

/* free a plain sequence */
void free_plain_sequence( plain_sequence * rna){
  free(rna->label);
  free(rna);
}


int are_Watson_Crick(char a, char b){
  switch (a){
  case 'A' : return (b=='U');
    break;
  case 'C' : return (b=='G');
    break;
  case 'U' : return (b=='A');
    break;
  case 'G' : return (b=='C');
    break;
  default : return 0;
    break;
 }
}

/* check whether two nucleotides can be paired: */
/* Watson-Crick (A-U, C-G)  or wobble G-U */
int are_complementary(char a, char b){
  switch (a){
  case 'A' : return (b=='U');
    break;
  case 'C' : return (b=='G');
    break;
  case 'U' : return ((b=='A')||(b=='G'));
    break;
  case 'G' : return ((b=='C')||(b=='U'));
    break;
  default : return 0;
    break;
 }
}
