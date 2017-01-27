#ifndef RNA_H_
#define RNA_H_


#define MIN(a,b)      ((a)>(b) ? (b):(a))
#define MAXI(a,b)      ((a)>(b) ? (a):(b))

typedef struct{int size;char * label;}plain_sequence; 

void display_plain_sequence(plain_sequence * rna);
plain_sequence * get_plain_sequence(char * inputFile,char * RNAname);
int are_complementary(char a, char b); 
int are_Watson_Crick(char a, char b);
void free_plain_sequence(plain_sequence * rna);
#endif /*RNA_H_*/
