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
