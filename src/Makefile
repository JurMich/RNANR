CC = gcc
#CFLAGS = -O3 -Wall
#-ansi removed from nextt line
CFLAGS =  -IViennaRNA -std=gnu99 -Wall -pedantic -pg -ggdb -g -O3 -fstack-protector-all 
LIBS =   -lRNA -lm -fopenmp

RNANR: main.o  rna.o  output_file.o base_pairs.o flat_structures.o locally_optimal_structures.o stack.o counting.o tree.o
	$(CC) $(CFLAGS) -o RNANR main.o rna.o output_file.o base_pairs.o flat_structures.o locally_optimal_structures.o stack.o counting.o tree.o $(LIBS)

main.o : main.c rna.h  base_pairs.h flat_structures.h flat_structures.h output_file.h

tree.o : tree.h tree.c

rna.o: rna.c rna.h 

output_file.o: output_file.h output_file.c

base_pairs.o: rna.h base_pairs.h base_pairs.c

flat_structures.o: rna.h base_pairs.h flat_structures.h flat_structures.c

stack.o: rna.h base_pairs.h flat_structures.h stack.h stack.c

counting.o : tree.h rna.h base_pairs.h flat_structures.h counting.h counting.c

locally_optimal_structures: rna.h base_pairs.h flat_structures.h stack.h locally_optimal_structures.h locally_optimal_structures.c
clean:
	rm -f *.o ; rm RNANR
