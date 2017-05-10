CC = gcc
#CFLAGS = -O3 -Wall
#-ansi removed from nextt line
CFLAGS =  -IViennaRNA -std=gnu99 -Wall -pedantic -pg -ggdb -g -O3 -fstack-protector-all 
LIBS =   -lRNA -lm -lmpfr -lgmp -fopenmp 

RNANR: main.o  rna.o  output_file.o base_pairs.o flat_structures.o locally_optimal_structures.o stack.o energies.o counting_double.o counting_mpfr.o sort.o
	$(CC) $(CFLAGS) -o RNANR main.o rna.o output_file.o base_pairs.o flat_structures.o locally_optimal_structures.o stack.o energies.o counting_double.o counting_mpfr.o sort.o $(LIBS)

main.o : main.c rna.h  base_pairs.h flat_structures.h energies.h counting_double.h counting_mpfr.h output_file.h

sort.o : flat_structures.h sort.h sort.c

rna.o: rna.c rna.h 

output_file.o: output_file.h output_file.c

base_pairs.o: rna.h base_pairs.h base_pairs.c

flat_structures.o: rna.h base_pairs.h flat_structures.h flat_structures.c

stack.o: rna.h base_pairs.h flat_structures.h stack.h stack.c

counting_double.o : sort.h rna.h base_pairs.h flat_structures.h energies.h counting_double.h counting_double.c tree.c counting.c

counting_mpfr.o : sort.h rna.h base_pairs.h flat_structures.h energies.h counting_mpfr.h counting_mpfr.c tree.c counting.c

energies.o : sort.h rna.h base_pairs.h flat_structures.h energies.h energies.c

locally_optimal_structures: rna.h base_pairs.h flat_structures.h stack.h locally_optimal_structures.h locally_optimal_structures.c
clean:
	rm -f *.o ; rm RNANR
