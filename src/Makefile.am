AM_CPPFLAGS = $(VRNA_CFLAGS)

bin_PROGRAMS = RNANR

LDADD = $(VRNA_LIBS)

RNANR_SOURCES = \
    main.c \
    base_pairs.c base_pairs.h \
    counting_double.c counting_double.h \
    energies.c energies.h \
    flat_structures.c flat_structures.h \
    locally_optimal_structures.c locally_optimal_structures.h \
    output_file.c output_file.h \
    rna.c rna.h \
    sort.c sort.h \
    stack.c stack.h \
    counting.inc \
    tree.inc

if RNANR_AM_SWITCH_MPFR
  RNANR_SOURCES += \
    counting_mpfr.c counting_mpfr.h

  LDADD += $(GMP_LIBS)
endif
