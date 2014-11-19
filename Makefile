# This Makefile to be run from the parent directory of src/ and include/
# Use icc for Intel compiler

#----------------------------------------------------------------------------

#all: src/main.cpp src/gamma.cpp src/output.cpp src/phi.cpp src/stdafx.cpp src/utility.cpp src/operator.cpp src/phiC1.cpp
#	icc -Iinclude -mkl src/main.cpp src/gamma.cpp src/output.cpp src/phi.cpp src/stdafx.cpp src/utility.cpp src/operator.cpp src/phiC1.cpp
#	mv a.out bin/msm

#all: src/interpolant/b_spline/b_spline.c src/interpolant/c1_spline/c1_spline.c src/interpolant/interpolant.c src/memory/memory.c src/output/output.c src/softener/even_powers/even_powers.c src/softener/softener.c src/tester/main.c
#	icc -Wall -Wshadow -Iinclude -mkl  src/interpolant/b_spline/b_spline.c src/interpolant/c1_spline/c1_spline.c src/interpolant/interpolant.c src/memory/memory.c src/output/output.c src/softener/even_powers/even_powers.c src/softener/softener.c src/tester/main.c
#	mv a.out bin/msm

#----------------------------------------------------------------------------

IDIR	= include
SDIR	= src
LDIR	= lib
ODIR	= bin
BDIR	= bin

LIBS	= -mkl

CC		= icc
CFLAGS	= -Wall -Wshadow -I$(IDIR)

_DEPS	= all.h b_spline.h c1_spline.h even_powers.h interpolant.h memory.h method.h msm.h naive.h output.h particle.h particle_collection.h polynomial.h simulation.h simulation_domain.h simulator.h softener.h stencil.h tester.h
DEPS	= $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ	= b_spline.o c1_spline.o even_powers.o interpolant.o main.o memory.o method.o msm.o naive.o output.o particle_collection.o polynomial.o simulation.o simulation_domain.o simulator.o softener.o stencil.o tester.o
OBJ		= $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.c $(IDIR)/$(DEP)
	$(CC) -c -o $@ $< $(CFLAGS)

msm: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	mv $@ $(BDIR)/$@

.PHONY: clean

clean:
	rm -rf bin/*.o bin/msm

# End of file