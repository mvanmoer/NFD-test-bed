PRG = pgm5

#FC = ifort 
FC = gfortran

#FCFLAGS = -free -openmp -O2 -o $@
#FCFLAGS = -ffree-form -fopenmp -O2 -o $@
FCFLAGS = -ffree-form -O2 -o $@

OBJECTS = initial_conditions.o diffusion.o \
	 boundary_conditions.o putfield.o \
	 pgm5.o 

# NOTE: might need dependencies to be .o instead of .mod depending
# on system/compiler.

all: $(PRG)

# default -- object files without dependencies
%.o: %.f90
	$(FC) -c $(FCFLAGS) $<

# .Fs
putfield.o: putfield.F
	$(FC) -c $(FCFLAGS) putfield.F

pgm5.o: pgm5.f90 diffusion.mod \
	initial_conditions.mod boundary_conditions.mod
	$(FC) -c $(FCFLAGS) $<

pgm5: $(OBJECTS)
	$(FC) $(FCFLAGS) $^

clean:
	rm -f $(OBJECTS) pgm5
	
deepclean:
	rm -f *.o pgm5 *~ *.exe *.stackdump *.mod gmeta 

.PHONY: all clean deepclean
