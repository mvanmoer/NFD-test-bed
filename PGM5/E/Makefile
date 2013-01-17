PRG = pgm5

#FC = ifort 
FC = gfortran

#FCFLAGS = -free -openmp -O2 -o $@
#FCFLAGS = -ffree-form -fopenmp -O2 -o $@
FCFLAGS = -ffree-form -O2 -o $@

OBJECTS = advection.o initial_conditions.o integration.o \
	 boundary_conditions.o putfield.o diffusion.o \
	 pgf_bouyancy.o pgm5.o 

# NOTE: might need dependencies to be .o instead of .mod depending
# on system/compiler.

all: $(PRG)

# default -- object files without dependencies
%.o: %.f90
	$(FC) -c $(FCFLAGS) $<

# .Fs
putfield.o: putfield.F
	$(FC) -c $(FCFLAGS) putfield.F

# .os that depend on modules
integration.o: integration.f90 advection.mod
	$(FC) -c $(FCFLAGS) $<

pgm5.o: pgm5.f90 integration.mod diffusion.mod \
	initial_conditions.mod boundary_conditions.mod \
	pgf_bouyancy.mod
	$(FC) -c $(FCFLAGS) $<

pgm5: $(OBJECTS)
	$(FC) $(FCFLAGS) $^

clean:
	rm -f $(OBJECTS) pgm5
	
deepclean:
	rm -f *.o pgm5 *~ *.exe *.stackdump *.mod gmeta 

.PHONY: all clean deepclean
