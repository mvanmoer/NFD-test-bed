PRG = pgm5

#FC = ifort 
FC = gfortran

#FCFLAGS = -free -openmp -O2 -o $@
#FCFLAGS = -ffree-form -fopenmp -O2 -o $@
FCFLAGS = -ffree-form -O3 -o $@

OBJECTS = advection.o initial_conditions.o integration.o \
	 boundary_conditions.o savedata.o diffusion.o \
	 pgf_bouyancy.o readinputparameters.o mesh.o pgm5.o 

# NOTE: might need dependencies to be .o instead of .mod depending
# on system/compiler.

all: $(PRG)

# default -- object files without dependencies
%.o: %.f90
	$(FC) -c $(FCFLAGS) $<

# .Fs
savedata.o: savedata.f03
	$(FC) -c $(FCFLAGS) savedata.f03

# .os that depend on modules
# need to add proper handling of mesh_type.mod, which is in quite a few files.
integration.o: integration.f90 advection.mod mesh_type.mod
	$(FC) -c $(FCFLAGS) $<

pgm5.o: pgm5.f90 integration.mod diffusion.mod \
	initial_conditions.mod boundary_conditions.mod \
	pgf_bouyancy.mod savedata.mod read_parameters_module.mod mesh_type.mod
	$(FC) -c $(FCFLAGS) $<

pgm5: $(OBJECTS)
	$(FC) $(FCFLAGS) $^

clean:
	rm -f $(OBJECTS) pgm5
	
deepclean:
	rm -f *.o pgm5 *~ *.exe *.stackdump *.mod gmeta 

.PHONY: all clean deepclean
