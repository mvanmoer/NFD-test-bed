# NFD-test-bed
Finite difference numerical fluid dynamics code for testing accelerator technologies.

This is a 3D nonlinear quasi-compressible flow simulation. Originally written for ATMS 502 Numerical Fluid Dynamics, Fall 2011, UIUC, taught by Brian Jewett. This has since evolved into a program for me to learn new technologies such as: CMake, OpenMP, OpenACC, MPI, JSON/FSON

##TODO:
- ~~build on NCSA's Blue Waters with PGI compilers~~
- reapply OpenMP (initial commit is the serial version)
- apply hybrid MPI/OpenACC
- optional vtkXMLWriterF.h routine
- benchmarking/profiling

##Build instructions
Requirements 
- Fortran 2003 compiler (for streaming writes)
- CMake 2.6+
- fson: https://github.com/josephalevin/fson -- makefile will require editing for non-gfortran

1. Build out of source, cd to build dir.
2. ccmake ../path/to/src
3. set build type and paths to fson mods and object file.
4. make

run as $ ./nfd ../path/to/config.json

Two sample JSON config files are included, see the fortran source for explanation. highres.json will go unstable at approximately ~4700 steps.

Output is RAW binary, which can be loaded directly into ParaView or other visualization packages.

##Blue Waters status for serial:
- PrgEnv-gnu
  + produces expected output
- PrgEnv-intel
  + produces expected output
- PrgEnv-pgi 
  + produces expected ouput, but must add -Mnopre, in particular for advection.f90.
- PrgEnv-cray
  + produces expected output, but must use Makefile-cray instead of CMake for time being

## Timings
https://docs.google.com/spreadsheets/d/1e2N0prgzqzjtAuzFxk4Ui5Zxp7Yx1EaBiR8f9oWMEzs/edit?usp=sharing
