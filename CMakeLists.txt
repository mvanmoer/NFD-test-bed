cmake_minimum_required(VERSION 2.8)

project (Serial Fortran)

set (PGM5_SRCS advection.f90 boundary_conditions.f90 diffusion.f90 initial_conditions.f90 integration.f90 mesh.f90 pgf_bouyancy.f90 pgm5.f90 savedata.f90)

set (FSON_LIB_PATH "NOT-DEFINED" CACHE PATH "docstring")
set (FSON_INC_PATH "NOT-DEFINED" CACHE PATH "docstring")
include_directories (${FSON_INC_PATH})
link_directories (${FSON_LIB_PATH})

set(EXTRA_LIBS ${EXTRA_LIBS} fson)

add_executable (pgm5  ${PGM5_SRCS})
target_link_libraries(pgm5 ${EXTRA_LIBS})
