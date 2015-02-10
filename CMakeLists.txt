cmake_minimum_required(VERSION 2.8)

project (Serial Fortran)

set (NFD_SRCS advection.f90 boundary_conditions.f90 diffusion.f90 initial_conditions.f90 integration.f90 mesh.f90 pgf_bouyancy.f90 nfd.f90 savedata.f90)

set (FSON_build_PATH "NOT-DEFINED" CACHE PATH "docstring")
set (FSON_dist_PATH "NOT-DEFINED" CACHE PATH "docstring")
include_directories (${FSON_build_PATH})
link_directories (${FSON_dist_PATH})

set(EXTRA_LIBS ${EXTRA_LIBS} fson)

add_executable (nfd ${NFD_SRCS})
target_link_libraries(nfd ${EXTRA_LIBS})
