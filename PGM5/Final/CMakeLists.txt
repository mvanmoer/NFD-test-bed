cmake_minimum_required(VERSION 2.8)

project (Final Fortran)

set (PGM5_SRCS advection.f90 boundary_conditions.f90 diffusion.f90 initial_conditions.f90 integration.f90 mesh.f90 pgf_bouyancy.f90 pgm5.f90 savedata.f03)



link_directories("/home/mvanmoer/builds/fson/dist")
set(EXTRA_LIBS ${EXTRA_LIBS} fson-lib.so)

add_executable (pgm5  ${PGM5_SRCS})
set_target_properties(pgm5 PROPERTIES CMAKE_Fotran_MODULE_DIRECTORY "/home/mvanmoer/builds/fson/build")
target_link_libraries(pgm5 ${EXTRALIBS})
