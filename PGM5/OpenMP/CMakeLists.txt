cmake_minimum_required(VERSION 2.8)

project (Final Fortran)

set (PGM5_SRCS advection.f90 boundary_conditions.f90 diffusion.f90 initial_conditions.f90 integration.f90 mesh.f90 pgf_bouyancy.f90 pgm5.f90 savedata.f90)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")
FIND_PACKAGE(OpenMP_Fortran)
SET_TARGET_PROPERTIES(${EXECUTABLE_NAME} PROPERTIES
						COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
						LINK_FLAGS "${OpenMP_Fortran_FLAGS}")

include_directories("/home/mvanmoer/builds/fson/build")
link_directories("/home/mvanmoer/builds/fson/dist")
set(EXTRA_LIBS ${EXTRA_LIBS} fson)

add_executable (pgm5  ${PGM5_SRCS})
target_link_libraries(pgm5 ${EXTRA_LIBS})
