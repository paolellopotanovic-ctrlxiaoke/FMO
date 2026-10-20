cmake_minimum_required(VERSION 2.8)

project(NPY-for-Fortran Fortran)

add_library(npy_o OBJECT src/npy.F90)

if(DEFINED unit)
  message(STATUS "Unit is defined as ${unit}")
  # add_definitions(-D__UNIT__=${unit})
  target_compile_definitions(npy_o PUBLIC NPY_UNIT=${unit})
else()
  target_compile_definitions(npy_o PUBLIC NPY_UNIT=100)
endif()

add_executable(example $<TARGET_OBJECTS:npy_o> src/example.f90)

add_custom_command(
  TARGET example POST_BUILD
  COMMAND cp -f $<TARGET_OBJECTS:npy_o> npy.o
  BYPRODUCTS npy.o
  VERBATIM)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -ffree-form -ffree-line-length-0 -std=f2003")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stand f2003")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mfreeform -Mallocatable=03 -Mstandard")
endif()
