cmake_minimum_required(VERSION 3.9)

project(gamess)
enable_language(Fortran)

file(GLOB SOURCES
  "sources/*.f"
  "sources/*.F"
)

file(GLOB SOURCES_F90
  "sources/*.F90"
)

add_compile_options(-fdefault-integer-8 -fsyntax-only)
add_compile_options(-fallow-argument-mismatch)
add_compile_options(-Wall -Wextra -Wconversion-extra -Wuse-without-only)
add_compile_options(-Werror=uninitialized -Werror=align-commons)
set_source_files_properties(${SOURCES_F90} PROPERTIES COMPILE_FLAGS -ffree-line-length-none)

add_executable(gamess ${SOURCES} ${SOURCES_F90})