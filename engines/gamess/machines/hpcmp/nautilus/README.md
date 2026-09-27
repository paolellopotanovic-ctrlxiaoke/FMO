# Nautilus

## Prerequisites

### Install Jinja 2

`pip3 install --user jinja2`

### Install your own OpenBLAS

```
wget https://github.com/xianyi/OpenBLAS/releases/download/v0.3.23/OpenBLAS-0.3.23.tar.gz
tar -xf OpenBLAS-0.3.23.tar.gz
cd OpenBLAS-0.3.23
make CC=gcc USE_THREAD=0 USE_LOCKING=1 INTERFACE64=1
make install PREFIX=$HOME/opt/openblas-0.3.23-gnu-8.5.0-serial
```

## Build combinations

### GNU + OpenBLAS + MPICH

```
module load mpich/gnu/3.3.2

./bin/create-install-info.py --nautilus --gfortran --fortran_version=8.5.0 \
--openblas --mathlib_path=$HOME/opt/openblas-0.3.23-gnu-8.5.0-serial/lib \
--mpich --mpi_path=/opt/scyld/mpich/3.3.2/gnu

make -j8
```

### Intel OneAPI + Intel MKL + Intel MPI

```
module load intel/compiler/latest intel/mkl/latest intel/mpi/latest

./bin/create-install-info.py --nautilus --oneapi-ifort --fortran_version=2022.1 \
--mkl \
--impi

make -j8
```

### AMD AOCC + AMD AOCL + MPICH

```
module load amd/aocc/4.0.0 amd/aocl/aocc/4.0 mpich/gnu/3.3.2

./bin/create-install-info.py --nautilus --aocc --fortran_version=4.0.0 \
--libflame --mathlib_path=$AOCL_ROOT/lib \
--mpich --mpi_path=/opt/scyld/mpich/3.3.2/gnu

make -j8
```

### AMD AOCC + AMD AOCL + OPENMPI

```
module load amd/aocc/4.0.0 amd/aocl/aocc/4.0 penguin/openmpi/4.1.4/aocc

./bin/create-install-info.py --nautilus --aocc --fortran_version=4.0.0 \
--libflame --mathlib_path=$AOCL_ROOT/lib \
--openmpi --mpi_path=$MPI_HOME

make -j8
```

### NVHPC + NVBLAS + MPICH

```
module load nvidia/nvhpc-nompi/22.3 mpich/gnu/3.3.2

./bin/create-install-info.py --nautilus --nvfortran \
--nvblas --mathlib_path=$NVHPC_ROOT/math_libs/lib64 \
--mpich --mpi_path=/opt/scyld/mpich/3.3.2/gnu

make -j8
```

### NVHPC + NVBLAS + OPENMPI

```
module load nvidia/nvhpc-nompi/22.3 penguin/openmpi/4.1.4/gcc

./bin/create-install-info.py --nautilus --nvfortran \
--nvblas --mathlib_path=$NVHPC_ROOT/math_libs/lib64 \
--openmpi --mpi_path=$MPI_HOME

make -j8
```
