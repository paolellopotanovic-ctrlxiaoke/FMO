# Building GAMESS on Perlmutter

Note: You may need to `pip install jinja` if the `./bin/create-install-info.py` file complains about not finding the jinja2 module.

## CPU threaded GAMESS

### GNU

#### GCC 10.3.0

```
module load PrgEnv-gnu/8.3.3
module load cpu/1.0
module load gcc/10.3.0
module unload perftools-base
./bin/create-install-info.py --perlmutter --openmp
make -j8
```

#### GCC 11.2.0 under PrgEnv-gnu/8.3.3

```
module load PrgEnv-gnu/8.3.3
module load cpu/1.0
module load gcc/11.2.0
module unload perftools-base
./bin/create-install-info.py --perlmutter --openmp
make -j8
```

#### GCC 11.2.0 under PrgEnv-gnu/8.4.0

```
module load PrgEnv-gnu/8.4.0
module load cpu/1.0
module load gcc/11.2.0
module unload perftools-base
./bin/create-install-info.py --perlmutter --openmp
make -j8
```

#### GCC 12.3.0

```
module load PrgEnv-gnu
module load cpu/1.0
module unload perftools-base
./bin/create-install-info.py --perlmutter --openmp
make -j8
```

### NVHPC

#### NVHPC 23.9

```
module load PrgEnv-nvhpc
module load cpu/1.0
module unload perftools-base
./bin/create-install-info.py --perlmutter --openmp
make -j8
```

### CPU OpenMP Test Cases (as of March 17, 2024)

The following test cases have been verified to work on Perlmutter CPU nodes.

```
OMP_NUM_THREADS=4  ./rungms-dev tests/rimp2grd-mpiomp/Br2.inp 00 4 4
OMP_NUM_THREADS=4  ./rungms-dev tests/rimp2grd-mpiomp/I2.inp 00 4 4
OMP_NUM_THREADS=4  ./rungms-dev tests/rimp2grd-mpiomp/I2NCH.inp 00 4 4
OMP_NUM_THREADS=32 ./rungms-dev tests/rimp2grd-mpiomp/benz_optimize.inp 00 4 4
OMP_NUM_THREADS=16 ./rungms-dev tests/rimp2grd-mpiomp/benz_optimize.inp 00 8 8
OMP_NUM_THREADS=32 ./rungms-dev tests/rimp2grd-mpiomp/benz.inp 00 4 4
OMP_NUM_THREADS=16 ./rungms-dev tests/rimp2grd-mpiomp/benz.inp 00 8 8
OMP_NUM_THREADS=32 ./rungms-dev tests/rimp2grd-mpiomp/cor.inp 00 4 4
OMP_NUM_THREADS=16 ./rungms-dev tests/rimp2grd-mpiomp/cor.inp 00 8 8
OMP_NUM_THREADS=32 ./rungms-dev tests/rimp2grd-mpiomp/c60.inp 00 4 4
OMP_NUM_THREADS=16 ./rungms-dev tests/rimp2grd-mpiomp/c60.inp 00 8 8
OMP_NUM_THREADS=8  ./rungms-dev tests/rimp2grd-mpiomp/c60.inp 00 16 16
OMP_NUM_THREADS=32 ./rungms-dev tests/rimp2grd-mpiomp/w2o2.hess.ecp.inp 00 4 4
OMP_NUM_THREADS=16 ./rungms-dev tests/rimp2grd-mpiomp/w2o2.hess.ecp.inp 00 8 8
```

## OpenMP Offloaded GAMESS (GPU accelerated)

### NVHPC

#### NVHPC 22.7

```
module swap cudatoolkit/12.2 cudatoolkit/11.7
module load PrgEnv-nvhpc/8.4.0
module load nvhpc/22.7
module load cray-libsci/23.05.1.4
module load cray-mpich/8.1.27
./bin/create-install-info.py --perlmutter --openmp --openmp-offload --cublas
make -j8

```

#### NVHPC 23.1

```
module load PrgEnv-nvhpc/8.4.0
module swap nvhpc/23.9 nvhpc/23.1
module swap cray-libsci/23.12.5 cray-libsci/23.09.1.1
module unload perftools-base
./bin/create-install-info.py --perlmutter --openmp --openmp-offload --cublas
make -j8

```

#### NVHPC 23.9

```
module load PrgEnv-nvhpc
module unload perftools-base
./bin/create-install-info.py --perlmutter --openmp --openmp-offload --cublas
make -j8
```

### OpenMP Offloaded Test Cases (as of March 17, 2024)

The following test cases have been verified to work on Perlmutter GPU nodes.

```
OMP_NUM_THREADS=16 ./tests/runtest.py --folder=efmo-gpu-offload -n 4 --pre="./bin/my_ipcrm" --post="./bin/my_ipcrm"

OMP_NUM_THREADS=16 ./rungms-dev tests/rimp2-gpu-offload/parallel/benz.inp 00 4 4
OMP_NUM_THREADS=16 ./rungms-dev tests/rimp2-gpu-offload/parallel/c60.inp 00 4 4
OMP_NUM_THREADS=16 ./rungms-dev tests/rimp2-gpu-offload/parallel/cor.inp 00 4 4

```

