# Building GAMESS on Polaris

Note: You may need to `pip install jinja` if the `./bin/create-install-info.py` file complains about not finding the jinja2 module.

```

## OpenMP Offloaded GAMESS (GPU accelerated)

### NVHPC

#### NVHPC 23.3

```
module swap nvhpc nvhpc/23.3
module load cray-libsci
module unload perftools-base
./bin/create-install-info.py --polaris --openmp-offload --cublas
make -j8
```

### OpenMP Offloaded Test Cases (as of March 17, 2024)

The following test cases have been verified to work on Polaris

```
OMP_NUM_THREADS=1 ./tests/runtest.py --supress -n 4 --filepath=efmo-gpu-offload/parallel/efmo_gafmo_gpu_FF.inp
OMP_NUM_THREADS=1 ./tests/runtest.py --supress -n 4 --filepath=efmo-gpu-offload/parallel/efmo_gafmo_gpu_FT.inp
OMP_NUM_THREADS=1 ./tests/runtest.py --supress -n 4 --filepath=efmo-gpu-offload/parallel/efmo_gafmo_gpu_TF.inp
OMP_NUM_THREADS=1 ./tests/runtest.py --supress -n 4 --filepath=efmo-gpu-offload/parallel/efmo_gafmo_gpu_TT.inp
OMP_NUM_THREADS=1 ./tests/runtest.py --supress -n 4 --filepath=efmo-gpu-offload/parallel/h2o_makefp_FT.inp
OMP_NUM_THREADS=1 ./tests/runtest.py --supress -n 4 --filepath=efmo-gpu-offload/parallel/h2o_makefp_TF.inp
OMP_NUM_THREADS=1 ./tests/runtest.py --supress -n 4 --filepath=efmo-gpu-offload/parallel/h2o_makefp_TT.inp
OMP_NUM_THREADS=1 ./tests/runtest.py --supress -n 4 --filepath=efmo-gpu-offload/parallel/h2o_rhf.inp
OMP_NUM_THREADS=1 ./tests/runtest.py --supress -n 4 --filepath=efmo-gpu-offload/parallel/meoh_makefp_FT.inp
OMP_NUM_THREADS=1 ./tests/runtest.py --supress -n 4 --filepath=efmo-gpu-offload/parallel/meoh_makefp_TF.inp
OMP_NUM_THREADS=1 ./tests/runtest.py --supress -n 4 --filepath=efmo-gpu-offload/parallel/meoh_makefp_TT.inp
OMP_NUM_THREADS=1 ./tests/runtest.py --supress -n 4 --filepath=efmo-gpu-offload/parallel/meoh_rhf.inp
OMP_NUM_THREADS=8 ./tests/runtest.py --supress -n 4 --filepath=rimp2-gpu-offload/parallel/benz.inp
OMP_NUM_THREADS=8 ./tests/runtest.py --supress -n 4 --filepath=rimp2-gpu-offload/parallel/cor.inp
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/Br2.inp 00 4 4
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/I2.inp 00 4 4
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/I2NCH.inp 00 4 4
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/benz.inp 00 4 4
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/benz_optimize.inp 00 4 4
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/cor.inp 00 4 4
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/w2o2.hess.ecp.inp 00 4 4
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/c60.inp 00 4 4

```

