# Building GAMESS on Frontier

Note: You may need to `pip install jinja` if the `./bin/create-install-info.py` file complains about not finding the jinja2 module.

## OpenMP Offloaded GAMESS (GPU accelerated)

### CRAY

#### CPE 23.03/CCE 15.0.1 + ROCM 5.3.0

```
module reset
module load git
module load cmake
module load craype-accel-amd-gfx90a
module load PrgEnv-cray
module load cpe/23.03 cce/15.0.1 rocm/5.3.0
module unload darshan-runtime perftools-base
module use /lustre/orion/world-shared/chm135/2024/frontier/hipfort-cpe23.03-cce15.0.1-rocm5.3.0/modulefiles/cce/15.0/rocm5.3
module load hipfort
./bin/create-install-info.py --frontier --hipfort --hipfort_path=${OLCF_HIPFORT_ROOT}
make -j8
```

#### CPE 23.09/CCE 16.0.1 + ROCM 5.7.1

```
module reset
module load git
module load cmake
module load craype-accel-amd-gfx90a
module load PrgEnv-cray
module load cpe/23.09 cce/16.0.1 rocm/5.7.1
module unload darshan-runtime perftools-base
module use /lustre/orion/world-shared/chm135/2024/frontier/hipfort-cpe23.09-cce16.0.1-rocm5.7.1/modulefiles/cce/16.0/rocm5.7
module load hipfort
./bin/create-install-info.py --frontier --hipfort --hipfort_path=${OLCF_HIPFORT_ROOT}
make -j8

```

### OpenMP Offloaded Test Cases (as of March 17, 2024)

The following test cases have been verified to work on Perlmutter GPU nodes.

```
OMP_NUM_THREADS=6 tests/runtest.py --folder=efmo-gpu-offload -n 8 --pre="./bin/my_ipcrm" --post="./bin/my_ipcrm"
OMP_NUM_THREADS=6 tests/runtest.py --folder=rimp2-gpu-offload -n 8 --pre="./bin/my_ipcrm" --post="./bin/my_ipcrm"
OMP_NUM_THREADS=6 ./rungms-dev ./tests/rimp2-gpu-offload/parallel/w150.moread.nprint.-5.inp 00 8 8
OMP_NUM_THREADS=6 ./rungms-dev ./tests/rimp2-gpu-offload/parallel/w175.moread.nprint.-5.inp 00 8 8
```
