#!/bin/bash

#PBS -N polaris-test
#PBS -j oe
#PBS -l select=1:system=polaris:ncpus=64:mpiprocs=8:ompthreads=8
#PBS -l place=scatter:excl
#PBS -l walltime=3:00:0
#PBS -q preemptable
#PBS -A XXXX
#PBS -l filesystems=home:grand

mpiexec --np 1 --ppn 1 $PWD/bin/my_ipcrm

## Load your module environment
#

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
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/Br2.inp 00 4 4 &> tests/rimp2grd-mpiomp/Br2.log
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/I2.inp 00 4 4 &> tests/rimp2grd-mpiomp/I2.log
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/I2NCH.inp 00 4 4 &> tests/rimp2grd-mpiomp/I2NCH.log
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/benz.inp 00 4 4 &> tests/rimp2grd-mpiomp/benz.log
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/benz_optimize.inp 00 4 4 &> tests/rimp2grd-mpiomp/benz_optimize.log
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/cor.inp 00 4 4 &> tests/rimp2grd-mpiomp/cor.log
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/w2o2.hess.ecp.inp 00 4 4 &> tests/rimp2grd-mpiomp/w2o2.hess.log
OMP_NUM_THREADS=8 ./rungms-dev tests/rimp2grd-mpiomp/c60.inp 00 4 4 &> tests/rimp2grd-mpiomp/c60.log

cd tests
./checkgms.py
cd ..
mpiexec --np 1 --ppn 1 $PWD/bin/my_ipcrm
