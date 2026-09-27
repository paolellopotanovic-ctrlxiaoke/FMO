#!/bin/bash
#SBATCH -A XXXXXX
#SBATCH -J gamess-gpu-offload-tests
#SBATCH -o %x-%j.out
#SBATCH -t 2:00:00
#SBATCH --exclusive
#SBATCH --threads-per-core=2
#SBATCH --gpus-per-node=8
#SBATCH --nodes=1

## Load module environment
#

cd $SLURM_SUBMIT_DIR

srun $PWD/bin/my_ipcrm
OMP_NUM_THREADS=6 tests/runtest.py --folder=efmo-gpu-offload -n 8 --pre="./bin/my_ipcrm" --post="./bin/my_ipcrm"
OMP_NUM_THREADS=6 tests/runtest.py --folder=rimp2-gpu-offload -n 8 --pre="./bin/my_ipcrm" --post="./bin/my_ipcrm"  --skip_file=perlmutter
OMP_NUM_THREADS=6 ./rungms-dev ./tests/rimp2-gpu-offload/parallel/w150.moread.nprint.-5.inp 00 8 8 &> ./tests/rimp2-gpu-offload/parallel/w150.moread.nprint.-5.log
OMP_NUM_THREADS=6 ./rungms-dev ./tests/rimp2-gpu-offload/parallel/w175.moread.nprint.-5.inp 00 8 8 &> ./tests/rimp2-gpu-offload/parallel/w175.moread.nprint.-5.log
srun $PWD/bin/my_ipcrm

cd tests
./checkgms.py
./checkgms.py -v -p

