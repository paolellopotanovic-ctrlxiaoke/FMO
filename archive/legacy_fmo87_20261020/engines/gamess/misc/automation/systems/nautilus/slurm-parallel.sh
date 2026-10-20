#!/bin/bash

#SBATCH --exclusive
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --partition=standard
#SBATCH --job-name=gamess-parallel-testing
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

export PYTHONUNBUFFERED=1

cd $SLURM_SUBMIT_DIR

rm -rf restart/* scratch/*

tests/runtest.py \
--skip_folder=libcchem,mpiomp,omprimp2,openmp,cim,standard,benchmark,neb,neo,vb2000,aambs \
--skip_file=elmom,-omp,cc-cct3-be2,tpa-abn,cim-mp2,2lde-mmff94-ffdata \
--folder=parallel \
--ncpus=8 \
--mpi --skip_log \
--pre="./bin/my_ipcrm" --post="./bin/my_ipcrm"

cd tests

./checkgms.py -v -p --fail_rename=".log.PARALLEL.FAILED" --skip_json_create

cd $SLURM_SUBMIT_DIR

srun ./bin/my_ipcrm
