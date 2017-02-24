#! /bin/bash
#SBATCH -J arted-native
#SBATCH -p mixed
#SBATCH -N 4
#SBATCH -t 02:00:00
#SBATCH -o joblog/arted.native.log

cd $SLURM_SUBMIT_DIR

module purge
module load intel intelmpi mkl

export OMP_STACKSIZE="2M"
export OMP_SCHEDULE="static"
export MIC_KMP_AFFINITY="scatter,granularity=fine"

TARGET_MODE=sc

./script/mpirun.symm -m ./bin/ARTED_${TARGET_MODE} < ./data/input_${TARGET_MODE}_Si.dat

