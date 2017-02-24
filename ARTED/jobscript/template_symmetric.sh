#! /bin/bash
#SBATCH -J arted-symmetric
#SBATCH -p mixed
#SBATCH -N 2
#SBATCH --ntasks-per-node=2
#SBATCH -t 02:00:00
#SBATCH -o joblog/arted.symmetric.log

cd $SLURM_SUBMIT_DIR

module purge
module load intel intelmpi mkl

export OMP_STACKSIZE="2M"
export OMP_SCHEDULE="static"
export MIC_OMP_NUM_THREADS=180
export MIC_KMP_AFFINITY="scatter,granularity=fine"

# load-balancer
export ARTED_ENABLE_LOAD_BALANCER=1
export ARTED_CPU_PPN=${SLURM_NTASKS_PER_NODE}
export ARTED_MIC_PPN=2
export ARTED_CPU_TASK_RATIO=0.78

TARGET_MODE=sc

./script/mpirun.symm -s ./bin/ARTED_${TARGET_MODE} < ./data/input_${TARGET_MODE}_Si.dat

