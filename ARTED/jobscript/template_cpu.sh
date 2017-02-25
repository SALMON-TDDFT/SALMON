#! /bin/bash
#SBATCH -J arted-cpu
#SBATCH -p mixed
#SBATCH -N 4
#SBATCH --ntasks-per-node=2
#SBATCH -t 02:00:00
#SBATCH -o joblog/arted.cpu.log

cd $SLURM_SUBMIT_DIR

module purge
module load intel intelmpi mkl

export OMP_STACKSIZE="2M"
export OMP_SCHEDULE="static"

TARGET_MODE=sc

./script/mpirun.symm -c ./bin/ARTED_${TARGET_MODE} < ./data/input_${TARGET_MODE}_Si.dat

