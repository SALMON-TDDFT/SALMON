#! /bin/bash
#PBS -S /bin/bash
#PBS -N arted-gpu
#PBS -q tcag
#PBS -l select=1:ncpus=20:mpiprocs=4:ompthreads=5
#PBS -l place=scatter
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -o joblog/arted.gpu.log

cd ${PBS_O_WORKDIR}

module purge
module load pgi/16.4 mvapich2/2.1_pgi_cuda-7.5.18

TARGET_MODE=sc
NUM_MPI_PROCS=`wc -l $PBS_NODEFILE | awk '{print $1}'`

mpirun_rsh -np ${NUM_MPI_PROCS} -hostfile ${PBS_NODEFILE} \
  MV2_NUM_PORTS=2 MV2_ENABLE_AFFINITY=1 MV2_USE_CUDA=1 OMP_NUM_THREADS=5 MV2_SHOW_CPU_BINDING=1 MV2_CPU_MAPPING=0-4:5-9:10-14:15-19 \
  numactl --localalloc ./script/select_GPU ./bin/ARTED_${TARGET_MODE}.pgi_acc < ./data/input_${TARGET_MODE}.dat
