#! /bin/bash
# generic execution file.
#
# Test environment
#   - GCC     4.4.7
#   - OpenMPI 1.10.3
#   - LAPACK  3.6.0
#
# Number of MPI process   = 2
# Number of OpenMP thread = 10
#

NUM_MPI_PROCS=2
export OMP_NUM_THREADS=10
export OMP_STACKSIZE="2M"
export OMP_SCHEDULE="static"

TARGET_MODE=sc
INPUT_FILE=./data/input_${TARGET_MODE}_Si.single.dat

mpirun \
  -np ${NUM_MPI_PROCS} \
  --map-by slot:PE=${OMP_NUM_THREADS} \
  --report-bindings \
  ./bin/ARTED_${TARGET_MODE}.cpu < ${INPUT_FILE}
