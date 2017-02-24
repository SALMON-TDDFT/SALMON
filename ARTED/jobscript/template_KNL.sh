#! /bin/bash
#PJM -L rscgrp=FLAT-QUADRANT,node=4
#PJM --mpi proc=4
#PJM --omp thread=256
#PJM -o joblog/arted.knl.log
#PJM -j

export I_MPI_PIN_DOMAIN=272

export KMP_AFFINITY=balanced,granularity=fine
export KMP_HW_SUBSET=64c,4t,2O
export OMP_NUM_THREADS=256

TARGET_MODE=sc

mpiexec.hydra -n ${PJM_MPI_PROC} numactl -m 1 ./bin/ARTED_${TARGET_MODE}.mic < ./data/input_${TARGET_MODE}_Si.dat
