#! /bin/bash
#PJM --rsc-list "node=8"
#PJM --rsc-list "elapse=04:00:00"
#PJM --mpi "rank-map-bynode"
#PJM -j
#PJM --stgin "./bin/ARTED_sc.cpu ./ARTED.cpu"
#PJM --stgin-dir  "./data ./data"
#PJM --stgout-dir "./data ./data"
#PJM -o joblog/arted.K.log

. /work/system/Env_base

TARGET_MODE=sc

mpiexec -stdin ./data/input_${TARGET_MODE}_Si.dat ./ARTED_${TARGET_MODE}.cpu
