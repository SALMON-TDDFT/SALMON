#! /bin/bash
#PJM -L "rscgrp=fx-small"
#PJM -L "node=8"
#PJM -L "elapse=04:00:00"
#PJM --mpi "proc=8,rank-map-bychip"
#PJM -j
#PJM -o joblog/arted.fx100.log

TARGET_MODE=sc

mpiexec -stdin ./data/input_${TARGET_MODE}_Si.dat numactl -i all ./bin/ARTED_${TARGET_MODE}.cpu
