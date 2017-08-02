export OMP_NUM_THREADS=24
date > C2H2_rt_response.out
mpiexec -np 1 ../../../salmon.cpu < C2H2_rt_response.inp >> C2H2_rt_response.out
date >> C2H2_rt_response.out
