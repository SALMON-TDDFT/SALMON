export OMP_NUM_THREADS=24
date > C2H2_rt_pulse.out
mpiexec -np 1 ../../../salmon.cpu < C2H2_rt_pulse.inp >> C2H2_rt_pulse.out
date >> C2H2_rt_pulse.out
