export OMP_NUM_THREADS=24
date > C2H2_gs.out
mpiexec -np 1 ../../../salmon.cpu < C2H2_gs.inp >> C2H2_gs.out
date >> C2H2_gs.out
