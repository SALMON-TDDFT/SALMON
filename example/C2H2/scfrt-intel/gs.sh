export OMP_NUM_THREADS=24
date > output_gs
mpiexec -np 1 ../../../salmon.cpu < C2H2_gs.inp >> output_gs
date >> output_gs
