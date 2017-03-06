export OMP_NUM_THREADS=24
date > output_rt
mpiexec -np 1 ../../../salmon.cpu < input_rt >> output_rt
date >> output_rt
