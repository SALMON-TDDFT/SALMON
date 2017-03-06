export OMP_NUM_THREADS=24
date > output_scf
mpiexec -np 1 ../../../salmon.cpu < input_scf >> output_scf
date >> output_scf
