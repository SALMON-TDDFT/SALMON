export OMP_NUM_THREADS=24
date > output_rt_response
mpiexec -np 1 ../../../salmon.cpu < C2H2_rt_response.inp >> output_rt_response
date >> output_rt_response
