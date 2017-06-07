export OMP_NUM_THREADS=24
date > output_rt_pulse
mpiexec -np 1 ../../../salmon.cpu < input_rt_pulse >> output_rt_pulse
date >> output_rt_pulse
