export OMP_NUM_THREADS=24
date > output_rt_pulse
mpiexec -np 1 ../../../salmon.cpu < C2H2_rt_pulse.inp >> output_rt_pulse
date >> output_rt_pulse
