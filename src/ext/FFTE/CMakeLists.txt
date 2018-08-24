set(SOURCES
    alltoall_1d.f90
    factor.f
    fft235.f
    kernel.f
    pzfft3dv_mod.f
   )

add_library(${EXT_FFTE_LIB} ${SOURCES})
add_dependencies(${SALMON_IO_LIB} ${SALMON_PARALLEL_LIB})
