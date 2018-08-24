set(SOURCES
    alltoall_1d.f90
    factor.f
    fft235.f
    kernel.f
    pzfft3dv_mod.f
   )

add_library(${EXT_FFTE_LIB} STATIC ${SOURCES})
add_dependencies(${EXT_FFTE_LIB} ${SALMON_IO_LIB} ${SALMON_PARALLEL_LIB})
