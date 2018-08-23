set(SOURCES
    alltoall_1d.f90
    factor.f
    fft235.f
    kernel.f
    pzfft3dv_mod.f
   )

add_library(${SALMON_FFTE_LIB} STATIC ${SOURCES})
add_dependencies(${SALMON_FFTE_LIB} ${SALMON_PARALLEL_LIB} ${SALMON_IO_LIB})
list(APPEND SALMON_LINK_LIBRARIES ${SALMON_FFTE_LIB})
set(SALMON_LINK_LIBRARIES ${SALMON_LINK_LIBRARIES} PARENT_SCOPE)
