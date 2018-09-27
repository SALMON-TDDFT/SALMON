###
### Compile fortran module sources
###
set(SOURCES
    global_variables.f90
    opt_variables.f90
    performance_analyzer.f90
    nvtx.f90
    )
#if (USE_NVTX)
#   set(SOURCES ${SOURCES} nvtx.f90)
#endif ()
add_library(${ARTED_MODULE_LIB} STATIC ${SOURCES})
add_dependencies(${ARTED_MODULE_LIB} ${SALMON_IO_LIB} ${SALMON_COMMON_LIB} ${SALMON_XC_LIB} ${SALMON_ATOM_LIB})
