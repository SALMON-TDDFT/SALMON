###
### Compile fortran module sources
###
set(SOURCES
    global_variables.f90
    timer.f90
    env_variables.f90
    opt_variables.f90
    performance_analyzer.f90
    env_variables_internal.c
    papi_wrap.c
    nvtx.f90
    misc_routines.f90
    backup_routines.f90
    )
#if (USE_NVTX)
#   set(SOURCES ${SOURCES} nvtx.f90)
#endif ()
add_library(${ARTED_MODULE_LIB} STATIC ${SOURCES})
add_dependencies(${ARTED_MODULE_LIB} ${SALMON_MODULE_LIB})
