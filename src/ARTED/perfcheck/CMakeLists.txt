set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FIX_Fortran_NOMAIN}")
set(TARGET stencilperf${TARGET_SUFFIX})
set(SOURCES
    wrap_variables.f90
    stencil_perf_check.f90
    main.c
    )

add_executable(${TARGET} EXCLUDE_FROM_ALL ${SOURCES})
add_dependencies(${TARGET} ${ARTED_MODULE_LIB} ${ARTED_COMMON_LIB} ${ARTED_RT_LIB} ${ARTED_STENCIL_LIB})
target_link_libraries(${TARGET} ${ARTED_MODULE_LIB} ${ARTED_COMMON_LIB} ${ARTED_RT_LIB} ${ARTED_STENCIL_LIB})

add_custom_target(stencilperf DEPENDS ${TARGET})
