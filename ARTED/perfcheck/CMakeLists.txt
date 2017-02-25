set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FIX_Fortran_NOMAIN}")
set(TARGET stencilperf${TARGET_SUFFIX})
set(SOURCES
    wrap_variables.f90
    stencil_perf_check.f90
    main.c
    )

add_executable(${TARGET} EXCLUDE_FROM_ALL ${SOURCES})
add_dependencies(${TARGET} ${MODULE_LIB} ${COMMON_LIB} ${RT_LIB} ${STENCIL_LIB})
target_link_libraries(${TARGET} ${MODULE_LIB} ${COMMON_LIB} ${RT_LIB} ${STENCIL_LIB})

add_custom_target(stencilperf DEPENDS ${TARGET})
