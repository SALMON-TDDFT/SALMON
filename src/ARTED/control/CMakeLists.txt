set(SOURCES
    control_sc.f90
    control_ms.f90
    inputfile.f90
    initialization.f90
    ground_state.f90
    )

add_library(${ARTED_CONTROL_LIB} STATIC ${SOURCES})
add_dependencies(${ARTED_CONTROL_LIB} ${SALMON_MODULE_LIB} ${ARTED_MODULE_LIB} ${ARTED_PREP_LIB} ${ARTED_GS_LIB} ${ARTED_RT_LIB} ${ARTED_COMMON_LIB} ${ARTED_STENCIL_LIB} ${ARTED_FDTD_LIB})
