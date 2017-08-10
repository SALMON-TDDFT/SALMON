set(SOURCES
    beam.f90
    FDTD.f90
    )

add_library(${ARTED_FDTD_LIB} ${SOURCES})
add_dependencies(${ARTED_FDTD_LIB} ${ARTED_MODULE_LIB} ${SALMON_MATH_LIB})
