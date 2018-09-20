set(SOURCES
    init_weyl.f90
    calc_weyl.f90
    finalize_weyl.f90
    )

add_library(${MAXWELL_WEYL_LIB} ${SOURCES})
add_dependencies(${MAXWELL_WEYL_LIB} ${SALMON_MAXWELL_LIB})
