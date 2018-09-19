set(SOURCES
    init_coulomb.f90
    calc_coulomb.f90
    finalize_coulomb.f90
    )

add_library(${MAXWELL_COULOMB_LIB} ${SOURCES})
add_dependencies(${MAXWELL_COULOMB_LIB} ${SALMON_MAXWELL_LIB})
