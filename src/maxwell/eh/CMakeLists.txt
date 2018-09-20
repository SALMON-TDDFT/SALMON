set(SOURCES
    init_eh.f90
    calc_eh.f90
    finalize_eh.f90
    )

add_library(${MAXWELL_EH_LIB} ${SOURCES})
add_dependencies(${MAXWELL_EH_LIB} ${SALMON_MAXWELL_LIB})
