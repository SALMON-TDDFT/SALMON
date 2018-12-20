set(SOURCES
    salmon_maxwell.f90
    classic_em.f90
    weyl_init.f90
    weyl_finalize.f90
    weyl_calc.f90
    coulomb_init.f90
    coulomb_finalize.f90
    coulomb_calc.f90
    eh_calc.f90
    eh_init.f90
    eh_finalize.f90
    )

add_library(${SALMON_MAXWELL_LIB} STATIC ${SOURCES})
add_dependencies(${SALMON_MAXWELL_LIB} ${SALMON_PARALLEL_LIB} ${SALMON_IO_LIB} ${SALMON_MISC_LIB} ${GCEED_LIB})
list(APPEND SALMON_LINK_LIBRARIES ${SALMON_MAXWELL_LIB})
set(SALMON_LINK_LIBRARIES ${SALMON_LINK_LIBRARIES} PARENT_SCOPE)
