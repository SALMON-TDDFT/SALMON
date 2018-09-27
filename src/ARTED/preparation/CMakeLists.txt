set(SOURCES
    init.f90
    init_wf.f90
    fd_coef.f90
    prep_ps.f90
    pp_postprocess.f90
    )

add_library(${ARTED_PREP_LIB} STATIC ${SOURCES})
add_dependencies(${ARTED_PREP_LIB} ${ARTED_MODULE_LIB})
