set(SOURCES
    init.f90
    init_wf.f90
    fd_coef.f90
    input_ps.f90
    prep_ps.f90
    )

add_library(${PREP_LIB} STATIC ${SOURCES})
