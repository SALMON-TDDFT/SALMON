set(SOURCES
    io_gs_wfn_k.f90
    Gram_Schmidt.f90
    diag.f90
    sp_energy.f90
    Occupation_Redistribution.f90
    Density_Update.f90
    CG.f90
    write_GS_data.f90
    Fermi_Dirac_distribution.f90
    )

add_library(${ARTED_GS_LIB} STATIC ${SOURCES})
add_dependencies(${ARTED_GS_LIB} ${SALMON_MODULE_LIB} ${ARTED_MODULE_LIB} ${ARTED_COMMON_LIB})
