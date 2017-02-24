set(SOURCES
    Gram_Schmidt.f90
    diag.f90
    sp_energy.f90
    Occupation_Redistribution.f90
    Density_Update.f90
    CG.f90
    write_GS_data.f90
    Fermi_Dirac_distribution.f90
    )

add_library(${GS_LIB} STATIC ${SOURCES})
