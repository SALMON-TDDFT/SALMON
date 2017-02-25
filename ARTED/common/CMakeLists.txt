set(SOURCES
    total_energy.f90
    ion_force.f90
    Exc_Cor.f90
    Hartree.f90
    Ylm_dYlm.f90
    psi_rho.f90
    hpsi.f90
    hpsi_stencil.f90
    preprocessor.f90
    )

add_library(${COMMON_LIB} STATIC ${SOURCES})
