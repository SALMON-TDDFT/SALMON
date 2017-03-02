set(SOURCES
    total_energy.f90
    ion_force.f90
    Exc_Cor.f90
    Hartree.f90
    Ylm_dYlm.f90
    psi_rho.f90
    hpsi.f90
    preprocessor.f90
    )

if (ENABLE_OPENACC)
  set(SOURCES ${SOURCES} acc/hpsi.f90)
endif ()

add_library(${COMMON_LIB} STATIC ${SOURCES})
