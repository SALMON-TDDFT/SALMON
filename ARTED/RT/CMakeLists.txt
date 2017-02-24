set(SOURCES
    dt_evolve.f90
    dt_evolve_hpsi.f90
    current.f90
    k_shift_wf.f90
    hpsi_RT.f90
    init_Ac.f90
    Fourier_tr.f90
    )

if (ENABLE_OPENACC)
  set(SOURCES ${SOURCES} acc/current.f90)
endif ()

add_library(${RT_LIB} STATIC ${SOURCES})
