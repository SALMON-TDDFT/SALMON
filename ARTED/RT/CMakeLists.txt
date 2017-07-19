set(SOURCES
    dt_evolve.f90
    current.f90
    k_shift_wf.f90
    init_Ac.f90
    Fourier_tr.f90
    )

set(SOURCES_ACC
    acc/current.f90
    acc/hamiltonian.f90
    )

set(SOURCES_OMP
    hamiltonian.f90
    )

if (ENABLE_OPENACC)
  set(SOURCES ${SOURCES} ${SOURCES_ACC})
else ()
  set(SOURCES ${SOURCES} ${SOURCES_OMP})
endif ()

add_library(${ARTED_RT_LIB} STATIC ${SOURCES})
add_dependencies(${ARTED_RT_LIB} ${SALMON_MODULE_LIB} ${ARTED_MODULE_LIB} ${ARTED_COMMON_LIB})
