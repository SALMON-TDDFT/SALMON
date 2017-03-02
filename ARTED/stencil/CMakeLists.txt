###
### Compile stencil computation sources
###

### Soruces for current and total_energy
set(SOURCES_EXPLICIT_VEC
    C/${SIMD_SET}/current.c
    C/${SIMD_SET}/total_energy.c
    )
set(SOURCES_COMPILER_VEC
    F90/total_energy.f90
    F90/current.f90
    )
set(SOURCES_ACC
    F90/acc/current.f90
    F90/acc/hpsi.f90
    )


### Sources for hpsi
set(SOURCE_HPSI_EXPLICIT_VEC   C/${SIMD_SET}/hpsi.c)
set(SOURCE_HPSI_COMPILER_VEC   F90/hpsi.f90)
set(SOURCE_HPSI_ORIGINAL       F90/hpsi_original.f90)

if (ENABLE_EXPLICIT_VEC)
  set(SOURCE_HPSI ${SOURCE_HPSI_EXPLICIT_VEC})
else ()
  set(SOURCE_HPSI ${SOURCE_HPSI_COMPILER_VEC})
endif ()

if (OPT_STENCIL)
else ()
  set(SOURCE_HPSI ${SOURCE_HPSI_ORIGINAL})
endif ()


###
if (ENABLE_EXPLICIT_VEC)
  set(SOURCES_OMP ${SOURCE_HPSI} ${SOURCES_EXPLICIT_VEC})
else ()
  set(SOURCES_OMP ${SOURCE_HPSI} ${SOURCES_COMPILER_VEC})
endif ()

if (ENABLE_OPENACC)
  set(SOURCES ${SOURCES_ACC} ${SOURCE_HPSI_COMPILER_VEC} ${SOURCES_COMPILER_VEC})
else ()
  set(SOURCES ${SOURCES_OMP})
endif ()

add_library(${STENCIL_LIB} STATIC ${SOURCES})
