###
### Compile stencil computation sources
###

### Soruces for current and total_energy
set(SOURCES_COMPILER_VEC
    F90/current.f90
    F90/total_energy.f90
    )
set(SOURCES_EXPLICIT_VEC
    C/${SIMD_SET}/current.c
    C/${SIMD_SET}/total_energy.c
    )

### Sources for hpsi
set(SOURCE_HPSI_ORIGINAL
    F90/hpsi_original.f90
    )
set(SOURCE_HPSI_F90_COMPILER_VEC
    F90/hpsi.f90
    )
set(SOURCE_HPSI_C_COMPILER_VEC
    C/hpsi.c
    )
set(SOURCE_HPSI_EXPLICIT_VEC
    C/${SIMD_SET}/hpsi.c
    )


###
if (OPT_STENCIL)
  if (STENCIL_WITH_C)
    if (ENABLE_EXPLICIT_VEC)
      set(SOURCES ${SOURCE_HPSI_EXPLICIT_VEC})
    else ()
      set(SOURCES ${SOURCE_HPSI_C_COMPILER_VEC})
    endif ()
  else ()
    set(SOURCES ${SOURCE_HPSI_F90_COMPILER_VEC})
  endif ()
else ()
  set(SOURCES ${SOURCE_HPSI_ORIGINAL})
endif ()

if (ENABLE_EXPLICIT_VEC)
  set(SOURCES ${SOURCES} ${SOURCES_EXPLICIT_VEC})
else ()
  set(SOURCES ${SOURCES} ${SOURCES_COMPILER_VEC})
endif ()

add_library(${STENCIL_LIB} STATIC ${SOURCES})
