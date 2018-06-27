set(SOURCES
    builtin_pz.f90
    builtin_pzm.f90
    builtin_pz_sp.f90
    builtin_pbe.f90
    builtin_tbmbj.f90
    salmon_xc.f90
   )

add_library(${SALMON_XC_LIB} STATIC ${SOURCES})
add_dependencies(${SALMON_XC_LIB} ${SALMON_PARALLEL_LIB} ${SALMON_IO_LIB} ${SALMON_MATH_LIB})
list(APPEND SALMON_LINK_LIBRARIES ${SALMON_XC_LIB})
set(SALMON_LINK_LIBRARIES ${SALMON_LINK_LIBRARIES} PARENT_SCOPE)
