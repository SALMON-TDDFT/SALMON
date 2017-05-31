set(GCEED_LIB      gceed)
set(MODULE_LIBDIR  modules)
set(MISC_LIBDIR    misc)
set(COMMON_LIBDIR  common)
set(SCF_LIBDIR     scf)
set(RT_LIBDIR      rt)

set(MODULE_LIB  ${GCEED_LIB}_${MODULE_LIBDIR})
set(MISC_LIB    ${GCEED_LIB}_${MISC_LIBDIR})
set(COMMON_LIB  ${GCEED_LIB}_${COMMON_LIBDIR})
set(SCF_LIB     ${GCEED_LIB}_${SCF_LIBDIR})
set(RT_LIB      ${GCEED_LIB}_${RT_LIBDIR})


### Build target
add_subdirectory(${MODULE_LIBDIR})
add_subdirectory(${MISC_LIBDIR})
add_subdirectory(${COMMON_LIBDIR})
add_subdirectory(${SCF_LIBDIR})
add_subdirectory(${RT_LIBDIR})
add_library(${GCEED_LIB} STATIC gceed.f90 read_input_gceed.f90)

add_dependencies(${MODULE_LIB} ${MISC_LIB})
add_dependencies(${COMMON_LIB} ${MODULE_LIB})
add_dependencies(${SCF_LIB}    ${MODULE_LIB})
add_dependencies(${RT_LIB}     ${MODULE_LIB} ${SCF_LIB})
add_dependencies(${GCEED_LIB}  ${MODULE_LIB} ${COMMON_LIB} ${SCF_LIB} ${RT_LIB})


set(LIBRARIES
    ${GCEED_LIB}
    ${RT_LIB}
    ${SCF_LIB}
    ${COMMON_LIB}
    ${MODULE_LIB}
    ${MISC_LIB}
    )
list(APPEND SALMON_LINK_LIBRARIES ${LIBRARIES})
set(SALMON_LINK_LIBRARIES ${SALMON_LINK_LIBRARIES} PARENT_SCOPE)
