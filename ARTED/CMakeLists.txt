set(TARGET_NAME   arted)
set(MODULE_LIB    modules)
set(PREP_LIB      preparation)
set(GS_LIB        GS)
set(RT_LIB        RT)
set(COMMON_LIB    common)
set(STENCIL_LIB   stencil)
set(FDTD_LIB      FDTD)
set(PERFCHECK_EXT perfcheck)
set(CONTROL_LIB   control)


### Build target
add_subdirectory(${MODULE_LIB})
add_subdirectory(${PREP_LIB})
add_subdirectory(${GS_LIB})
add_subdirectory(${RT_LIB})
add_subdirectory(${COMMON_LIB})
add_subdirectory(${STENCIL_LIB})
add_subdirectory(${PERFCHECK_EXT})
add_subdirectory(${FDTD_LIB})
add_subdirectory(${CONTROL_LIB})

add_library(${TARGET_NAME} STATIC main.f90)

add_dependencies(${TARGET_NAME} ${CONTROL_LIB})
add_dependencies(${CONTROL_LIB} 
                 ${MODULE_LIB}
                 ${PREP_LIB}
                 ${GS_LIB}
                 ${RT_LIB}
                 ${COMMON_LIB}
                 ${STENCIL_LIB}
                 ${FDTD_LIB}
                 )
add_dependencies(${PREP_LIB}    ${MODULE_LIB})
add_dependencies(${GS_LIB}      ${MODULE_LIB})
add_dependencies(${RT_LIB}      ${MODULE_LIB})
add_dependencies(${COMMON_LIB}  ${MODULE_LIB})
add_dependencies(${STENCIL_LIB} ${MODULE_LIB})
add_dependencies(${FDTD_LIB}    ${MODULE_LIB})

set(LIBRARIES
    ${CONTROL_LIB}
    ${PREP_LIB}
    ${GS_LIB}
    ${RT_LIB}
    ${COMMON_LIB}
    ${STENCIL_LIB}
    ${MODULE_LIB}
    ${FDTD_LIB}
    )

target_link_libraries(${TARGET_NAME} ${LIBRARIES})

list(APPEND SALMON_LINK_LIBRARIES ${TARGET_NAME})
set(SALMON_LINK_LIBRARIES ${SALMON_LINK_LIBRARIES} PARENT_SCOPE)
