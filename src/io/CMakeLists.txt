set(SOURCES
    inputoutput.f90
    salmon_file.f90
    salmon_global.f90
   )

add_library(${SALMON_IO_LIB} STATIC ${SOURCES})
add_dependencies(${SALMON_IO_LIB} ${SALMON_PARALLEL_LIB})
list(APPEND SALMON_LINK_LIBRARIES ${SALMON_IO_LIB})
set(SALMON_LINK_LIBRARIES ${SALMON_LINK_LIBRARIES} PARENT_SCOPE)
