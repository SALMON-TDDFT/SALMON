set(SOURCES
    salmon_file.f90
    inputoutput.f90
    salmon_global.f90
    salmon_parallel.f90
    salmon_math.f90
)

if (USE_MPI)
  set(SOURCES ${SOURCES} salmon_communication.f90)
else ()
  set(SOURCES ${SOURCES} salmon_communication_dummy.f90)
endif ()

add_library(${SALMON_MODULE_LIB} STATIC ${SOURCES})
add_dependencies(${SALMON_MODULE_LIB} ${SALMON_CORE_LIB})
list(APPEND SALMON_LINK_LIBRARIES ${SALMON_MODULE_LIB})
set(SALMON_LINK_LIBRARIES ${SALMON_LINK_LIBRARIES} PARENT_SCOPE)
