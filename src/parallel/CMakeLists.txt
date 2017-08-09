set(SOURCES
    salmon_parallel.f90
)

if (USE_MPI)
  set(SOURCES ${SOURCES} salmon_communication.f90)
else ()
  set(SOURCES ${SOURCES} salmon_communication_dummy.f90)
endif ()

add_library(${SALMON_PARALLEL_LIB} STATIC ${SOURCES})
list(APPEND SALMON_LINK_LIBRARIES ${SALMON_PARALLEL_LIB})
set(SALMON_LINK_LIBRARIES ${SALMON_LINK_LIBRARIES} PARENT_SCOPE)
