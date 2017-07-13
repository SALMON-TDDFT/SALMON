set(TARGET_NAME   salmon_module)

set(SOURCES
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

add_library(${TARGET_NAME} STATIC ${SOURCES})
list(APPEND SALMON_LINK_LIBRARIES ${TARGET_NAME})
set(SALMON_LINK_LIBRARIES ${SALMON_LINK_LIBRARIES} PARENT_SCOPE)
