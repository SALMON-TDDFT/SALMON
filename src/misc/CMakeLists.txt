set(SOURCES
    timer.f90
    misc_routines.f90
    backup_routines.f90
    unusedvar.f90
   )

add_library(${SALMON_MISC_LIB} STATIC ${SOURCES})
list(APPEND SALMON_LINK_LIBRARIES ${SALMON_MISC_LIB})
set(SALMON_LINK_LIBRARIES ${SALMON_LINK_LIBRARIES} PARENT_SCOPE)
