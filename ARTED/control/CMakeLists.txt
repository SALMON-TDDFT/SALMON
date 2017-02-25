set(SOURCES
    control_sc.f90
    control_ms.f90
    inputfile.f90
    )

add_library(${CONTROL_LIB} STATIC ${SOURCES})
