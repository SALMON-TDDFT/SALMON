add_subdirectory("./ARTED")
add_subdirectory("./GCEED")
add_subdirectory("./parallel")
add_subdirectory("./io")
add_subdirectory("./math")
#add_subdirectory("./gs")
#add_subdirectory("./rt")
#add_subdirectory("./atom")
add_subdirectory("./xc")
#add_subdirectory("./poisson")
#add_subdirectory("./maxwell")
add_subdirectory("./common")
add_subdirectory("./misc")

# Executable binary
add_executable(${TARGET_NAME} main.f90)
target_link_libraries(${TARGET_NAME} ${SALMON_LINK_LIBRARIES} ${EXTERNAL_LIBS})
add_dependencies(${TARGET_NAME} ${SALMON_LINK_LIBRARIES})

# Install binary
install(TARGETS ${TARGET_NAME} DESTINATION "bin")
