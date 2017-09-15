include(${CMAKE_SOURCE_DIR}/cmakefiles/create_test.cmake)

find_package(MPI REQUIRED) # find executable file for running MPI application.

set(TEST_MPI_NPROC   4)
set(TEST_MPI_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${TEST_MPI_NPROC} ${CMAKE_BINARY_DIR}/${TARGET_NAME}")
message(STATUS "testing MPI: ${TEST_MPI_COMMAND}")

add_custom_target("check_test_requirements" COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TARGET_NAME})
add_custom_target("test_preparations" ALL
                  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "generate ctest scripts"
                  COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/pseudo" "${CMAKE_CURRENT_BINARY_DIR}/pseudo")

# test directories
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/001_bulk_gs_rt")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/101_C2H2_gs")
