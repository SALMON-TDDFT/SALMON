include(${CMAKE_SOURCE_DIR}/cmakefiles/create_test.cmake)

if (USE_MPI)
  find_package(MPI REQUIRED)
  set(TEST_MPI_NPROC 4)
  set(TEST_COMMAND "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${TEST_MPI_NPROC} ${CMAKE_BINARY_DIR}/${TARGET_NAME}")
else ()
  set(TEST_COMMAND "${CMAKE_BINARY_DIR}/${TARGET_NAME}")
endif ()
message(STATUS "testing command: ${TEST_COMMAND}")

add_custom_target("check_test_requirements" COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS ${TARGET_NAME})
add_custom_target("test_preparations" ALL
                  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "generate ctest scripts"
                  COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/pseudo" "${CMAKE_CURRENT_BINARY_DIR}/pseudo")

# test directories
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/001_bulk_gs")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/011_bulk_ms")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/101_C2H2_gs")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/102_C2H2_rt_response")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/103_C2H2_rt_pulse")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/111_bulk_Si_gs_dp")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/112_bulk_Si_rt_response_dp")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/113_bulk_Si_rt_pulse_dp")
