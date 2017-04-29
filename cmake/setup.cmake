list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/platform)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

string(COMPARE EQUAL "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" IN_SOURCE_BUILD)
if (${IN_SOURCE_BUILD})
  message(WARNING "CMake generates the build configurations within source directories (the in-source build). We recommend the out-of-source build.")
endif ()
