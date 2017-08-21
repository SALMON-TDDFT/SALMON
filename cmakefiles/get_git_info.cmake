# [Ref] CMake: Use Git Branch & Commit Details in Project
#       http://xit0.org/2013/04/cmake-use-git-branch-and-commit-details-in-project/

find_package(Git)

if (GIT_FOUND)
  # Get the current working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    RESULT_VARIABLE GIT_BRANCH_FOUND
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Get the latest abbreviated commit hash of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 --format=%H
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    RESULT_VARIABLE GIT_COMMIT_HASH_FOUND
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  if (${GIT_BRANCH_FOUND} EQUAL 0)
    message(STATUS "Git: ${GIT_COMMIT_HASH} in ${GIT_BRANCH}")
  else ()
    set(GIT_FOUND FALSE)
  endif ()
else ()
  message(STATUS "Git not found.")
endif ()
