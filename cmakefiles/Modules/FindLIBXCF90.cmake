#
#  Copyright 2018 SALMON developers
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#-----------------------------------------------------------------------------------------

# FindLIBXCF90.cmake
# ----------
# This module finds a installed binary and Fortran90 includes of 
# Libxc (http://www.tddft.org/programs/libxc/) package.
#
# This module sets the following variables:
#  LIBXCF90_INCLUDE_DIRS
#  LIBXCF90_LIBRARIES
#  LIBXCF90_DEFINITIONS
#  LIBXCF90_FOUND

find_package(PkgConfig)

# Detect libxc setting
pkg_check_modules(PC_LIBXCF90 QUIET libxc)

set(_libxc_hints 
  ${PC_LIBXCF90_INCLUDEDIR} 
  ${PC_LIBXCF90_INCLUDE_DIRS}
)

# Search-path priority
set(_libxc_paths 
  ~/program            # Uemoto's environment :P
  ~/local              #
  /usr/local/Celler    # Homebrew
  /opt/local           # MacPorts
  ~/Library/Frameworks # MacOS
  /Library/Frameworks  #
  /sw                  # Flink
  /usr                 # General linux package manager
  /opt                 #
)

# Search 'include' directory containing fortran90 module file
find_path(
  LIBXCF90_INCLUDE_DIR 
  xc_f90_types_m.mod HINTS ${_libxc_hints} PATHS ${_libxc_paths}
)

# Search libxcf90 (fortran90 library)
find_library(
  LIBXCF90_LIBRARY
  NAMES xcf90 HINTS ${_libxc_hints} PATHS ${_libxc_paths}
)

# Search libxc (general libxc library)
find_library(
  LIBXC_LIBRARY
  NAMES xc HINTS ${_libxc_hints} PATHS ${_libxc_paths}
)

# Show error messages
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  LIBXCF90 DEFAULT_MSG
  LIBXCF90_INCLUDE_DIR LIBXCF90_LIBRARY LIBXCF90_INCLUDE_DIR
)


mark_as_advanced(LIBXCF90_INCLUDE_DIR LIBXCF90_LIBRARY LIBXC_LIBRARY)

# Store results
set(LIBXCF90_LIBRARIES ${LIBXCF90_LIBRARY} ${LIBXC_LIBRARY})
set(LIBXCF90_INCLUDE_DIRS ${LIBXCF90_INCLUDE_DIR})
set(LIBXCF90_DEFINITIONS ${PC_LIBXCF90_CFLAGS_OTHER})
