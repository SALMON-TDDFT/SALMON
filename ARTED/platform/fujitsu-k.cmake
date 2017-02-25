### Fujitsu Compiler for K computer @RIKEN AICS
set(TARGET_SUFFIX               ".cpu")

set(ARCH                        "")
set(SIMD_SET                    "HPC_ACE")
set(OPENMP_FLAGS                "-Kopenmp")
set(LAPACK_FLAGS                "-SSL2")
set(ADDITIONAL_MACRO            "")
set(ADDITIONAL_OPTIMIZE_FLAGS   "")

set(Fortran_FLAGS_General       "-Cpp -Kocl,nooptmsg")
set(C_FLAGS_General             "-Kocl,nooptmsg")

set(CMAKE_Fortran_COMPILER      "mpifrtpx")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -Kfast")
set(CMAKE_C_COMPILER            "mpifccpx")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g")
set(CMAKE_C_FLAGS_RELEASE       "-O3 -Kfast")


########
# Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Cross-compiling for Fujitsu SPARC64VIIIfx")
set(CMAKE_SYSTEM_PROCESSOR    "s64fx")
set(CMAKE_Fortran_COMPILER_ID "Fujitsu" CACHE STRING "Fujitsu MPI Fortran cross-compiler" FORCE)
set(CMAKE_C_COMPILER_ID       "Fujitsu" CACHE STRING "Fujitsu MPI C cross-compiler" FORCE)
set(CMAKE_Fortran_MODDIR_FLAG "-M ")
