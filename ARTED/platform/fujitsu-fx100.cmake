### Fujitsu Compiler, FX100 system
set(TARGET_SUFFIX               ".cpu")

set(ARCH                        "-KHPC_ACE2")
set(SIMD_SET                    "HPC_ACE2")
set(OPENMP_FLAGS                "-Kopenmp")
set(LAPACK_FLAGS                "-SSL2")
set(ADDITIONAL_MACRO            "")
set(ADDITIONAL_OPTIMIZE_FLAGS   "")

set(Fortran_FLAGS_General       "-Cpp -Kocl,nooptmsg")
set(C_FLAGS_General             "-Kocl,nooptmsg -Xg -std=gnu99")

set(CMAKE_Fortran_COMPILER      "mpifrtpx")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -Kfast,simd=1")
set(CMAKE_C_COMPILER            "mpifccpx")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g")
set(CMAKE_C_FLAGS_RELEASE       "-O3 -Kfast,simd=1")

set(ENABLE_REDUCE_FOR_MANYCORE ON)


########
# Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Cross-compiling for Fujitsu SPARC64 XIfx")
set(CMAKE_SYSTEM_PROCESSOR    "s64fx")
set(CMAKE_Fortran_COMPILER_ID "Fujitsu" CACHE STRING "Fujitsu MPI Fortran cross-compiler" FORCE)
set(CMAKE_C_COMPILER_ID       "Fujitsu" CACHE STRING "Fujitsu MPI C cross-compiler" FORCE)
set(CMAKE_Fortran_MODDIR_FLAG "-M ")
