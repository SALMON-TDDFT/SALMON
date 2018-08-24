### Intel Compiler for Haswell, Broadwell...
set(TARGET_SUFFIX               ".cpu")

set(ARCH                        "-xCORE-AVX2")
set(SIMD_SET                    "AVX")
set(OPENMP_FLAGS                "-qopenmp")
set(LAPACK_FLAGS                "-mkl=parallel")
set(ScaLAPACK_FLAGS             "-mkl=cluster")
set(ADDITIONAL_MACRO            "")
set(ADDITIONAL_OPTIMIZE_FLAGS   "-ansi-alias -fno-alias")

set(Fortran_FLAGS_General       "-fpp -nogen-interface -std03 -warn all -diag-disable 6477,7025")
set(C_FLAGS_General             "-Wall -restrict")

set(CMAKE_Fortran_COMPILER      "mpiifort")
set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
set(CMAKE_C_COMPILER            "mpiicc")
set(CMAKE_C_FLAGS_DEBUG         "-O2 -g")
set(CMAKE_C_FLAGS_RELEASE       "-O3")

set(USE_MPI             ON)
set(REDUCE_FOR_MANYCORE ON)


########
# CMake Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Cross-compiling for Intel Haswell, Broadwell... (AVX2)")
set(CMAKE_SYSTEM_PROCESSOR "avx2")
