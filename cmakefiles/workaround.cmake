####
# for Intel Fortran
if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  set(FIX_Fortran_NOMAIN "-nofor-main")
endif ()

####
# for PGI OpenACC
#   In several systems, CMake may links not-parallelizable OpenACC libraries explicitly, but this will be causes linker error.
#     (`parallelizable` means that API can be called within the multiple OpenMP threads concurrently.)
#   PGI compiler implicitly links the libraries when enabling OpenACC feature.
#   We remove OpenACC libraries from variable `CMAKE_C_IMPLICIT_LINK_LIBRARIES`.
set(_OPENACC_LIBRARIES "accapi;accg;accn;accg2;dl;cudadevice;")
foreach(_accLib ${_OPENACC_LIBRARIES})
  list(REMOVE_ITEM CMAKE_C_IMPLICIT_LINK_LIBRARIES ${_accLib})
endforeach()
