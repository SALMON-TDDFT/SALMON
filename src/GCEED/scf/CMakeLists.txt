set(SOURCES
    add_in_broyden.f90
    broyden.f90
    calc_dos.f90
    calc_occupation.f90
    calc_pdos.f90
    calc_rho_in.f90
    check_dos_pdos.f90
    convert_input_scf.f90
    copy_density.f90
    deallocate_sendrecv.f90
    eigen_subdiag_${FILE_MATHLIB}.f90
    gram_schmidt.f90
    prep_ini.f90
    real_space_dft.f90
    rmmdiis_eigen.f90
    rmmdiis.f90
    set_numcpu_scf.f90
    simple_mixing.f90
    structure_opt.f90
    subdgemm_lapack.f90
    subspace_diag.f90
    ybcg.f90
    )

add_library(${GCEED_SCF_LIB} ${SOURCES})
add_dependencies(${GCEED_SCF_LIB} ${SALMON_MODULE_LIB} ${GCEED_MODULE_LIB})
