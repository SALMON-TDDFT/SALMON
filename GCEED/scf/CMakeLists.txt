set(SOURCES
    calc_dos.f90
    calc_occupation.f90
    calc_pdos.f90
    calc_rho_in.f90
    copy_density.f90
    deallocate_sendrecv.f90
    eigen_subdiag_${FILE_MATHLIB}.f90
    gram_schmidt.f90
    prep_ini.f90
    read_input_scf.f90
    real_space_dft.f90
    rmmdiis_eigen.f90
    rmmdiis.f90
    simple_mixing.f90
    structure_opt.f90
    subdgemm_lapack.f90
    subspace_diag.f90
    ybcg.f90
    )

add_library(${SCF_LIB} ${SOURCES})
