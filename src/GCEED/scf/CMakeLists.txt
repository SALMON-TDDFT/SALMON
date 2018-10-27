set(SOURCES
    buffer_broyden_ns.f90
    calc_dos.f90
    calc_occupation.f90
    calc_pdos.f90
    check_dos_pdos.f90
    convert_input_scf.f90
    copy_density.f90
    deallocate_sendrecv.f90
    eigen_subdiag_${FILE_MATHLIB}.f90
    eigen_subdiag_periodic_${FILE_MATHLIB}.f90
    inner_product5.f90
    gram_schmidt.f90
    gram_schmidt_periodic.f90
    gscg.f90
    gscg_periodic.f90
    prep_ini.f90
    real_space_dft.f90
    rmmdiis_eigen.f90
    rmmdiis.f90
    set_numcpu_scf.f90
    simple_mixing.f90
    subdgemm_lapack.f90
    subspace_diag.f90
    subspace_diag_periodic.f90
    total_energy_periodic_scf.f90
    total_energy_periodic_scf_esp.f90
    write_eigen.f90
    ybcg.f90
    ybcg_periodic.f90
    )

add_library(${GCEED_SCF_LIB} ${SOURCES})
add_dependencies(${GCEED_SCF_LIB} ${GCEED_MODULE_LIB})
