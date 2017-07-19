set(SOURCES
    allocate_mat.f90
    allocate_psl.f90
    calc_density.f90
    calc_invA_lapack.f90
    change_order.f90
    copyVlocal.f90
    copy_psi_mesh.f90
    deallocate_mat.f90
    gradient2.f90
    gradient.f90
    hpsi2.f90
    init_sendrecv.f90
    inner_product.f90
    laplacian2.f90
    new_world.f90
    readbox_rt.f90
    read_pslfile.f90
    rmmdiis_eigen_lapack.f90
    scf_data.f90
    sendrecv.f90
    sendrecv_groupob.f90
    sendrecv_groupob_ngp.f90
    sendrecvh.f90
    total_energy.f90
    writebox_rt.f90
    )
  #  share_mesh_1d_old.f90

add_library(${GCEED_MODULE_LIB} ${SOURCES})
add_dependencies(${GCEED_MODULE_LIB} ${SALMON_MODULE_LIB} ${GCEED_MISC_LIB})
