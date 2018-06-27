set(SOURCES
    add_polynomial.f90
    aft_fourier.f90
    calc_current.f90
    calc_env_trigon.f90
    calcAext.f90
    calcEstatic.f90
    calcVbox.f90
    check_ae_shape.f90
    convert_input_rt.f90
    dip.f90
    gradient_ex.f90
    hpsi_groupob.f90
    initA.f90
    projection.f90
    read_rt.f90
    real_time_dft.f90
    set_numcpu_rt.f90
    set_numcpu_rt_sp.f90
    set_vonf_sd.f90
    taylor_coe.f90
    taylor${HPSI_TEST_SUFFIX}.f90
    time_evolution_step.f90
    total_energy_groupob.f90
    xc_fast.f90
    )

add_library(${GCEED_RT_LIB} ${SOURCES})
add_dependencies(${GCEED_RT_LIB} ${GCEED_MODULE_LIB} ${GCEED_SCF_LIB})
