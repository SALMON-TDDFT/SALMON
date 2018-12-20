!
!  Copyright 2017 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module inputoutput
  use salmon_global
  implicit none
!Physical constant
  real(8),parameter :: au_time_fs = 0.02418884326505d0
  real(8),parameter :: au_energy_ev = 27.21138505d0
  real(8),parameter :: au_length_aa = 0.52917721067d0


  integer :: fh_variables_log
  integer :: fh_namelist
!  integer, parameter :: fh_atomic_spiecies = 902
  integer :: fh_atomic_coor
  integer :: fh_reentrance
  integer :: fh_atomic_red_coor
  logical :: if_nml_coor, if_nml_red_coor

  integer :: inml_calculation
  integer :: inml_control
  integer :: inml_units
  integer :: inml_parallel
  integer :: inml_system
  integer :: inml_pseudo
  integer :: inml_functional
  integer :: inml_rgrid
  integer :: inml_kgrid
  integer :: inml_tgrid
  integer :: inml_propagation
  integer :: inml_scf
  integer :: inml_emfield
  integer :: inml_multiscale
  integer :: inml_maxwell
  integer :: inml_analysis
  integer :: inml_hartree
  integer :: inml_ewald
  integer :: inml_opt
  integer :: inml_md
  integer :: inml_misc
  integer :: inml_group_fundamental
  integer :: inml_group_parallel
  integer :: inml_group_hartree
  integer :: inml_group_file
  integer :: inml_group_others

!! === old variables: will be removed after some time
  integer :: inml_group_function
!  integer :: inml_units
!  integer :: inml_control
!  integer :: inml_system
  integer :: inml_incident
!  integer :: inml_propagation
!  integer :: inml_rgrid
!  integer :: inml_kgrid
  integer :: inml_tstep
  integer :: inml_electrons
!  integer :: inml_pseudo
  integer :: inml_response
!  integer :: inml_multiscale
  integer :: inml_group_atom
!! === old variables: will be removed after some time

!Input/Output units
  integer :: iflag_unit_time
  integer,parameter :: ntype_unit_time_au = 0
  integer,parameter :: ntype_unit_time_fs = 1
  real(8) :: utime_to_au, utime_from_au

  integer :: iflag_unit_length
  integer,parameter :: ntype_unit_length_au = 0
  integer,parameter :: ntype_unit_length_aa = 1
  real(8) :: ulength_to_au, ulength_from_au

  integer :: iflag_unit_energy
  integer,parameter :: ntype_unit_energy_au = 0
  integer,parameter :: ntype_unit_energy_ev = 1
  real(8) :: uenergy_to_au, uenergy_from_au

  integer :: iflag_unit_charge
  integer,parameter :: ntype_unit_charge_au = 0
  real(8) :: ucharge_to_au, ucharge_from_au



  type unit_t
     character(32) :: name
     real(8)       :: conv
  end type unit_t

  type(unit_t) :: t_unit_energy
  type(unit_t) :: t_unit_energy_inv
  type(unit_t) :: t_unit_time
  type(unit_t) :: t_unit_time_inv
  type(unit_t) :: t_unit_current
  type(unit_t) :: t_unit_ac
  type(unit_t) :: t_unit_elec

contains
  subroutine read_input
    implicit none

    call read_stdin
    call read_input_common ! Should be renamed properly later
    if(restart_option == 'restart')return
    call read_atomic_coordinates
    call dump_input_common ! Should be renamed properly later
    call check_bad_input

  end subroutine read_input


  subroutine read_stdin
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    use salmon_file, only: get_filehandle
    implicit none

    integer :: cur 
    integer :: ret = 0
    character(100) :: buff, text
        

    if (comm_is_root(nproc_id_global)) then
      fh_namelist = get_filehandle()
      open(fh_namelist, file='.namelist.tmp', status='replace')
!      open(fh_atomic_spiecies, file='.atomic_spiecies.tmp', status='replace')
      if_nml_coor =.false.
      fh_atomic_coor = get_filehandle()
      open(fh_atomic_coor, file='.atomic_coor.tmp', status='replace')
      if_nml_red_coor = .false.
      fh_atomic_red_coor = get_filehandle()
      open(fh_atomic_red_coor, file='.atomic_red_coor.tmp', status='replace')
      fh_reentrance = get_filehandle()
      open(fh_reentrance, file='.reenetrance.tmp', status='replace')
      
      cur = fh_namelist
      do while (.true.)
        read(*, '(a)', iostat=ret) buff
        if (ret < 0) then
          exit
        else
          text = trim(adjustl(buff))
          ! Comment lines
          if (text(1:1) == '!') cycle
!          ! Beginning of 'atomic_species' part
!          if (text == '&atomic_spiecies') then
!            cur = fh_atomic_spiecies
!            cycle
!          end if
          ! Beginning of 'atomic_positions' part
          if (text == '&atomic_coor') then
            cur = fh_atomic_coor
            if_nml_coor =.true. 
            cycle
          end if
          if (text == '&atomic_red_coor') then
            cur = fh_atomic_red_coor
            if_nml_red_coor = .true.
            cycle
          end if
          ! Beginning of 'atomic_species' part
          if (text == '&reentrance') then
            cur = fh_reentrance
            cycle
          end if
          ! End of 'atomic_(spiecies|positions)' part
          if ((text == '/') .and. (cur /= fh_namelist)) then
            cur = fh_namelist
            cycle
          end if
          
          write(cur, '(a)') text
        end if
      end do
      close(fh_namelist)
      close(fh_atomic_coor)
      close(fh_atomic_red_coor)
!      close(fh_atomic_spiecies)
      close(fh_reentrance)
    end if

!    call comm_sync_all()


    return
  end subroutine read_stdin

  subroutine read_input_common
    use salmon_parallel
    use salmon_communication
    use salmon_file, only: get_filehandle
    implicit none
    integer :: ii
    real(8) :: norm

    namelist/calculation/ &
      & theory, &
      & calc_mode, &
      & use_ehrenfest_md, &
      & use_adiabatic_md, &
      & use_ms_maxwell, &
      & use_force, &
      & use_geometry_opt

    namelist/control/ &
      & restart_option, &
      & backup_frequency, &
      & time_shutdown, &
      & sysname, &
      & directory, &
      & dump_filename, &
      & modify_gs_wfn_k, &  !changed from modify_initial_guess
      & read_gs_wfn_k,   &  !changed from read_initial_guess
      & read_rt_wfn_k,   &
      & write_gs_wfn_k,  &
      & write_rt_wfn_k,  &
      & read_gs_wfn_k_ms,  &
      & read_rt_wfn_k_ms,  &
      & write_gs_wfn_k_ms, &
      & write_rt_wfn_k_ms


    namelist/units/ &
      & unit_system

    namelist/parallel/ &
      & domain_parallel, &
      & nproc_k, &
      & nproc_ob, &
      & nproc_domain, &
      & nproc_domain_s, &
      & num_datafiles_in, &
      & num_datafiles_out

    namelist/system/ &
      & iperiodic, &
      & ispin, &
      & al, &
      & isym, &
      & crystal_structure, &
      & nstate, &
      & nstate_spin, &
      & nelec, &
      & nelec_spin, &
      & temperature, &
      & temperature_k, &
      & nelem, &
      & natom, &
      & file_atom_coor, &
      & file_atom_red_coor

    namelist/pseudo/ &
      & pseudo_file, &
      & lmax_ps, &
      & lloc_ps, &
      & izatom, &
      & psmask_option, &
      & alpha_mask, &
      & gamma_mask, &
      & eta_mask

    namelist/functional/ &
      & xc, &
      & cname, &
      & xname, &
#ifdef SALMON_USE_LIBXC
      & alibx, &
      & alibc, &
      & alibxc, &
#endif
      & cval, &
      & no_update_func


    namelist/rgrid/ &
      & dl, &
      & num_rgrid

    namelist/kgrid/ &
      & num_kgrid, &
      & file_kw

    namelist/tgrid/ &
      & nt, &
      & dt

    namelist/propagation/ &
      & n_hamil, &
      & propagator

    namelist/scf/ &
      & amin_routine, &
      & ncg, &
      & amixing, &
      & rmixrate, &
      & nmemory_mb, &
      & alpha_mb, &
      & fsset_option, &
      & nfsset_start, &
      & nfsset_every, &
      & nscf, &
      & ngeometry_opt, &
      & subspace_diagonalization, &
      & convergence, &
      & threshold, &
      & threshold_norm_rho, &
      & threshold_norm_pot, &
      & omp_loop, &
      & skip_gsortho, &
      & iditer_notemperature, &
      & gscg

    namelist/emfield/ &
      & trans_longi, &
      & ae_shape1, &
      & e_impulse, &
      & amplitude1, &
      & rlaser_int_wcm2_1, &
      & pulse_tw1, &
      & omega1, &
      & epdir_re1, &
      & epdir_im1, &
      & phi_cep1, &
      & ae_shape2, &
      & amplitude2, &
      & rlaser_int_wcm2_2, &
      & pulse_tw2, &
      & omega2, &
      & epdir_re2, &
      & epdir_im2, &
      & phi_cep2, &
      & t1_t2, &
      & t1_delay, &
      & quadrupole, &
      & quadrupole_pot, &
      & alocal_laser , &
      & rlaserbound_sta , &
      & rlaserbound_end , &
      & nump , &
      & vecp , &
      & coop , &
      & radp_diele

    namelist/multiscale/ &
      & fdtddim, &
      & twod_shape, &
      & nx_m, &
      & ny_m, &
      & nz_m, &
      & hx_m, &
      & hy_m, &
      & hz_m, &
      & nksplit, &
      & nxysplit, &
      & nxvacl_m, &
      & nxvacr_m, &
      & nx_origin_m, &
      & ny_origin_m, &
      & nz_origin_m, &
      & file_macropoint, &
      & num_macropoint,  &
      & set_ini_coor_vel,&
      & nmacro_write_group

    namelist/maxwell/ &
      & al_em, &
      & dl_em, &
      & dt_em, &
      & nt_em, &
      & iboundary, &
      & wave_input, &
      & ek_dir1, &
      & source_loc1, &
      & ek_dir2, &
      & source_loc2, &
      & iobs_num_em, &
      & iobs_samp_em, &
      & obs_loc_em, &
      & shape_file, &
      & imedia_num, &
      & type_media, &
      & epsilon, &
      & rmu, &
      & sigma, &
      & omega_p_d, &
      & gamma_d, &
      & smooth_d, &
      & weight_d, &
      & wf_em

    namelist/analysis/ &
      & projection_option, &
      & projection_decomp, &
      & nenergy, &
      & de, &
      & out_psi, &
      & out_dos, &
      & out_pdos, &
      & out_dns, &
      & out_dos_start, &
      & out_dos_end, &
      & iout_dos_nenergy, &
      & out_dos_smearing, &
      & out_dos_method, &
      & out_dos_fshift, &
      & out_old_dns, &
      & out_dns_rt, &
      & out_dns_rt_step, &
      & out_dns_trans, &
      & out_dns_trans_energy, &
      & out_elf, &
      & out_elf_rt, &
      & out_elf_rt_step, &
      & out_estatic_rt, &
      & out_estatic_rt_step, &
      & out_rvf_rt, &
      & out_rvf_rt_step, &
      & out_tm, &
      & out_projection_step, &
      & out_ms_step, &
      & format3d, &
      & numfiles_out_3d, &
      & timer_process

    namelist/hartree/ &
      & meo, &
      & num_pole_xyz

    namelist/ewald/ &
      & newald, &
      & aewald

    namelist/opt/ &
      & cg_alpha_ini, &     !not use now if flag_use_grad_wf_on_force=.T.
      & cg_alpha_up, &
      & cg_alpha_down, &
      & convrg_scf_force, &
      & convrg_scf_ene, &
      & convrg_opt_fmax, &
      & convrg_opt_ene      !not use now if flag_use_grad_wf_on_force=.T.

    namelist/md/ &
      & ensemble, &
      & thermostat, &
      & step_velocity_scaling, &
      & step_update_ps, &
      & step_update_ps2,&
      & temperature0_ion, &
      & set_ini_velocity, &
      & file_ini_velocity, &
      & file_set_shake, &
      & thermostat_tau, &
      & friction, &
      & stop_system_momt

    namelist/misc/ &
      & fourier

    namelist/group_fundamental/ &
      & iditerybcg, &
      & iditer_nosubspace_diag, &
      & ntmg, &
      & idisnum, &
      & iwrite_projection, &
      & itwproj, &
      & iwrite_projnum, &
      & itcalc_ene 

    namelist/group_parallel/ &
      & isequential, &
      & imesh_s_all, &
      & iflag_comm_rho

    namelist/group_hartree/ &
      & hconv, &
      & lmax_meo

    namelist/group_file/ &
      & ic, &
      & oc, &
      & ic_rt, &
      & oc_rt

    namelist/group_others/ &
      & iparaway_ob, &
      & iscf_order, &
      & iswitch_orbital_mesh, &
      & iflag_psicube, &
      & lambda1_diis, &
      & lambda2_diis, &
      & file_ini, &
      & iparaway_ob, &
      & num_projection, &
      & iwrite_projection_ob, &
      & iwrite_projection_k, &
      & filename_pot, &
      & iwrite_external, &
      & iflag_dip2, &
      & iflag_intelectron, &
      & num_dip2, &
      & dip2boundary, &
      & dip2center, &
      & iflag_fourier_omega, &
      & num_fourier_omega, &
      & fourier_omega, &
      & itotntime2, &
      & iwdenoption, &
      & iwdenstep, &
      & iflag_estatic


!! == default for &unit ==
    unit_system='au'
!! =======================

    if (comm_is_root(nproc_id_global)) then
      open(fh_namelist, file='.namelist.tmp', status='old')
      read(fh_namelist, nml=units, iostat=inml_units)
      rewind(fh_namelist)
      close(fh_namelist)
    end if

    call comm_bcast(unit_system,nproc_group_global)
    
    select case(unit_system)
    case('au','a.u.','A_eV_fs')
      continue
    case default
      stop 'invalid unit_system'
    end select

    select case(unit_system)
    case('au','a.u.')
      unit_time='au'
      unit_length='au'
      unit_energy='au'
      unit_charge='au'
    case('A_eV_fs')
      unit_time='fs'
      unit_length='AA'
      unit_energy='eV'
      unit_charge='au'
    end select

    call initialize_inputoutput_units

!! == default for &calculation 
    theory           = 'TDDFT'
    calc_mode        = 'none'
    use_ehrenfest_md = 'n'
    use_adiabatic_md = 'n'
    use_ms_maxwell   = 'n'
    use_force        = 'n'
    use_geometry_opt = 'n'
!! == default for &control
    restart_option   = 'new'
    backup_frequency = 0
    time_shutdown    = -1d0
    sysname          = 'default'
    directory        = './'
    dump_filename    = 'default'
    modify_gs_wfn_k  = 'n'
    read_gs_wfn_k    = 'n'
    read_rt_wfn_k    = 'n'
    write_gs_wfn_k   = 'n'
    write_rt_wfn_k   = 'n'
    read_gs_wfn_k_ms = 'n'
    read_rt_wfn_k_ms = 'n'
    write_gs_wfn_k_ms= 'n'
    write_rt_wfn_k_ms= 'n'

!! == default for &parallel
    domain_parallel   = 'n'
    nproc_k           = 0
    nproc_ob          = 0
    nproc_domain      = 0
    nproc_domain_s    = 0
    num_datafiles_in  = 1
    num_datafiles_out = 1
!! == default for &system
    iperiodic          = 0
    ispin              = 0
    al                 = 0d0
    isym               = 1
    crystal_structure  = 'none'
    nstate             = 0
    nstate_spin(:)     = 0
    nelec              = 0
    nelec_spin (:)     = 0
    temperature        = -1d0
    temperature_k      = -1d0
    nelem              = 0
    natom              = 0
    file_atom_coor          = 'none'
    file_atom_red_coor          = 'none'
!! == default for &pseudo
    pseudo_file     = 'none'
    lmax_ps       = -1
    lloc_ps       = -1
    izatom        = -1
    psmask_option = 'n'
    alpha_mask    = 0.8d0
    gamma_mask    = 1.8d0
    eta_mask      = 15d0
!! == default for &functional
    xc   = 'none'
    ! xcname = 'PZ'
    xname = 'none'
    cname = 'none'
    alibx = 'none'
    alibc = 'none'
    alibxc = 'none'
    cval = -1d0
    no_update_func = 'n'
!! == default for &rgrid
    dl        = 0d0
    num_rgrid = 0
!! == default for &kgrid
    num_kgrid = 1
    file_kw   =  'none'
!! == default for &tgrid
    nt = 0
    dt = 0
!! == default for &propagation
    n_hamil = 4
    propagator = 'middlepoint'
!! == default for &scf
    amin_routine  = 'cg'
    ncg           = 5
    amixing       = 'broyden'
    rmixrate      = 0.5d0
    nmemory_mb    = 8
    alpha_mb      = 0.75d0
    fsset_option  = 'n'
    nfsset_start  = 75
    nfsset_every  = 25
    nscf          = 0
    ngeometry_opt = 1
    subspace_diagonalization = 'y'
    convergence   = 'rho_dne'
    threshold     = 1d-17
    threshold_norm_rho = -1d0/ulength_from_au**6                     ! a.u., -1 a.u. = -45.54 AA**(-6)
    threshold_norm_pot = -1d0/ulength_from_au**6*uenergy_from_au**2  ! a.u., -1 a.u. = -33.72d4 AA**(-6)*eV**2
    omp_loop      = 'k'
    skip_gsortho  = 'n'
    iditer_notemperature = 10
    gscg          = 'y'

!! == default for &emfield
    trans_longi    = 'tr'
    ae_shape1      = 'none'
    e_impulse = 1d-2*uenergy_from_au/ulength_from_au*utime_from_au ! a.u.
    amplitude1     = 0d0
    rlaser_int_wcm2_1 = -1d0
    pulse_tw1      = 0d0
    omega1         = 0d0
    epdir_re1      = (/1d0,0d0,0d0/)
    epdir_im1      = 0d0
    phi_cep1       = 0d0
    ae_shape2      = 'none'
    amplitude2     = 0d0
    rlaser_int_wcm2_2 = -1d0
    pulse_tw2      = 0d0
    omega2         = 0d0
    epdir_re2      = (/1d0,0d0,0d0/)
    epdir_im2      = 0d0
    phi_cep2       = 0d0
    t1_t2          = 0d0
    t1_delay       = 0d0
    quadrupole     = 'n'
    quadrupole_pot = ''
    alocal_laser    = 'n'
    rlaserbound_sta(1) = -1.d7*ulength_from_au ! a.u.
    rlaserbound_sta(2) = -1.d7*ulength_from_au ! a.u.
    rlaserbound_sta(3) = -1.d7*ulength_from_au ! a.u.
    rlaserbound_end(1) =  1.d7*ulength_from_au ! a.u.
    rlaserbound_end(2) =  1.d7*ulength_from_au ! a.u.
    rlaserbound_end(3) =  1.d7*ulength_from_au ! a.u.
    nump       = 0
    vecp       = 0d0
    coop       = 0d0
    radp_diele = 2d0 ! a.u.
    
!! == default for &multiscale
    fdtddim    = '1d'
    twod_shape = 'periodic'
    nx_m       = 1
    ny_m       = 1
    nz_m       = 1
    hx_m       = 0d0
    hy_m       = 0d0
    hz_m       = 0d0
    nksplit    = 0
    nxysplit   = 0
    nxvacl_m   = 1
    nxvacr_m   = 0
    nx_origin_m = 1
    ny_origin_m = 1
    nz_origin_m = 1
    file_macropoint = ''
    set_ini_coor_vel= 'n'
    nmacro_write_group=-1

!! == default for &maxwell
    al_em(:)        = 0d0
    dl_em(:)        = 0d0
    dt_em           = 0d0
    nt_em           = 0
    iboundary(:,:)  = 0
    wave_input      = 'none'
    ek_dir1(:)      = 0d0
    source_loc1(:)  = 0d0
    ek_dir2(:)      = 0d0
    source_loc2(:)  = 0d0
    iobs_num_em     = 0
    iobs_samp_em    = 1
    obs_loc_em(:,:) = 0d0
    shape_file      = 'none'
    imedia_num      = 0
    type_media(:)   = 'vacuum'
    epsilon(:)      = 1d0
    rmu(:)          = 1d0
    sigma(:)        = 0d0
    omega_p_d(:)    = 0d0
    gamma_d(:)      = 0d0
    smooth_d        = 'n'
    weight_d        = 0.5d0
    wf_em           = 'y'

!! == default for &analysis
    projection_option   = 'no'
    projection_decomp   = 'n'
    nenergy             = 1000
    de                  = (0.01d0/au_energy_ev)*uenergy_from_au  ! eV
    out_psi             = 'n'
    out_dos             = 'n'
    out_dos_start       = -1.d10 / au_energy_ev * uenergy_from_au
    out_dos_end         = +1.d10 / au_energy_ev * uenergy_from_au
    iout_dos_nenergy    = 601
    out_dos_smearing    = 0.1d0 / au_energy_ev * uenergy_from_au
    out_dos_method      = 'gaussian'
    out_dos_fshift      = 'n'
    out_pdos            = 'n'
    out_dns             = 'n'
    out_old_dns         = 'n'
    out_dns_rt          = 'n'
    out_dns_rt_step     = 50
    out_dns_trans       = 'n'
    out_dns_trans_energy= 1.55d0 / au_energy_ev * uenergy_from_au  ! eV

    out_elf             = 'n'
    out_elf_rt          = 'n'
    out_elf_rt_step     = 50
    out_estatic_rt      = 'n'
    out_estatic_rt_step = 50
    out_rvf_rt          = 'n'
    out_rvf_rt_step     = 10
    out_tm              = 'n'
    out_projection_step = 100
    out_ms_step      = 100
    format3d            = 'cube'
    numfiles_out_3d     = 1
    timer_process       = 'n'

!! == default for &hartree
    meo          = 3
    num_pole_xyz = 0
!! == default for &ewald
    newald = 4
    aewald = 0.5d0
!! == default for &opt
    cg_alpha_ini    =  0.8d0 !not use now
    cg_alpha_up     =  1.3d0
    cg_alpha_down   =  0.5d0
    convrg_scf_force= -1d0
    convrg_scf_ene  = -1d0
    convrg_opt_fmax =  1d-3
    convrg_opt_ene  =  1d-6  !not use now
!! == default for &md
    ensemble              = 'nve'
    thermostat            = 'nose-hoover'
    step_velocity_scaling = -1
    step_update_ps        = 10
    step_update_ps2       = 1
    temperature0_ion      = 298.15d0
    set_ini_velocity      = 'n'
    file_ini_velocity     = 'none'
    file_set_shake        = 'none'
    thermostat_tau        =  41.34d0  !=1fs: just test value
    friction              =  0d0
    stop_system_momt      = 'n'
!! == default for &misc
    fourier               = 'ft'
!! == default for &group_fundamental
    iditerybcg             = 20
    iditer_nosubspace_diag = 10
    ntmg                   = 1
    idisnum                = (/1,2/)
    iwrite_projection      = 0
    itwproj                = -1
    iwrite_projnum         = 0
    itcalc_ene             = 1
!! == default for &group_parallel
    isequential    = 2
    imesh_s_all    = 1
    iflag_comm_rho = 1
!! == default for &group_hartree
    hconv    = 1.d-15*uenergy_from_au**2*ulength_from_au**3 ! a.u., 1.d-15 a.u. = ! 1.10d-13 eV**2*AA**3
    lmax_meo = 4
!! == default for &group_file
    ic    = 0
    oc    = 1
    ic_rt = 0
    oc_rt = 0
!! == default for &group_others
    iparaway_ob = 2
    iscf_order  = 1
    iswitch_orbital_mesh = 0
    iflag_psicube        = 0
    lambda1_diis         = 0.5d0
    lambda2_diis         = 0.3d0
    file_ini             = 'file_ini'
    num_projection       = 1
    do ii=1,200
      iwrite_projection_ob(ii) = ii
    end do
    iwrite_projection_k(1:200) = 1
    filename_pot               = 'pot'
    iwrite_external            = 0
    iflag_dip2                 = 0
    iflag_intelectron          = 0
    num_dip2                   = 1
    dip2boundary(1:100)        = 0.d0*ulength_from_au ! a.u.
    dip2center(1:100)          = 0.d0*ulength_from_au ! a.u.
    iflag_fourier_omega        = 0
    num_fourier_omega          = 1
    fourier_omega(1:200)       = 0.d0*uenergy_from_au ! a.u.
    itotntime2                 = 0
    iwdenoption                = 0
    iwdenstep                  = 0
    iflag_estatic              = 0


    if (comm_is_root(nproc_id_global)) then
      fh_namelist = get_filehandle()
      open(fh_namelist, file='.namelist.tmp', status='old')

      read(fh_namelist, nml=calculation, iostat=inml_calculation)
      rewind(fh_namelist)

      read(fh_namelist, nml=control, iostat=inml_control)
      rewind(fh_namelist)

      read(fh_namelist, nml=parallel, iostat=inml_parallel)
      rewind(fh_namelist)

      read(fh_namelist, nml=system, iostat=inml_system)
      rewind(fh_namelist)

      read(fh_namelist, nml=pseudo, iostat=inml_pseudo)
      rewind(fh_namelist)

      read(fh_namelist, nml=functional, iostat=inml_functional)
      rewind(fh_namelist)

      read(fh_namelist, nml=rgrid, iostat=inml_rgrid)
      rewind(fh_namelist)

      read(fh_namelist, nml=kgrid, iostat=inml_kgrid)
      rewind(fh_namelist)

      read(fh_namelist, nml=tgrid, iostat=inml_tgrid)
      rewind(fh_namelist)

      read(fh_namelist, nml=propagation, iostat=inml_propagation)
      rewind(fh_namelist)

      read(fh_namelist, nml=scf, iostat=inml_scf)
      rewind(fh_namelist)

      read(fh_namelist, nml=emfield, iostat=inml_emfield)
      rewind(fh_namelist)

      read(fh_namelist, nml=multiscale, iostat=inml_multiscale)
      rewind(fh_namelist)

      read(fh_namelist, nml=maxwell, iostat=inml_maxwell)
      rewind(fh_namelist)

      read(fh_namelist, nml=analysis, iostat=inml_analysis)
      rewind(fh_namelist)

      read(fh_namelist, nml=hartree, iostat=inml_hartree)
      rewind(fh_namelist)

      read(fh_namelist, nml=ewald, iostat=inml_ewald)
      rewind(fh_namelist)

      read(fh_namelist, nml=opt, iostat=inml_opt)
      rewind(fh_namelist)

      read(fh_namelist, nml=md, iostat=inml_md)
      rewind(fh_namelist)

      read(fh_namelist, nml=misc, iostat=inml_misc)
      rewind(fh_namelist)

      read(fh_namelist, nml=group_fundamental, iostat=inml_group_fundamental)
      rewind(fh_namelist)

      read(fh_namelist, nml=group_parallel, iostat=inml_group_parallel)
      rewind(fh_namelist)

      read(fh_namelist, nml=group_hartree, iostat=inml_group_hartree)
      rewind(fh_namelist)

      read(fh_namelist, nml=group_file, iostat=inml_group_file)
      rewind(fh_namelist)

      read(fh_namelist, nml=group_others, iostat=inml_group_others)
      rewind(fh_namelist)

      close(fh_namelist)
    end if

! Broad cast
!! == bcast for &calculation
    call comm_bcast(theory          ,nproc_group_global)
    call comm_bcast(calc_mode       ,nproc_group_global)
    call comm_bcast(use_ehrenfest_md,nproc_group_global)
    call comm_bcast(use_adiabatic_md,nproc_group_global)
    call comm_bcast(use_ms_maxwell  ,nproc_group_global)
    call comm_bcast(use_force       ,nproc_group_global)
    call comm_bcast(use_geometry_opt,nproc_group_global)
!! == bcast for &control
    call comm_bcast(restart_option  ,nproc_group_global)
    call comm_bcast(backup_frequency,nproc_group_global)
    call comm_bcast(time_shutdown   ,nproc_group_global)
    call comm_bcast(sysname         ,nproc_group_global)
    call comm_bcast(directory       ,nproc_group_global)
    if(directory(len_trim(directory):len_trim(directory)).ne.'/') &
    &  directory = trim(directory)//'/'
    call comm_bcast(dump_filename   ,nproc_group_global)
    call comm_bcast(modify_gs_wfn_k ,nproc_group_global)
    call comm_bcast(read_gs_wfn_k   ,nproc_group_global)
    call comm_bcast(read_rt_wfn_k   ,nproc_group_global)
    call comm_bcast(write_gs_wfn_k  ,nproc_group_global)
    call comm_bcast(write_rt_wfn_k  ,nproc_group_global)
    call comm_bcast(read_gs_wfn_k_ms ,nproc_group_global)
    call comm_bcast(read_rt_wfn_k_ms ,nproc_group_global)
    call comm_bcast(write_gs_wfn_k_ms,nproc_group_global)
    call comm_bcast(write_rt_wfn_k_ms,nproc_group_global)


!! == bcast for &parallel
    call comm_bcast(domain_parallel  ,nproc_group_global)
    call comm_bcast(nproc_k          ,nproc_group_global)
    call comm_bcast(nproc_ob         ,nproc_group_global)
    call comm_bcast(nproc_domain     ,nproc_group_global)
    call comm_bcast(nproc_domain_s   ,nproc_group_global)
    call comm_bcast(num_datafiles_in ,nproc_group_global)
    call comm_bcast(num_datafiles_out,nproc_group_global)
!! == bcast for &system
    call comm_bcast(iperiodic,nproc_group_global)
    call comm_bcast(ispin    ,nproc_group_global)
    call comm_bcast(al       ,nproc_group_global)
    al = al * ulength_to_au
    call comm_bcast(isym               ,nproc_group_global)
    call comm_bcast(crystal_structure  ,nproc_group_global)
    call comm_bcast(nstate             ,nproc_group_global)
    call comm_bcast(nstate_spin        ,nproc_group_global)
    call comm_bcast(nelec              ,nproc_group_global)
    call comm_bcast(nelec_spin         ,nproc_group_global)
    call comm_bcast(temperature        ,nproc_group_global)
    call comm_bcast(temperature_k      ,nproc_group_global)
    call comm_bcast(nelem              ,nproc_group_global)
    call comm_bcast(natom              ,nproc_group_global)
    call comm_bcast(file_atom_coor     ,nproc_group_global)
    call comm_bcast(file_atom_red_coor ,nproc_group_global)
!! == bcast for &pseudo
    call comm_bcast(pseudo_file  ,nproc_group_global)
    call comm_bcast(lmax_ps      ,nproc_group_global)
    call comm_bcast(lloc_ps      ,nproc_group_global)
    call comm_bcast(izatom       ,nproc_group_global)
    call comm_bcast(psmask_option,nproc_group_global)
    call comm_bcast(alpha_mask   ,nproc_group_global)
    call comm_bcast(gamma_mask   ,nproc_group_global)
    call comm_bcast(eta_mask     ,nproc_group_global)
!! == bcast for &functional

#ifdef SALMON_USE_LIBXC
    if (alibxc .ne. 'none') xc = 'libxc:' // trim(alibxc)
    if (alibx .ne. 'none') xname = 'libxc:' // trim(alibx)
    if (alibc .ne. 'none') cname = 'libxc:' // trim(alibc)
#endif
    call comm_bcast(xc           ,nproc_group_global)
    call comm_bcast(cname        ,nproc_group_global)
    call comm_bcast(xname        ,nproc_group_global)
#ifdef SALMON_USE_LIBXC
    call comm_bcast(alibxc       ,nproc_group_global)
    call comm_bcast(alibx        ,nproc_group_global)
    call comm_bcast(alibc        ,nproc_group_global)
#endif
    call comm_bcast(cval         ,nproc_group_global)
    call comm_bcast(no_update_func,nproc_group_global)
!! == bcast for &rgrid
    call comm_bcast(dl,nproc_group_global)
    dl = dl * ulength_to_au
    call comm_bcast(num_rgrid,nproc_group_global)
!! == bcast for &kgrid
    call comm_bcast(num_kgrid,nproc_group_global)
    call comm_bcast(file_kw  ,nproc_group_global)
!! == bcast for &tgrid
    call comm_bcast(nt,nproc_group_global)
    call comm_bcast(dt,nproc_group_global)
    dt = dt * utime_to_au
!! == bcast for &propagation
    call comm_bcast(n_hamil   ,nproc_group_global)
    call comm_bcast(propagator,nproc_group_global)
!! == bcast for &scf
    call comm_bcast(amin_routine            ,nproc_group_global)
    call comm_bcast(ncg                     ,nproc_group_global)
    call comm_bcast(amixing                 ,nproc_group_global)
    call comm_bcast(rmixrate                ,nproc_group_global)
    call comm_bcast(nmemory_mb              ,nproc_group_global)
    call comm_bcast(alpha_mb                ,nproc_group_global)
    call comm_bcast(fsset_option            ,nproc_group_global)
    call comm_bcast(nfsset_start            ,nproc_group_global)
    call comm_bcast(nfsset_every            ,nproc_group_global)
    call comm_bcast(nscf                    ,nproc_group_global)
    call comm_bcast(ngeometry_opt           ,nproc_group_global)
    call comm_bcast(subspace_diagonalization,nproc_group_global)
    call comm_bcast(convergence             ,nproc_group_global)
    call comm_bcast(threshold               ,nproc_group_global)
    call comm_bcast(threshold_norm_rho      ,nproc_group_global)
    threshold_norm_rho = threshold_norm_rho / (ulength_to_au)**6
    call comm_bcast(threshold_norm_pot      ,nproc_group_global)
    threshold_norm_pot = threshold_norm_pot * (uenergy_to_au)**2 / (ulength_to_au)**6
    call comm_bcast(omp_loop                ,nproc_group_global)
    call comm_bcast(skip_gsortho            ,nproc_group_global)
    call comm_bcast(iditer_notemperature    ,nproc_group_global)
    call comm_bcast(gscg                    ,nproc_group_global)
!! == bcast for &emfield
    call comm_bcast(trans_longi,nproc_group_global)
    call comm_bcast(ae_shape1  ,nproc_group_global)
    call comm_bcast(e_impulse,nproc_group_global)
    e_impulse = e_impulse *uenergy_to_au/ulength_to_au*utime_to_au
    call comm_bcast(amplitude1 ,nproc_group_global)
    amplitude1 = amplitude1*(uenergy_to_au/ulength_to_au/ucharge_to_au)
    call comm_bcast(rlaser_int_wcm2_1,nproc_group_global)
    call comm_bcast(pulse_tw1  ,nproc_group_global)
    pulse_tw1 = pulse_tw1 * utime_to_au
    call comm_bcast(omega1,nproc_group_global)
    omega1 = omega1 * uenergy_to_au
    call comm_bcast(epdir_re1 ,nproc_group_global)
    call comm_bcast(epdir_im1 ,nproc_group_global)
    norm = sqrt(sum(epdir_re1(:)**2+epdir_im1(:)**2))
    epdir_re1 = epdir_re1 / norm
    epdir_im1 = epdir_im1 / norm
    call comm_bcast(phi_cep1  ,nproc_group_global)
    call comm_bcast(ae_shape2 ,nproc_group_global)
    call comm_bcast(amplitude2,nproc_group_global)
    amplitude2 = amplitude2*(uenergy_to_au/ulength_to_au/ucharge_to_au)
    call comm_bcast(rlaser_int_wcm2_2,nproc_group_global)
    call comm_bcast(pulse_tw2  ,nproc_group_global)
    pulse_tw2 = pulse_tw2 * utime_to_au
    call comm_bcast(omega2,nproc_group_global)
    omega2 = omega2 * uenergy_to_au
    call comm_bcast(epdir_re2,nproc_group_global)
    call comm_bcast(epdir_im2,nproc_group_global)
    norm = sqrt(sum(epdir_re2(:)**2+epdir_im2(:)**2))
    epdir_re2 = epdir_re2 / norm
    epdir_im2 = epdir_im2 / norm
    call comm_bcast(phi_cep2 ,nproc_group_global)
    call comm_bcast(t1_t2    ,nproc_group_global)
    t1_t2 = t1_t2 * utime_to_au
    call comm_bcast(t1_delay ,nproc_group_global)
    t1_delay = t1_delay * utime_to_au
    call comm_bcast(quadrupole    ,nproc_group_global)
    call comm_bcast(quadrupole_pot,nproc_group_global)
    call comm_bcast(alocal_laser  ,nproc_group_global)
    call comm_bcast(rlaserbound_sta,nproc_group_global)
    rlaserbound_sta = rlaserbound_sta * ulength_to_au
    call comm_bcast(rlaserbound_end,nproc_group_global)
    rlaserbound_end = rlaserbound_end * ulength_to_au
    call comm_bcast(nump,nproc_group_global)
    call comm_bcast(vecp,nproc_group_global)
    vecp = vecp * ulength_to_au
    call comm_bcast(coop,nproc_group_global)
    coop = coop * ulength_to_au
    call comm_bcast(radp_diele,nproc_group_global)
    radp_diele = radp_diele * ulength_to_au
    
!! == bcast for &multiscale
    call comm_bcast(fdtddim   ,nproc_group_global)
    call comm_bcast(twod_shape,nproc_group_global)
    call comm_bcast(nx_m      ,nproc_group_global)
    call comm_bcast(ny_m      ,nproc_group_global)
    call comm_bcast(nz_m      ,nproc_group_global)
    call comm_bcast(hx_m      ,nproc_group_global)
    hx_m = hx_m * ulength_to_au
    call comm_bcast(hy_m      ,nproc_group_global)
    hy_m = hy_m * ulength_to_au
    call comm_bcast(hz_m      ,nproc_group_global)
    hz_m = hz_m * ulength_to_au
    call comm_bcast(nksplit   ,nproc_group_global)
    call comm_bcast(nxysplit  ,nproc_group_global)
    call comm_bcast(nxvacl_m  ,nproc_group_global)
    call comm_bcast(nxvacr_m  ,nproc_group_global)
    call comm_bcast(nx_origin_m,nproc_group_global)
    call comm_bcast(ny_origin_m,nproc_group_global)
    call comm_bcast(nz_origin_m,nproc_group_global)
    call comm_bcast(file_macropoint, nproc_group_global)
    call comm_bcast(num_macropoint,  nproc_group_global)
    call comm_bcast(set_ini_coor_vel,nproc_group_global)
    call comm_bcast(nmacro_write_group,nproc_group_global)

!! == bcast for &maxwell
    call comm_bcast(al_em        ,nproc_group_global)
    al_em = al_em * ulength_to_au
    call comm_bcast(dl_em        ,nproc_group_global)
    dl_em = dl_em * ulength_to_au
    call comm_bcast(dt_em        ,nproc_group_global)
    dt_em = dt_em * utime_to_au
    call comm_bcast(nt_em        ,nproc_group_global)
    call comm_bcast(iboundary    ,nproc_group_global)
    call comm_bcast(wave_input,nproc_group_global)
    call comm_bcast(ek_dir1      ,nproc_group_global)
    call comm_bcast(source_loc1  ,nproc_group_global)
    source_loc1 = source_loc1 * ulength_to_au
    call comm_bcast(ek_dir2      ,nproc_group_global)
    call comm_bcast(source_loc2  ,nproc_group_global)
    source_loc2 = source_loc2 * ulength_to_au
    call comm_bcast(iobs_num_em  ,nproc_group_global)
    call comm_bcast(iobs_samp_em ,nproc_group_global)
    call comm_bcast(obs_loc_em   ,nproc_group_global)
    obs_loc_em = obs_loc_em * ulength_to_au
    call comm_bcast(shape_file   ,nproc_group_global)
    call comm_bcast(imedia_num   ,nproc_group_global)
    call comm_bcast(type_media   ,nproc_group_global)
    call comm_bcast(epsilon      ,nproc_group_global)
    call comm_bcast(rmu          ,nproc_group_global)
    call comm_bcast(sigma        ,nproc_group_global)
    call comm_bcast(omega_p_d    ,nproc_group_global)
    omega_p_d = omega_p_d * uenergy_to_au
    call comm_bcast(gamma_d      ,nproc_group_global)
    gamma_d = gamma_d * uenergy_to_au
    call comm_bcast(smooth_d     ,nproc_group_global)
    call comm_bcast(weight_d     ,nproc_group_global)
    call comm_bcast(wf_em        ,nproc_group_global)
    
!! == bcast for &analysis
    call comm_bcast(projection_option,nproc_group_global)
    call comm_bcast(projection_decomp,nproc_group_global)
    call comm_bcast(nenergy          ,nproc_group_global)
    call comm_bcast(de               ,nproc_group_global)
    de = de * uenergy_to_au
    call comm_bcast(out_psi             ,nproc_group_global)
    call comm_bcast(out_dos             ,nproc_group_global)
    call comm_bcast(out_dos_start       ,nproc_group_global)
    out_dos_start = out_dos_start * uenergy_to_au
    call comm_bcast(out_dos_end         ,nproc_group_global)
    out_dos_end = out_dos_end * uenergy_to_au
    call comm_bcast(iout_dos_nenergy    ,nproc_group_global)
    call comm_bcast(out_dos_smearing    ,nproc_group_global)
    out_dos_smearing = out_dos_smearing * uenergy_to_au
    call comm_bcast(out_dos_method      ,nproc_group_global)
    call comm_bcast(out_dos_fshift      ,nproc_group_global)
    call comm_bcast(out_pdos            ,nproc_group_global)
    call comm_bcast(out_dns             ,nproc_group_global)
    call comm_bcast(out_old_dns         ,nproc_group_global)
    call comm_bcast(out_dns_rt          ,nproc_group_global)
    call comm_bcast(out_dns_rt_step     ,nproc_group_global)
    call comm_bcast(out_dns_trans       ,nproc_group_global)
    call comm_bcast(out_dns_trans_energy,nproc_group_global)
    out_dns_trans_energy = out_dns_trans_energy * uenergy_to_au
    call comm_bcast(out_elf            ,nproc_group_global)
    call comm_bcast(out_elf_rt         ,nproc_group_global)
    call comm_bcast(out_elf_rt_step    ,nproc_group_global)
    call comm_bcast(out_estatic_rt     ,nproc_group_global)
    call comm_bcast(out_estatic_rt_step,nproc_group_global)
    call comm_bcast(out_rvf_rt         ,nproc_group_global)
    call comm_bcast(out_rvf_rt_step    ,nproc_group_global)
    call comm_bcast(out_tm             ,nproc_group_global)
    call comm_bcast(out_projection_step,nproc_group_global)
    call comm_bcast(out_ms_step     ,nproc_group_global)
    call comm_bcast(format3d           ,nproc_group_global)
    call comm_bcast(numfiles_out_3d    ,nproc_group_global)
    call comm_bcast(timer_process      ,nproc_group_global)

!! == bcast for &hartree
    call comm_bcast(meo         ,nproc_group_global)
    call comm_bcast(num_pole_xyz,nproc_group_global)
!! == bcast for &ewald
    call comm_bcast(newald,nproc_group_global)
    call comm_bcast(aewald,nproc_group_global)
!! == bcast for &opt
    call comm_bcast(cg_alpha_ini     ,nproc_group_global)
    call comm_bcast(cg_alpha_up      ,nproc_group_global)
    call comm_bcast(cg_alpha_down    ,nproc_group_global)
    call comm_bcast(convrg_scf_force ,nproc_group_global)
    call comm_bcast(convrg_scf_ene   ,nproc_group_global)
    call comm_bcast(convrg_opt_fmax  ,nproc_group_global)
    call comm_bcast(convrg_opt_ene   ,nproc_group_global)
!! == bcast for &md
    call comm_bcast(ensemble               ,nproc_group_global)
    call comm_bcast(thermostat             ,nproc_group_global)
    call comm_bcast(step_velocity_scaling  ,nproc_group_global)
    call comm_bcast(step_update_ps         ,nproc_group_global)
    call comm_bcast(step_update_ps2        ,nproc_group_global)
    call comm_bcast(temperature0_ion       ,nproc_group_global)
    call comm_bcast(set_ini_velocity       ,nproc_group_global)
    call comm_bcast(file_ini_velocity      ,nproc_group_global)
    call comm_bcast(file_set_shake         ,nproc_group_global)
    call comm_bcast(thermostat_tau         ,nproc_group_global)
    thermostat_tau = thermostat_tau * utime_to_au
    call comm_bcast(friction               ,nproc_group_global)
    call comm_bcast(stop_system_momt       ,nproc_group_global)
!! == bcast for &misc
    call comm_bcast(fourier                ,nproc_group_global)
!! == bcast for &group_fundamental
    call comm_bcast(iditerybcg            ,nproc_group_global)
    call comm_bcast(iditer_nosubspace_diag,nproc_group_global)
    call comm_bcast(ntmg                  ,nproc_group_global)
    call comm_bcast(idisnum               ,nproc_group_global)
    call comm_bcast(iwrite_projection     ,nproc_group_global)
    call comm_bcast(itwproj               ,nproc_group_global)
    call comm_bcast(iwrite_projnum        ,nproc_group_global)
    call comm_bcast(itcalc_ene            ,nproc_group_global)
!! == bcast for &group_parallel
    call comm_bcast(isequential   ,nproc_group_global)
    call comm_bcast(imesh_s_all   ,nproc_group_global)
    call comm_bcast(iflag_comm_rho,nproc_group_global)
!! == bcast for &group_hartree
    call comm_bcast(hconv   ,nproc_group_global)
    hconv = hconv * (uenergy_to_au)**2 * (ulength_to_au)**3 
    call comm_bcast(lmax_meo,nproc_group_global)
!! == bcast for &group_file
    call comm_bcast(ic   ,nproc_group_global)
    call comm_bcast(oc   ,nproc_group_global)
    call comm_bcast(ic_rt,nproc_group_global)
    call comm_bcast(oc_rt,nproc_group_global)
!! == bcast for &group_others
    call comm_bcast(iparaway_ob         ,nproc_group_global)
    call comm_bcast(iscf_order          ,nproc_group_global)
    call comm_bcast(iswitch_orbital_mesh,nproc_group_global)
    call comm_bcast(iflag_psicube       ,nproc_group_global)
    call comm_bcast(lambda1_diis        ,nproc_group_global)
    call comm_bcast(lambda2_diis        ,nproc_group_global)
    call comm_bcast(file_ini            ,nproc_group_global)
    call comm_bcast(iparaway_ob         ,nproc_group_global)
    call comm_bcast(num_projection      ,nproc_group_global)
    call comm_bcast(iwrite_projection_ob,nproc_group_global)
    call comm_bcast(iwrite_projection_k ,nproc_group_global)
    call comm_bcast(filename_pot        ,nproc_group_global)
    call comm_bcast(iwrite_external     ,nproc_group_global)
    call comm_bcast(iflag_dip2          ,nproc_group_global)
    call comm_bcast(iflag_intelectron   ,nproc_group_global)
    call comm_bcast(num_dip2            ,nproc_group_global)
    call comm_bcast(dip2boundary        ,nproc_group_global)
    dip2boundary  = dip2boundary  * ulength_to_au
    call comm_bcast(dip2center          ,nproc_group_global)
    dip2center    = dip2center    * ulength_to_au
    call comm_bcast(iflag_fourier_omega ,nproc_group_global)
    call comm_bcast(num_fourier_omega   ,nproc_group_global)
    call comm_bcast(fourier_omega       ,nproc_group_global)
    fourier_omega = fourier_omega * uenergy_to_au
    call comm_bcast(itotntime2          ,nproc_group_global)
    call comm_bcast(iwdenoption         ,nproc_group_global)
    call comm_bcast(iwdenstep           ,nproc_group_global)
    call comm_bcast(iflag_estatic       ,nproc_group_global)

  end subroutine read_input_common

  subroutine read_atomic_coordinates
    use salmon_parallel
    use salmon_communication
    use salmon_file, only: get_filehandle
    character(256) :: filename_tmp,char_atom
    integer :: icount,i
    logical :: if_error, if_cartesian


    if (comm_is_root(nproc_id_global)) then


      if_error = .false.
      if_cartesian = .true.
      iflag_atom_coor = ntype_atom_coor_none
      icount = 0
      if(file_atom_coor /= 'none')then
        icount = icount + 1
        if_cartesian = .true.
        filename_tmp = trim(file_atom_coor)
        iflag_atom_coor = ntype_atom_coor_cartesian
      end if

      if(file_atom_red_coor /= 'none')then
        icount = icount + 1
        if_cartesian = .false.
        filename_tmp = trim(file_atom_coor)
        iflag_atom_coor = ntype_atom_coor_reduced
      end if

      if(if_nml_coor)then
        icount = icount + 1
        if_cartesian = .true.
        filename_tmp = '.atomic_coor.tmp'
        iflag_atom_coor = ntype_atom_coor_cartesian
      end if

      if(if_nml_red_coor)then
        icount = icount + 1
        if_cartesian = .false.
        filename_tmp = '.atomic_red_coor.tmp'
        iflag_atom_coor = ntype_atom_coor_reduced
      end if

    end if


    call comm_bcast(icount,nproc_group_global)
    call comm_bcast(if_cartesian,nproc_group_global)
    call comm_bcast(iflag_atom_coor,nproc_group_global)

    if(0 < natom .and. icount/=1)then
       if (comm_is_root(nproc_id_global))then
         write(*,"(A)")'Error in input: The following inputs are incompatible.'
         write(*,"(A)")'file_atom_coor, file_atom_red_coor, &atomic_coor, and &atomic_red_coor.'
       end if
       call end_parallel
       stop
    end if

    if( (.not.if_cartesian) .and. iperiodic == 0)then
       if (comm_is_root(nproc_id_global))then
         write(*,"(A)")'Error in input: Reduced coordinate is invalid for isolated systems.'
       end if
       call end_parallel
       stop
    end if

    allocate(atom_name(natom))
    allocate(rion(3,natom), rion_red(3,natom),kion(natom), flag_geo_opt_atom(natom))
    rion = 0d0
    rion_red = 0d0
    kion = 0
    flag_geo_opt_atom = 'n'
    
    if (0 < natom) then
      
      if (comm_is_root(nproc_id_global))then
        fh_atomic_coor = get_filehandle()
        open(fh_atomic_coor, file=filename_tmp, status='old')
        select case(iflag_atom_coor)
        case(ntype_atom_coor_cartesian)
           do i=1, natom
              if(use_geometry_opt == 'y')then
                 read(fh_atomic_coor, *) char_atom, rion(:,i), kion(i), flag_geo_opt_atom(i)
              else
                 read(fh_atomic_coor, *) char_atom, rion(:,i), kion(i)
              end if
              atom_name(i) = char_atom
           end do
           rion = rion*ulength_to_au
        case(ntype_atom_coor_reduced)
           do i=1, natom
              if(use_geometry_opt == 'y')then
                 read(fh_atomic_coor, *) char_atom, rion_red(:,i), kion(i), flag_geo_opt_atom(i)
              else
                 read(fh_atomic_coor, *) char_atom, rion_red(:,i), kion(i)
              end if
              atom_name(i) = char_atom
           end do
        end select
        close(fh_atomic_coor)
        
      end if

      call comm_bcast(rion,nproc_group_global)
      call comm_bcast(rion_red,nproc_group_global)
      call comm_bcast(kion,nproc_group_global)
      call comm_bcast(flag_geo_opt_atom,nproc_group_global)
      call comm_bcast(atom_name,nproc_group_global)
    end if ! if 0 < natom


  end subroutine read_atomic_coordinates

  subroutine initialize_inputoutput_units
    implicit none


! Unit for time
    select case(unit_time)
    case('au','a.u.')
      utime_to_au   = 1d0
      utime_from_au = 1d0
      iflag_unit_time = ntype_unit_time_au
    case('fs','femtosecond')
      utime_to_au   = 1d0/au_time_fs
      utime_from_au = au_time_fs
      iflag_unit_time = ntype_unit_time_fs
    case default
      stop "Invalid unit for time."
    end select

! Unit for length
    select case(unit_length)
    case('au','a.u.')
      ulength_to_au   = 1d0
      ulength_from_au = 1d0
      iflag_unit_length = ntype_unit_length_au
    case('AA','angstrom','Angstrom')
      ulength_to_au   = 1d0/au_length_aa
      ulength_from_au = au_length_aa
      iflag_unit_length = ntype_unit_length_aa
    case default
      stop "Invalid unit for length."
    end select

! Unit for energy
    select case(unit_energy)
    case('au','a.u.')
      uenergy_to_au   = 1d0
      uenergy_from_au = 1d0
      iflag_unit_energy = ntype_unit_energy_au
    case('ev','eV')
      uenergy_to_au   = 1d0/au_energy_ev
      uenergy_from_au = au_energy_ev
      iflag_unit_energy = ntype_unit_energy_ev
    case default
      stop "Invalid unit for energy."
    end select

! Unit for charge
    select case(unit_charge)
    case('au','a.u.')
      ucharge_to_au   = 1d0
      ucharge_from_au = 1d0
      iflag_unit_charge = ntype_unit_charge_au
    case default
      stop "Invalid unit for charge."
    end select

!! prepare type(unit_t) :: t_unit_energy,t_unit_energy_inv
    t_unit_energy%conv = uenergy_from_au
    t_unit_energy_inv%conv = 1d0/uenergy_from_au
    if(iflag_unit_energy == ntype_unit_energy_ev)then
      t_unit_energy%name     = 'eV'
      t_unit_energy_inv%name = '1/eV'
    else 
      t_unit_energy%name     = 'a.u.'
      t_unit_energy_inv%name = 'a.u.'
      t_unit_energy%conv = 1d0
      t_unit_energy_inv%conv = 1d0
    end if

!! prepare type(unit_t) :: t_unit_time,t_unit_time_inv
    t_unit_time%conv = utime_from_au
    t_unit_time_inv%conv = 1d0/utime_from_au
    if(iflag_unit_time == ntype_unit_time_fs)then
      t_unit_time%name     = 'fs'
      t_unit_time_inv%name = '1/fs'
    else 
      t_unit_time%name     = 'a.u.'
      t_unit_time_inv%name = 'a.u.'
      t_unit_time%conv = 1d0
      t_unit_time_inv%conv = 1d0
    end if

!! prepare type(unit_t) :: t_unit_current
    t_unit_current%conv = (ulength_from_au/utime_from_au)/ulength_from_au**3
    if(iflag_unit_time == ntype_unit_time_fs .and. &
       iflag_unit_length == ntype_unit_length_aa &
         )then
      t_unit_current%name  = '1/fs*Angstrom^2'
    else 
      t_unit_current%name  = 'a.u.'
      t_unit_current%conv  = 1d0
    end if

!! prepare type(unit_t) :: t_unit_ac
    t_unit_ac%conv = utime_from_au*uenergy_from_au/ulength_from_au/ucharge_from_au
    if(iflag_unit_time == ntype_unit_time_fs .and. &
       iflag_unit_length == ntype_unit_length_aa .and. &
       iflag_unit_energy == ntype_unit_energy_ev .and. &
       iflag_unit_charge == ntype_unit_charge_au &
         )then
      t_unit_ac%name     = 'fs*V/Angstrom'
    else 
      t_unit_ac%name     = 'a.u.'
      t_unit_ac%conv     = 1d0
    end if

    !! prepare type(unit_t) :: t_unit_elec
    if(iflag_unit_time == ntype_unit_time_fs .and. &
       iflag_unit_length == ntype_unit_length_aa .and. &
       iflag_unit_energy == ntype_unit_energy_ev .and. &
       iflag_unit_charge == ntype_unit_charge_au &
         )then
      t_unit_elec%name     = 'V/Angstrom'
      t_unit_elec%conv     = 51.42206707d0
    else 
      t_unit_elec%name     = 'a.u.'
      t_unit_elec%conv     = 1d0
    end if


  end subroutine initialize_inputoutput_units

  subroutine dump_input_common
    use salmon_parallel
    use salmon_communication
    use salmon_file, only: get_filehandle
    use misc_routines, only: create_directory
    implicit none
    integer :: i,ierr_nml
    ierr_nml = 0

    if (comm_is_root(nproc_id_global)) then

      fh_variables_log = get_filehandle()
      open(fh_variables_log,file='variables.log')

      if(inml_calculation >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'calculation', inml_calculation
      write(fh_variables_log, '("#",4X,A,"=",A)') 'theory', theory
      write(fh_variables_log, '("#",4X,A,"=",A)') 'calc_mode', calc_mode
      write(fh_variables_log, '("#",4X,A,"=",A)') 'use_ehrenfest_md', use_ehrenfest_md
      write(fh_variables_log, '("#",4X,A,"=",A)') 'use_adiabatic_md', use_adiabatic_md
      write(fh_variables_log, '("#",4X,A,"=",A)') 'use_ms_maxwell', use_ms_maxwell
      write(fh_variables_log, '("#",4X,A,"=",A)') 'use_force', use_force
      write(fh_variables_log, '("#",4X,A,"=",A)') 'use_geometry_opt', use_geometry_opt

      if(inml_control >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'control', inml_control
      write(fh_variables_log, '("#",4X,A,"=",A)') 'restart_option', restart_option
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'backup_frequency', backup_frequency
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'time_shutdown', time_shutdown
      write(fh_variables_log, '("#",4X,A,"=",A)') 'sysname', trim(sysname)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'directory', trim(directory)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'dump_filename', trim(dump_filename)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'modify_gs_wfn_k', trim(modify_gs_wfn_k)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'read_gs_wfn_k', trim(read_gs_wfn_k)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'read_rt_wfn_k', trim(read_rt_wfn_k)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'write_gs_wfn_k', trim(write_gs_wfn_k)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'write_rt_wfn_k', trim(write_rt_wfn_k)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'read_gs_wfn_k_ms', trim(read_gs_wfn_k_ms)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'read_rt_wfn_k_ms', trim(read_rt_wfn_k_ms)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'write_gs_wfn_k_ms', trim(write_gs_wfn_k_ms)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'write_rt_wfn_k_ms', trim(write_rt_wfn_k_ms)


      if(inml_units >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'units', inml_units
      write(fh_variables_log, '("#",4X,A,"=",A)') 'unit_system', unit_system

      if(inml_parallel >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'parallel', inml_parallel
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nproc_k', nproc_k
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nproc_ob', nproc_ob
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nproc_domain(1)', nproc_domain(1)
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nproc_domain(2)', nproc_domain(2)
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nproc_domain(3)', nproc_domain(3)
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nproc_domain_s(1)', nproc_domain_s(1)
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nproc_domain_s(2)', nproc_domain_s(2)
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nproc_domain_s(3)', nproc_domain_s(3)
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'num_datafiles_in', num_datafiles_in
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'num_datafiles_out', num_datafiles_out

      if(inml_system >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'system', inml_system
      write(fh_variables_log, '("#",4X,A,"=",I1)') 'iperiodic', iperiodic
      write(fh_variables_log, '("#",4X,A,"=",I1)') 'ispin', ispin
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'al(1)', al(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'al(2)', al(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'al(3)', al(3)
      write(fh_variables_log, '("#",4X,A,"=",I1)') 'isym', isym
      write(fh_variables_log, '("#",4X,A,"=",A)') 'crystal_structure', crystal_structure
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nstate', nstate
      write(fh_variables_log, '("#",4X,A,"=",I4,2x,I4)') 'nstate_spin(1:2)', nstate_spin
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nelec', nelec
      write(fh_variables_log, '("#",4X,A,"=",I4,2x,I4)') 'nelec_spin(1:2)', nelec_spin
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'temperature', temperature
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'temperature_k', temperature_k
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nelem', nelem
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'natom', natom
      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_atom_coor', trim(file_atom_coor)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_atom_red_coor', trim(file_atom_red_coor)

      if(inml_pseudo >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'pseudo', inml_pseudo

      do i = 1,nelem
        write(fh_variables_log, '("#",4X,A,I2,A,"=",A)') 'pseudo_file(',i,')', trim(pseudo_file(i))
        write(fh_variables_log, '("#",4X,A,I2,A,"=",I4)') 'lmax_ps(',i,')', lmax_ps(i)
        write(fh_variables_log, '("#",4X,A,I2,A,"=",I4)') 'lloc_ps(',i,')', lloc_ps(i)
        write(fh_variables_log, '("#",4X,A,I2,A,"=",I4)') 'izatom(',i,')', izatom(i)
      end do
      write(fh_variables_log, '("#",4X,A,"=",A)') 'psmask_option', psmask_option
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'alpha_mask', alpha_mask
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'gamma_mask', gamma_mask
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'eta_mask', eta_mask

      if(inml_functional >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'functional', inml_functional
      write(fh_variables_log, '("#",4X,A,"=",A)') 'xc', trim(xc)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'xname', trim(xname)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'cname', trim(cname)
#ifdef SALMON_USE_LIBXC
      write(fh_variables_log, '("#",4X,A,"=",A)') 'alibxc', trim(alibxc)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'alibx', trim(alibx)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'alibc', trim(alibc)
#endif
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cval', cval
      write(fh_variables_log, '("#",4X,A,"=",A)') 'no_update_func', no_update_func

      if(inml_rgrid >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'rgrid', inml_rgrid
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl(1)', dl(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl(2)', dl(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl(3)', dl(3)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_rgrid(1)', num_rgrid(1)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_rgrid(2)', num_rgrid(2)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_rgrid(3)', num_rgrid(3)

      if(inml_kgrid >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'kgrid', inml_kgrid
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_kgrid(1)', num_kgrid(1)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_kgrid(2)', num_kgrid(2)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_kgrid(3)', num_kgrid(3)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_kw', trim(file_kw)

      if(inml_tgrid >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'tgrid', inml_tgrid
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'nt', nt
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dt', dt

      if(inml_scf >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#",4X,A,"=",A)') 'amin_routine', amin_routine
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'scf', inml_scf
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'ncg', ncg
      write(fh_variables_log, '("#",4X,A,"=",A)') 'amixing', amixing
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rmixrate', rmixrate
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'nmemory_mb', nmemory_mb
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'alpha_mb', alpha_mb
      write(fh_variables_log, '("#",4X,A,"=",A)') 'fsset_option', fsset_option
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'nfsset_start', nfsset_start
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'nfsset_every', nfsset_every
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'nscf', nscf
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'ngeometry_opt', ngeometry_opt
      write(fh_variables_log, '("#",4X,A,"=",A)') 'subspace_diagonalization', subspace_diagonalization
      write(fh_variables_log, '("#",4X,A,"=",A)') 'convergence', convergence
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'threshold', threshold
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'threshold_norm_rho', threshold_norm_rho
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'threshold_norm_pot', threshold_norm_pot
      write(fh_variables_log, '("#",4X,A,"=",A)') 'omp_loop', omp_loop
      write(fh_variables_log, '("#",4X,A,"=",A)') 'skip_gsortho', skip_gsortho
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'iditer_notemperature', iditer_notemperature
      write(fh_variables_log, '("#",4X,A,"=",A)') 'gscg', gscg

      if(inml_emfield >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'emfield', inml_emfield
      write(fh_variables_log, '("#",4X,A,"=",A)') 'trans_longi', trans_longi
      write(fh_variables_log, '("#",4X,A,"=",A)') 'ae_shape1', ae_shape1
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'e_impulse', e_impulse
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'amplitude1', amplitude1
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaser_int_wcm2_1', rlaser_int_wcm2_1
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'pulse_tw1', pulse_tw1
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'omega1', omega1
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_re1(1)', epdir_re1(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_re1(2)', epdir_re1(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_re1(3)', epdir_re1(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_im1(1)', epdir_im1(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_im1(2)', epdir_im1(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_im1(3)', epdir_im1(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'phi_cep1', phi_cep1
      write(fh_variables_log, '("#",4X,A,"=",A)') 'ae_shape2', ae_shape2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'amplitude2', amplitude2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaser_int_wcm2_2', rlaser_int_wcm2_2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'pulse_tw2', pulse_tw2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'omega2', omega2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_re2(1)', epdir_re2(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_re2(2)', epdir_re2(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_re2(3)', epdir_re2(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_im2(1)', epdir_im2(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_im2(2)', epdir_im2(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'epdir_im2(3)', epdir_im2(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'phi_cep2', phi_cep2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 't1_t2', t1_t2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 't1_delay', t1_delay
      write(fh_variables_log, '("#",4X,A,"=",A)') 'quadrupole', quadrupole
      write(fh_variables_log, '("#",4X,A,"=",A)') 'quadrupole_pot', quadrupole_pot
      write(fh_variables_log, '("#",4X,A,"=",A)') 'alocal_laser', alocal_laser
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaserbound_sta(1)', rlaserbound_sta(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaserbound_sta(2)', rlaserbound_sta(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaserbound_sta(3)', rlaserbound_sta(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaserbound_end(1)', rlaserbound_end(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaserbound_end(2)', rlaserbound_end(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaserbound_end(3)', rlaserbound_end(3)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nump', nump
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'vecp(1,1)', vecp(1,1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'vecp(2,1)', vecp(2,1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'vecp(3,1)', vecp(3,1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'vecp(1,2)', vecp(1,2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'vecp(2,2)', vecp(2,2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'vecp(3,2)', vecp(3,2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'coop(1,1)', coop(1,1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'coop(2,1)', coop(2,1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'coop(3,1)', coop(3,1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'coop(1,2)', coop(1,2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'coop(2,2)', coop(2,2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'coop(3,2)', coop(3,2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'radp_diele', radp_diele

      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'multiscale', inml_multiscale
      write(fh_variables_log, '("#",4X,A,"=",A)') 'fdtddim', fdtddim
      write(fh_variables_log, '("#",4X,A,"=",A)') 'twod_shape', twod_shape
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nx_m', nx_m
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'ny_m', ny_m
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nz_m', nz_m
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'hx_m', hx_m
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'hy_m', hy_m
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'hz_m', hz_m
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nksplit', nksplit
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nxysplit', nxysplit
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nxvacl_m', nxvacl_m
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nxvacr_m', nxvacr_m
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nx_origin_m', nx_origin_m
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'ny_origin_m', ny_origin_m
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nz_origin_m', nz_origin_m
      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_macropoint', trim(file_macropoint)
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'num_macropoint', num_macropoint
      write(fh_variables_log, '("#",4X,A,"=",A)') 'set_ini_coor_vel', set_ini_coor_vel
      write(fh_variables_log, '("#",4X,A,"=",I5)') 'nmacro_write_group', nmacro_write_group

      if(inml_maxwell >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'maxwell', inml_maxwell
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'al_em(1)', al_em(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'al_em(2)', al_em(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'al_em(3)', al_em(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl_em(1)', dl_em(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl_em(2)', dl_em(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl_em(3)', dl_em(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dt_em', dt_em
      write(fh_variables_log, '("#",4X,A,"=",I6)')     'nt_em', nt_em
      write(fh_variables_log, '("#",4X,A,"=",I6)')     'iboundary(1,1)', iboundary(1,1)
      write(fh_variables_log, '("#",4X,A,"=",I6)')     'iboundary(1,2)', iboundary(1,2)
      write(fh_variables_log, '("#",4X,A,"=",I6)')     'iboundary(2,1)', iboundary(2,1)
      write(fh_variables_log, '("#",4X,A,"=",I6)')     'iboundary(2,2)', iboundary(2,2)
      write(fh_variables_log, '("#",4X,A,"=",I6)')     'iboundary(3,1)', iboundary(3,1)
      write(fh_variables_log, '("#",4X,A,"=",I6)')     'iboundary(3,2)', iboundary(3,2)
      write(fh_variables_log, '("#",4X,A,"=",A)')      'wave_input', wave_input
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'ek_dir1(1)', ek_dir1(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'ek_dir1(2)', ek_dir1(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'ek_dir1(3)', ek_dir1(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'source_loc1(1)', source_loc1(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'source_loc1(2)', source_loc1(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'source_loc1(3)', source_loc1(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'ek_dir2(1)', ek_dir2(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'ek_dir2(2)', ek_dir2(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'ek_dir2(3)', ek_dir2(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'source_loc2(1)', source_loc2(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'source_loc2(2)', source_loc2(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'source_loc2(3)', source_loc2(3)
      write(fh_variables_log, '("#",4X,A,"=",I6)')     'iobs_num_em', iobs_num_em
      write(fh_variables_log, '("#",4X,A,"=",I6)')     'iobs_samp_em', iobs_samp_em
      if(iobs_num_em==0) then
        write(fh_variables_log, '("#",4X,A,"=",3ES14.5)') 'obs_loc_em', obs_loc_em(1,:)
      else
        do i = 1,iobs_num_em
          write(fh_variables_log, '("#",4X,A,I3,A,"=",3ES14.5)') 'obs_loc_em(',i,',:)', obs_loc_em(i,:)
        end do
      end if
      write(fh_variables_log, '("#",4X,A,"=",A)')      'shape_file', trim(shape_file)
      write(fh_variables_log, '("#",4X,A,"=",I6)')     'imedia_num', imedia_num
      do i = 0,imedia_num
        write(fh_variables_log, '("#",4X,A,I3,A,"=",A)')      'type_media(',i,')', type_media(i)
        write(fh_variables_log, '("#",4X,A,I3,A,"=",ES12.5)') 'epsilon(',i,')', epsilon(i)
        write(fh_variables_log, '("#",4X,A,I3,A,"=",ES12.5)') 'rmu(',i,')', rmu(i)
        write(fh_variables_log, '("#",4X,A,I3,A,"=",ES12.5)') 'sigma(',i,')', sigma(i)
        write(fh_variables_log, '("#",4X,A,I3,A,"=",ES12.5)') 'omega_p_d(',i,')', omega_p_d(i)
        write(fh_variables_log, '("#",4X,A,I3,A,"=",ES12.5)') 'gamma_d(',i,')', gamma_d(i)
      end do
      write(fh_variables_log, '("#",4X,A,"=",A)')      'smooth_d', smooth_d
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'weight_d', weight_d
      write(fh_variables_log, '("#",4X,A,"=",A)')      'wf_em', wf_em

      if(inml_analysis >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'analysis', inml_analysis
      write(fh_variables_log, '("#",4X,A,"=",A)') 'projection_option', projection_option
      write(fh_variables_log, '("#",4X,A,"=",A)') 'projection_decomp', projection_decomp
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'nenergy', nenergy
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'de', de
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_psi', out_psi
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_dos', out_dos
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'out_dos_start', out_dos_start
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'out_dos_end', out_dos_end
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iout_dos_nenergy', iout_dos_nenergy
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'out_dos_smearing', out_dos_smearing
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_dos_method', out_dos_method
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_dos_fshift', out_dos_fshift
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_pdos', out_pdos
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_dns', out_dns
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_old_dns', out_old_dns
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_dns_rt', out_dns_rt
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_dns_rt_step', out_dns_rt_step
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_dns_trans', out_dns_trans
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'out_dns_trans_energy', out_dns_trans_energy
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_elf', out_elf
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_elf_rt', out_elf_rt
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_elf_rt_step', out_elf_rt_step
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_estatic_rt', out_estatic_rt
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_estatic_rt_step', out_estatic_rt_step
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_rvf_rt', out_rvf_rt
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_rvf_rt_step', out_rvf_rt_step
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_tm', out_tm
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_projection_step', out_projection_step
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_ms_step', out_ms_step
      write(fh_variables_log, '("#",4X,A,"=",A)') 'format3d', format3d
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'numfiles_out_3d', numfiles_out_3d
      write(fh_variables_log, '("#",4X,A,"=",A)') 'timer_process', timer_process

      if(inml_hartree >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'hartree', inml_hartree
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'meo', meo
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_pole_xyz(1)', num_pole_xyz(1)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_pole_xyz(2)', num_pole_xyz(2)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_pole_xyz(3)', num_pole_xyz(3)

      if(inml_ewald >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'ewald', inml_ewald
      write(fh_variables_log, '("#",4X,A,"=",I3)') 'newald', newald
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'aewald', aewald

      if(inml_opt >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'opt', inml_opt
     !write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cg_alpha_ini', cg_alpha_ini !not use now
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cg_alpha_up', cg_alpha_up
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cg_alpha_down', cg_alpha_down
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'convrg_scf_force', convrg_scf_force
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'convrg_scf_ene', convrg_scf_ene
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'convrg_opt_fmax',convrg_opt_fmax
     !write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'convrg_opt_ene', convrg_opt_ene !not use now
      if(inml_md >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'md', inml_md
      write(fh_variables_log, '("#",4X,A,"=",A)') 'ensemble', ensemble
      write(fh_variables_log, '("#",4X,A,"=",A)') 'thermostat', thermostat
      write(fh_variables_log, '("#",4X,A,"=",I8)') 'step_velocity_scaling', step_velocity_scaling
      write(fh_variables_log, '("#",4X,A,"=",I8)') 'step_update_ps', step_update_ps
      write(fh_variables_log, '("#",4X,A,"=",I8)') 'step_update_ps2', step_update_ps2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'temperature0_ion', temperature0_ion
      write(fh_variables_log, '("#",4X,A,"=",A)') 'set_ini_velocity', set_ini_velocity
      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_ini_velocity', trim(file_ini_velocity)
!      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_set_shake', trim(file_set_shake)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'thermostat_tau', thermostat_tau
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'friction', friction
      write(fh_variables_log, '("#",4X,A,"=",A)') 'stop_system_momt', stop_system_momt
      
      if(inml_misc >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'misc', inml_misc
      write(fh_variables_log, '("#",4X,A,"=",A)') 'fourier', fourier

      if(inml_group_fundamental >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'group_fundamental', inml_group_fundamental
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iditerybcg', iditerybcg
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iditer_nosubspace_diag', iditer_nosubspace_diag
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'ntmg', ntmg
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'idisnum(1)', idisnum(1)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'idisnum(2)', idisnum(2)
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iwrite_projection', iwrite_projection
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'itwproj', itwproj
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iwrite_projnum', iwrite_projnum
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'itcalc_ene', itcalc_ene

      if(inml_group_parallel >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'group_parallel', inml_group_parallel
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'isequential', isequential
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'imesh_s_all', imesh_s_all
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iflag_comm_rho', iflag_comm_rho

      if(inml_group_hartree >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'group_hartree', inml_group_hartree
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'hconv', hconv
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'lmax_meo', lmax_meo

      if(inml_group_file >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'group_file', inml_group_file
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'ic', ic
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'oc', oc
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'ic_rt', ic_rt
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'oc_rt', oc_rt

      if(inml_group_others >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'group_others', inml_group_others
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iparaway_ob', iparaway_ob
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iscf_order', iscf_order
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iswitch_orbital_mesh', iswitch_orbital_mesh
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iflag_psicube', iflag_psicube
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'lambda1_diis', lambda1_diis
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'lambda2_diis', lambda2_diis
      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_ini', file_ini
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iparaway_ob', iparaway_ob
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'num_projection', num_projection
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'num_projection', num_projection
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iwrite_projection_ob(1)', iwrite_projection_ob(1)
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iwrite_projection_ob(2)', iwrite_projection_ob(2)
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iwrite_projection_k(1)', iwrite_projection_k(1)
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iwrite_projection_k(2)', iwrite_projection_k(2)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'filename_pot', filename_pot
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iwrite_external', iwrite_external
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iflag_dip2', iflag_dip2
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iflag_intelectron', iflag_intelectron
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'num_dip2', num_dip2
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dip2boundary(1)', dip2boundary(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dip2boundary(2)', dip2boundary(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dip2center(1)', dip2center(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dip2center(2)', dip2center(2)
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iflag_fourier_omega', iflag_fourier_omega
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'num_fourier_omega', num_fourier_omega
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'fourier_omega(1)', fourier_omega(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'fourier_omega(2)', fourier_omega(2)
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'itotntime2', itotntime2
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iwdenoption', iwdenoption
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'iwdenstep', iwdenstep
      write(fh_variables_log, '("#",4X,A,"=",I2)') 'iflag_estatic', iflag_estatic


      select case(iflag_atom_coor)
      case(ntype_atom_coor_cartesian)
        write(fh_variables_log, '("#namelist: ",A)') 'atom_coor'
        do i = 1,natom
          write(fh_variables_log, '("#",4X,A,I2,A,"=",3ES14.5)') 'rion(',i,')', rion(1:3,i)
        end do
      case(ntype_atom_coor_reduced)
        write(fh_variables_log, '("#namelist: ",A)') 'atom_red_coor'
        do i = 1,natom
          write(fh_variables_log, '("#",4X,A,I2,A,"=",3ES14.5)') 'rion_red(',i,')', rion_red(1:3,i)
        end do
      case default
      end select

      close(fh_variables_log)

    end if

    call comm_bcast(ierr_nml,nproc_group_global)
    if(ierr_nml > 0)then
      if (comm_is_root(nproc_id_global)) write(*,"(I4,2x,A)")ierr_nml,'error(s) in input.'
      call end_parallel
      stop
    end if

    !(create output directory)
    if (comm_is_root(nproc_id_global)) then
       if(directory(1:3).ne."./ ") call create_directory(directory)
    endif

  end subroutine dump_input_common

  subroutine check_bad_input
    use salmon_parallel
    use salmon_communication
    implicit none

    !! Add wrong input keyword or wrong/unavailable input combinations here
    !! (now only a few)

    if(iperiodic==0.or.(iperiodic==3.and.domain_parallel=='y')) then
      select case(convergence)
      case('rho_dne')
        continue
      case('norm_rho','norm_rho_dng')
        if(threshold_norm_rho<-1.d-12)then
          if (comm_is_root(nproc_id_global)) then
            write(*,*) 'set threshold_norm_rho when convergence is norm_rho or norm_rho_dng.'
          endif
          call end_parallel
        end if
      case('norm_pot','norm_pot_dng')
        if(threshold_norm_pot<-1.d-12)then
          if (comm_is_root(nproc_id_global)) then
            write(*,*) 'set threshold_norm_pot when convergence is norm_pot or norm_rho_pot.'
          endif
          call end_parallel
        end if
      case default
        if (comm_is_root(nproc_id_global)) then
          write(*,*) 'check a keyword of convergence.'
        endif
        call end_parallel
      end select

    else if(iperiodic==3.and.domain_parallel=='n') then
      if(convergence.ne.'rho_dne') call stop_by_bad_input2('iperiodic','convergence')
      if(abs(t1_delay).ge.1d-10)then
         if(index(ae_shape1,'Acos')==0) call stop_by_bad_input2('t1_delay','ae_shape1')
      endif
    endif
      
    select case(fourier)
    case('ft','FT','ffte','FFTE')
      continue
    case default
      if (comm_is_root(nproc_id_global)) then
        write(*,*) "keyword 'fourier' must be 'FT' or 'FFTE'"
      end if
      call end_parallel
    end select

  end subroutine check_bad_input

  subroutine stop_by_bad_input2(inp1,inp2,inp3)
    use salmon_parallel
    use salmon_communication
    implicit none
    character(*) :: inp1
    character(*) :: inp2
    character(*),optional :: inp3
    if (comm_is_root(nproc_id_global)) then
      write(*,*) ' Bad input combination: '
      if(present(inp3))then
        write(*,*) ' check keywords of ',trim(inp1),' and ',trim(inp2),' and ',trim(inp3)
      else
        write(*,*) ' check keywords of ',trim(inp1),' and ',trim(inp2)
      end if
    endif
    call end_parallel
    stop
  end subroutine stop_by_bad_input2


end module inputoutput
