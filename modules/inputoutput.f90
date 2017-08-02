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
  integer :: inml_analysis
  integer :: inml_hartree
  integer :: inml_ewald
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

contains
  subroutine read_input
    implicit none

    call read_stdin
    call read_input_common ! Should be renamed properly later
    if(restart_option == 'restart')return
    call read_atomic_coordinates
    call dump_input_common ! Should be renamed properly later

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

    namelist/calculation/ &
      & calc_mode, &
      & use_ehrenfest_md, &
      & use_ms_maxwell, &
      & use_force, &
      & use_geometry_opt

    namelist/control/ &
      & restart_option, &
      & backup_frequency, &
      & time_shutdown, &
      & sysname, &
      & directory, &
      & dump_filename

    namelist/units/ &
      & unit_time, &
      & unit_length, &
      & unit_energy, &
      & unit_charge

    namelist/parallel/ &
      & domain_parallel, &
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
      & nelem, &
      & natom, &
      & file_atom_coor, &
      & file_atom_red_coor

    namelist/pseudo/ &
      & pseudo_file, &
      & Lmax_ps, &
      & Lloc_ps, &
      & iZatom, &
      & psmask_option, &
      & alpha_mask, &
      & gamma_mask, &
      & eta_mask

    namelist/functional/ &
      & xc, &
      & cval

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
      & threshold_pot

    namelist/emfield/ &
      & trans_longi, &
      & ae_shape1, &
      & e_impulse, &
      & amplitude1, &
      & rlaser_int1, &
      & pulse_tw1, &
      & omega1, &
      & epdir_re1, &
      & epdir_im1, &
      & phi_cep1, &
      & ae_shape2, &
      & amplitude2, &
      & rlaser_int2, &
      & pulse_tw2, &
      & omega2, &
      & epdir_re2, &
      & epdir_im2, &
      & phi_cep2, &
      & t1_t2, &
      & quadrupole, &
      & quadrupole_pot, &
      & alocal_laser , &
      & rlaserbound_sta , &
      & rlaserbound_end

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
      & nxvacr_m

    namelist/analysis/ &
      & projection_option, &
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
      & out_elf, &
      & out_dns_rt, &
      & out_dns_rt_step, &
      & out_elf_rt, &
      & out_elf_rt_step, &
      & out_estatic_rt, &
      & out_estatic_rt_step, &
      & format3d, &
      & numfiles_out_3d, &
      & timer_process

    namelist/hartree/ &
      & meo, &
      & num_pole_xyz

    namelist/ewald/ &
      & newald, &
      & aewald

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
    unit_time='au'
    unit_length='au'
    unit_energy='au'
    unit_charge='au'
!! =======================

    if (comm_is_root(nproc_id_global)) then
      open(fh_namelist, file='.namelist.tmp', status='old')
      read(fh_namelist, nml=units, iostat=inml_units)
      rewind(fh_namelist)
      close(fh_namelist)
    end if

    call comm_bcast(unit_time,  nproc_group_global)
    call comm_bcast(unit_length,nproc_group_global)
    call comm_bcast(unit_energy,nproc_group_global)
    call comm_bcast(unit_charge,nproc_group_global)

    call initialize_inputoutput_units


!! == default for &calculation 
    calc_mode        = 'none'
    use_ehrenfest_md = 'n'
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
!! == default for &parallel
    domain_parallel   = 'n'
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
    nelem              = 0
    natom              = 0
    file_atom_coor          = 'none'
    file_atom_red_coor          = 'none'
!! == default for &pseudo
    pseudo_file     = 'none'
    Lmax_ps       = -1
    Lloc_ps       = -1
    iZatom        = -1
    psmask_option = 'n'
    alpha_mask    = 0.8d0
    gamma_mask    = 1.8d0
    eta_mask      = 15d0
!! == default for &functional
    xc   = 'PZ'
    cval = -1d0
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
    convergence   = 'rho_dng'
    threshold     = 1d-17/ulength_from_au**3  ! a.u., 1d-17 a.u. = 6.75d-17 AA**(-3)
    threshold_pot = -1d0*uenergy_from_au**2*uenergy_from_au**3  ! a.u., -1 a.u. = -1.10d2 eV**2*AA**3
!! == default for &emfield
    trans_longi    = 'tr'
    ae_shape1      = 'none'
    e_impulse = 1d-2*uenergy_from_au/ulength_from_au*utime_from_au ! a.u.
    amplitude1     = 0d0
    rlaser_int1    = -1d0
    pulse_tw1      = 0d0
    omega1         = 0d0
    epdir_re1      = (/1d0,0d0,0d0/)
    epdir_im1      = 0d0
    phi_cep1       = 0d0
    ae_shape2      = 'none'
    amplitude2     = 0d0
    rlaser_int2    = -1d0
    pulse_tw2      = 0d0
    omega2         = 0d0
    epdir_re2      = (/1d0,0d0,0d0/)
    epdir_im2      = 0d0
    phi_cep2       = 0d0
    t1_t2          = 0d0
    quadrupole     = 'n'
    quadrupole_pot = ''
    alocal_laser    = 'n'
    rlaserbound_sta(1) = -1.d7*ulength_from_au ! a.u.
    rlaserbound_sta(2) = -1.d7*ulength_from_au ! a.u.
    rlaserbound_sta(3) = -1.d7*ulength_from_au ! a.u.
    rlaserbound_end(1) = 1.d7*ulength_from_au ! a.u.
    rlaserbound_end(2) = 1.d7*ulength_from_au ! a.u.
    rlaserbound_end(3) = 1.d7*ulength_from_au ! a.u.
!! == default for &multiscale
    fdtddim    = '1d'
    twod_shape = 'periodic'
    nx_m       = 1
    ny_m       = 1
    nz_m       = 1
    hx_m       = 1d0
    hy_m       = 1d0
    hz_m       = 1d0
    nksplit    = 1
    nxysplit   = 1
    nxvacl_m   = 0
    nxvacr_m   = 0
!! == default for &analysis
    projection_option   = 'no'
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
    out_elf             = 'n'
    out_dns_rt          = 'n'
    out_dns_rt_step     = 50
    out_elf_rt          = 'n'
    out_elf_rt_step     = 50
    out_estatic_rt      = 'n'
    out_estatic_rt_step = 50
    format3d            = 'cube'
    numfiles_out_3d     = 1
    timer_process       = 'n'
!! == default for &hartree
    meo          = 3
    num_pole_xyz = 0
!! == default for &ewald
    newald = 4
    aewald = 0.5d0
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

      read(fh_namelist, nml=analysis, iostat=inml_analysis)
      rewind(fh_namelist)

      read(fh_namelist, nml=hartree, iostat=inml_hartree)
      rewind(fh_namelist)

      read(fh_namelist, nml=ewald, iostat=inml_ewald)
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
    call comm_bcast(calc_mode       ,nproc_group_global)
    call comm_bcast(use_ehrenfest_md,nproc_group_global)
    call comm_bcast(use_ms_maxwell  ,nproc_group_global)
    call comm_bcast(use_force       ,nproc_group_global)
    call comm_bcast(use_geometry_opt,nproc_group_global)
!! == bcast for &control
    call comm_bcast(restart_option  ,nproc_group_global)
    call comm_bcast(backup_frequency,nproc_group_global)
    call comm_bcast(time_shutdown   ,nproc_group_global)
    call comm_bcast(sysname         ,nproc_group_global)
    call comm_bcast(directory       ,nproc_group_global)
    call comm_bcast(dump_filename   ,nproc_group_global)
!! == bcast for &parallel
    call comm_bcast(domain_parallel  ,nproc_group_global)
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
    call comm_bcast(nelem              ,nproc_group_global)
    call comm_bcast(natom              ,nproc_group_global)
    call comm_bcast(file_atom_coor     ,nproc_group_global)
    call comm_bcast(file_atom_red_coor ,nproc_group_global)
!! == bcast for &pseudo
    call comm_bcast(pseudo_file  ,nproc_group_global)
    call comm_bcast(Lmax_ps      ,nproc_group_global)
    call comm_bcast(Lloc_ps      ,nproc_group_global)
    call comm_bcast(iZatom       ,nproc_group_global)
    call comm_bcast(psmask_option,nproc_group_global)
    call comm_bcast(alpha_mask   ,nproc_group_global)
    call comm_bcast(gamma_mask   ,nproc_group_global)
    call comm_bcast(eta_mask     ,nproc_group_global)
!! == bcast for &functional
    call comm_bcast(xc  ,nproc_group_global)
    call comm_bcast(cval,nproc_group_global)
!! == bcast for &rgrid
    call comm_bcast(dl,nproc_group_global)
    dl = dl * ulength_to_au
    call comm_bcast(num_rgrid,nproc_group_global)
!! == bcast for &kgrid
    call comm_bcast(num_kgrid,nproc_group_global)
    call comm_bcast(file_kw  ,nproc_group_global)
!! == bcast for &kgrid
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
    threshold = threshold / (ulength_to_au)**3
    call comm_bcast(threshold_pot           ,nproc_group_global)
    threshold_pot = threshold_pot * (uenergy_to_au)**2 * (ulength_to_au)**3 
!! == bcast for &emfield
    call comm_bcast(trans_longi,nproc_group_global)
    call comm_bcast(ae_shape1  ,nproc_group_global)
    call comm_bcast(e_impulse,nproc_group_global)
    e_impulse = e_impulse *uenergy_to_au/ulength_to_au*utime_to_au
    call comm_bcast(amplitude1 ,nproc_group_global)
    amplitude1 = amplitude1*(uenergy_to_au/ulength_to_au/ucharge_to_au)
    call comm_bcast(rlaser_int1,nproc_group_global)
    call comm_bcast(pulse_tw1  ,nproc_group_global)
    pulse_tw1 = pulse_tw1 * utime_to_au
    call comm_bcast(omega1,nproc_group_global)
    omega1 = omega1 * uenergy_to_au
    call comm_bcast(epdir_re1 ,nproc_group_global)
    call comm_bcast(epdir_im1 ,nproc_group_global)
    call comm_bcast(phi_cep1  ,nproc_group_global)
    call comm_bcast(ae_shape2 ,nproc_group_global)
    call comm_bcast(amplitude2,nproc_group_global)
    amplitude2 = amplitude2*(uenergy_to_au/ulength_to_au/ucharge_to_au)
    call comm_bcast(rlaser_int2,nproc_group_global)
    call comm_bcast(pulse_tw2  ,nproc_group_global)
    pulse_tw2 = pulse_tw2 * utime_to_au
    call comm_bcast(omega2,nproc_group_global)
    omega2 = omega2 * uenergy_to_au
    call comm_bcast(epdir_re2,nproc_group_global)
    call comm_bcast(epdir_im2,nproc_group_global)
    call comm_bcast(phi_cep2 ,nproc_group_global)
    call comm_bcast(t1_t2    ,nproc_group_global)
    t1_t2 = t1_t2 * utime_to_au
    call comm_bcast(quadrupole    ,nproc_group_global)
    call comm_bcast(quadrupole_pot,nproc_group_global)
    call comm_bcast(alocal_laser  ,nproc_group_global)
    call comm_bcast(rlaserbound_sta,nproc_group_global)
    rlaserbound_sta = rlaserbound_sta * ulength_to_au
    call comm_bcast(rlaserbound_end,nproc_group_global)
    rlaserbound_end = rlaserbound_end * ulength_to_au
!! == bcast for &multiscale
    call comm_bcast(fdtddim   ,nproc_group_global)
    call comm_bcast(twod_shape,nproc_group_global)
    call comm_bcast(nx_m      ,nproc_group_global)
    call comm_bcast(ny_m      ,nproc_group_global)
    call comm_bcast(nz_m      ,nproc_group_global)
    call comm_bcast(hx_m      ,nproc_group_global)
    hx_m = hx_m * ulength_to_au
    call comm_bcast(hy_m,nproc_group_global)
    hy_m = hy_m * ulength_to_au
    call comm_bcast(hz_m,nproc_group_global)
    hz_m = hz_m * ulength_to_au
    call comm_bcast(nksplit ,nproc_group_global)
    call comm_bcast(nxysplit,nproc_group_global)
    call comm_bcast(nxvacl_m,nproc_group_global)
    call comm_bcast(nxvacr_m,nproc_group_global)
!! == bcast for &analysis
    call comm_bcast(projection_option,nproc_group_global)
    call comm_bcast(nenergy          ,nproc_group_global)
    call comm_bcast(de               ,nproc_group_global)
    de = de * uenergy_to_au
    call comm_bcast(out_psi            ,nproc_group_global)
    call comm_bcast(out_dos            ,nproc_group_global)
    call comm_bcast(out_dos_start      ,nproc_group_global)
    out_dos_start = out_dos_start * uenergy_to_au
    call comm_bcast(out_dos_end        ,nproc_group_global)
    out_dos_end = out_dos_end * uenergy_to_au
    call comm_bcast(iout_dos_nenergy   ,nproc_group_global)
    call comm_bcast(out_dos_smearing   ,nproc_group_global)
    out_dos_smearing = out_dos_smearing * uenergy_to_au
    call comm_bcast(out_dos_method     ,nproc_group_global)
    call comm_bcast(out_dos_fshift     ,nproc_group_global)
    call comm_bcast(out_pdos           ,nproc_group_global)
    call comm_bcast(out_dns            ,nproc_group_global)
    call comm_bcast(out_elf            ,nproc_group_global)
    call comm_bcast(out_dns_rt         ,nproc_group_global)
    call comm_bcast(out_dns_rt_step    ,nproc_group_global)
    call comm_bcast(out_elf_rt         ,nproc_group_global)
    call comm_bcast(out_elf_rt_step    ,nproc_group_global)
    call comm_bcast(out_estatic_rt     ,nproc_group_global)
    call comm_bcast(out_estatic_rt_step,nproc_group_global)
    call comm_bcast(format3d           ,nproc_group_global)
    call comm_bcast(numfiles_out_3d    ,nproc_group_global)
    call comm_bcast(timer_process      ,nproc_group_global)
!! == bcast for &hartree
    call comm_bcast(meo         ,nproc_group_global)
    call comm_bcast(num_pole_xyz,nproc_group_global)
!! == bcast for &ewald
    call comm_bcast(newald,nproc_group_global)
    call comm_bcast(aewald,nproc_group_global)
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

    if(icount/=1)then
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

    allocate(Rion(3,natom), Rion_red(3,natom),Kion(natom), flag_geo_opt_atom(natom))
    Rion = 0d0
    Rion_red = 0d0
    Kion = 0
    flag_geo_opt_atom = 'n'

    if (comm_is_root(nproc_id_global))then
      open(fh_atomic_coor, file=filename_tmp, status='old')
      select case(iflag_atom_coor)
      case(ntype_atom_coor_cartesian)
         do i=1, natom
            if(use_geometry_opt == 'y')then
               read(fh_atomic_coor, *) char_atom, Rion(:,i), Kion(i), flag_geo_opt_atom(i)
            else
               read(fh_atomic_coor, *) char_atom, Rion(:,i), Kion(i)
            end if
         end do
         Rion = Rion*ulength_to_au
      case(ntype_atom_coor_reduced)
         do i=1, natom
            if(use_geometry_opt == 'y')then
               read(fh_atomic_coor, *) char_atom, Rion_red(:,i), Kion(i), flag_geo_opt_atom(i)
            else
               read(fh_atomic_coor, *) char_atom, Rion_red(:,i), Kion(i)
            end if
         end do
      end select
      close(fh_atomic_coor)
    end if

    call comm_bcast(Rion,nproc_group_global)
    call comm_bcast(Rion_red,nproc_group_global)
    call comm_bcast(Kion,nproc_group_global)
    call comm_bcast(flag_geo_opt_atom,nproc_group_global)


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




  end subroutine initialize_inputoutput_units

  subroutine dump_input_common
    use salmon_parallel
    use salmon_communication
    use salmon_file, only: get_filehandle
    implicit none
    integer :: i,ierr_nml
    ierr_nml = 0

    if (comm_is_root(nproc_id_global)) then

      fh_variables_log = get_filehandle()
      open(fh_variables_log,file='variables.log')

      if(inml_calculation >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'calculation', inml_calculation
      write(fh_variables_log, '("#",4X,A,"=",A)') 'calc_mode', calc_mode
      write(fh_variables_log, '("#",4X,A,"=",A)') 'use_ehrenfest_md', use_ehrenfest_md
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

      if(inml_units >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'units', inml_units
      write(fh_variables_log, '("#",4X,A,"=",A)') 'unit_time', unit_time
      write(fh_variables_log, '("#",4X,A,"=",A)') 'unit_length', unit_length
      write(fh_variables_log, '("#",4X,A,"=",A)') 'unit_energy', unit_energy
      write(fh_variables_log, '("#",4X,A,"=",A)') 'unit_charge', unit_charge

      if(inml_parallel >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'parallel', inml_parallel
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
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'nelem', nelem
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'natom', natom
      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_atom_coor', trim(file_atom_coor)
      write(fh_variables_log, '("#",4X,A,"=",A)') 'file_atom_red_coor', trim(file_atom_red_coor)

      if(inml_pseudo >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'pseudo', inml_pseudo

      do i = 1,nelem
        write(fh_variables_log, '("#",4X,A,I2,A,"=",A)') 'pseudo_file(',i,')', trim(pseudo_file(i))
        write(fh_variables_log, '("#",4X,A,I2,A,"=",I4)') 'Lmax_ps(',i,')', Lmax_ps(i)
        write(fh_variables_log, '("#",4X,A,I2,A,"=",I4)') 'Lloc_ps(',i,')', Lloc_ps(i)
        write(fh_variables_log, '("#",4X,A,I2,A,"=",I4)') 'iZatom(',i,')', iZatom(i)
      end do
      write(fh_variables_log, '("#",4X,A,"=",A)') 'psmask_option', psmask_option
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'alpha_mask', alpha_mask
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'gamma_mask', gamma_mask
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'eta_mask', eta_mask

      if(inml_functional >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'functional', inml_functional
      write(fh_variables_log, '("#",4X,A,"=",A)') 'xc', trim(xc)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'cval', cval

      if(inml_rgrid >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'rgrid', inml_rgrid
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl(1)', dl(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl(2)', dl(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'dl(3)', dl(3)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_rgrid(1)', num_rgrid(1)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_rgrid(2)', num_rgrid(2)
      write(fh_variables_log, '("#",4X,A,"=",I4)') 'num_rgrid(3)', num_rgrid(3)

      if(inml_kgrid >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I1)') 'kgrid', inml_kgrid
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
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'threshold_pot', threshold_pot

      if(inml_emfield >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'emfield', inml_emfield
      write(fh_variables_log, '("#",4X,A,"=",A)') 'trans_longi', trans_longi
      write(fh_variables_log, '("#",4X,A,"=",A)') 'ae_shape1', ae_shape1
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'e_impulse', e_impulse
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'amplitude1', amplitude1
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaser_int1', rlaser_int1
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
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaser_int2', rlaser_int2
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
      write(fh_variables_log, '("#",4X,A,"=",A)') 'quadrupole', quadrupole
      write(fh_variables_log, '("#",4X,A,"=",A)') 'quadrupole_pot', quadrupole_pot
      write(fh_variables_log, '("#",4X,A,"=",A)') 'alocal_laser', alocal_laser
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaserbound_sta(1)', rlaserbound_sta(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaserbound_sta(2)', rlaserbound_sta(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaserbound_sta(3)', rlaserbound_sta(3)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaserbound_end(1)', rlaserbound_end(1)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaserbound_end(2)', rlaserbound_end(2)
      write(fh_variables_log, '("#",4X,A,"=",ES12.5)') 'rlaserbound_end(3)', rlaserbound_end(3)

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

      if(inml_analysis >0)ierr_nml = ierr_nml +1
      write(fh_variables_log, '("#namelist: ",A,", status=",I3)') 'analysis', inml_analysis
      write(fh_variables_log, '("#",4X,A,"=",A)') 'projection_option', projection_option
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
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_elf', out_elf
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_dns_rt', out_dns_rt
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_dns_rt_step', out_dns_rt_step
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_elf_rt', out_elf_rt
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_elf_rt_step', out_elf_rt_step
      write(fh_variables_log, '("#",4X,A,"=",A)') 'out_estatic_rt', out_estatic_rt
      write(fh_variables_log, '("#",4X,A,"=",I6)') 'out_estatic_rt_step', out_estatic_rt_step
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
          write(fh_variables_log, '("#",4X,A,I2,A,"=",3ES14.5)') 'Rion(',i,')', Rion(1:3,i)
        end do
      case(ntype_atom_coor_reduced)
        write(fh_variables_log, '("#namelist: ",A)') 'atom_red_coor'
        do i = 1,natom
          write(fh_variables_log, '("#",4X,A,I2,A,"=",3ES14.5)') 'Rion_red(',i,')', Rion_red(1:3,i)
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


  end subroutine dump_input_common
    
end module inputoutput
