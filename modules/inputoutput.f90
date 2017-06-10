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


  integer, parameter :: fh_namelist = 901
  integer, parameter :: fh_atomic_spiecies = 902
  integer, parameter :: fh_atomic_coor = 903
  integer, parameter :: fh_reentrance = 904
  integer, parameter :: fh_atomic_red_coor = 905
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
  integer :: inml_linear_response
  integer :: inml_multiscale
  integer :: inml_analysis
  integer :: inml_hartree
  integer :: inml_ewald


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
  real(8) :: utime_to_au, utime_from_au
  real(8) :: ulength_to_au, ulength_from_au
  real(8) :: uenergy_to_au, uenergy_from_au
  real(8) :: ucharge_to_au, ucharge_from_au

contains


  subroutine read_stdin
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    implicit none

    integer :: cur = fh_namelist
    integer :: ret = 0
    character(100) :: buff, text


    
    if (comm_is_root(nproc_id_global)) then
      open(fh_namelist, file='.namelist.tmp', status='replace')
!      open(fh_atomic_spiecies, file='.atomic_spiecies.tmp', status='replace')
      if_nml_coor =.false. 
      open(fh_atomic_coor, file='.atomic_coor.tmp', status='replace')
      if_nml_red_coor = .false.
      open(fh_atomic_red_coor, file='.atomic_red_coor.tmp', status='replace')
      open(fh_reentrance, file='.reenetrance.tmp', status='replace')
      
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
    implicit none

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
      & directory

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
      & nelec, &
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
      & ncg, &
      & nmemory_mb, &
      & alpha_mb, &
      & fsset_option, &
      & nfsset_start, &
      & nfsset_every, &
      & nscf, &
      & ngeometry_opt, &
      & subspace_diagonalization, &
      & cmixing, &
      & rmixrate, &
      & convergence, &
      & threshold

    namelist/emfield/ &
      & trans_longi, &
      & ae_shape1, &
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

    namelist/linear_response/ &
      & e_impulse

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
      & out_elf, &
      & out_dns_rt, &
      & out_dns_rt_step, &
      & out_elf_rt, &
      & out_elf_rt_step, &
      & out_estatic_rt, &
      & out_estatic_rt_step, &
      & format3d

    namelist/hartree/ &
      & meo, &
      & num_pole_xyz

    namelist/ewald/ &
      & newald, &
      & aewald


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
    nelec              = 0
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
    cval = 1d0
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
    ncg           = 5
    nmemory_mb    = 8
    alpha_mb      = 0.75d0
    fsset_option  = 'n'
    nfsset_start  = 75
    nfsset_every  = 25
    nscf          = 0
    ngeometry_opt = 1
    subspace_diagonalization = 'y'
    cmixing       = 'broyden'
    rmixrate      = 0.5d0
    convergence   = 'rho'
    threshold     = 1d-6
!! == default for &emfield
    trans_longi    = 'tr'
    ae_shape1      = 'none'
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
    rlaserbound_sta(1) = -1.d7/au_length_aa*ulength_from_au
    rlaserbound_sta(2) = -1.d7/au_length_aa*ulength_from_au
    rlaserbound_sta(3) = -1.d7/au_length_aa*ulength_from_au
    rlaserbound_end(1) = 1.d7/au_length_aa*ulength_from_au
    rlaserbound_end(2) = 1.d7/au_length_aa*ulength_from_au
    rlaserbound_end(3) = 1.d7/au_length_aa*ulength_from_au
!! == default for &linear_response
    e_impulse = 1d-2*uenergy_from_au/ulength_from_au*utime_from_au ! a.u.
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
!! == default for &hartree
    meo          = 3
    num_pole_xyz = 0
!! == default for &ewald
    newald = 4
    aewald = 0.5d0

    if (comm_is_root(nproc_id_global)) then
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

      read(fh_namelist, nml=linear_response, iostat=inml_linear_response)
      rewind(fh_namelist)

      read(fh_namelist, nml=multiscale, iostat=inml_multiscale)
      rewind(fh_namelist)

      read(fh_namelist, nml=analysis, iostat=inml_analysis)
      rewind(fh_namelist)

      read(fh_namelist, nml=hartree, iostat=inml_hartree)
      rewind(fh_namelist)

      read(fh_namelist, nml=ewald, iostat=inml_ewald)
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
    call comm_bcast(nelec              ,nproc_group_global)
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
    call comm_bcast(ncg                     ,nproc_group_global)
    call comm_bcast(nmemory_mb              ,nproc_group_global)
    call comm_bcast(alpha_mb                ,nproc_group_global)
    call comm_bcast(fsset_option            ,nproc_group_global)
    call comm_bcast(nfsset_start            ,nproc_group_global)
    call comm_bcast(nfsset_every            ,nproc_group_global)
    call comm_bcast(nscf                    ,nproc_group_global)
    call comm_bcast(ngeometry_opt           ,nproc_group_global)
    call comm_bcast(subspace_diagonalization,nproc_group_global)
    call comm_bcast(cmixing                 ,nproc_group_global)
    call comm_bcast(rmixrate                ,nproc_group_global)
    call comm_bcast(convergence             ,nproc_group_global)
    call comm_bcast(threshold               ,nproc_group_global)
!! == bcast for &emfield
    call comm_bcast(trans_longi,nproc_group_global)
    call comm_bcast(ae_shape1  ,nproc_group_global)
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
    call comm_bcast(rlaserbound_end,nproc_group_global)
!! == bcast for &linear_response
    call comm_bcast(e_impulse,nproc_group_global)
    e_impulse = e_impulse *uenergy_to_au/ulength_to_au*utime_to_au
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
!! == bcast for &hartree
    call comm_bcast(meo         ,nproc_group_global)
    call comm_bcast(num_pole_xyz,nproc_group_global)
!! == bcast for &ewald
    call comm_bcast(newald,nproc_group_global)
    call comm_bcast(aewald,nproc_group_global)

  end subroutine read_input_common

  subroutine read_atomic_coordinates
    use salmon_parallel
    use salmon_communication
    character(256) :: filename_tmp,char_atom
    integer :: icount,i
    logical :: if_error, if_cartesian


    if (comm_is_root(nproc_id_global)) then


      if_error = .false.
      icount = 0
      if(file_atom_coor /= 'none')then
        icount = icount + 1
        if_cartesian = .true.
        filename_tmp = trim(file_atom_coor)
      end if

      if(file_atom_red_coor /= 'none')then
        icount = icount + 1
        if_cartesian = .false.
        filename_tmp = trim(file_atom_coor)
      end if

      if(if_nml_red_coor)then
        icount = icount + 1
        if_cartesian = .true.
        filename_tmp = '.atomic_coor.tmp'
      end if

      if(if_nml_red_coor)then
        icount = icount + 1
        if_cartesian = .false.
        filename_tmp = '.atomic_red_coor.tmp'
      end if

      if(icount /= 1)if_error = .true.
    end if

    call comm_bcast(if_error,nproc_group_global)
    if(if_error)then
       if (comm_is_root(nproc_id_global))then
         write(*,"(I4,2x,A)")'Error in input: The following inputs are incompatible.'
         write(*,"(I4,2x,A)")'file_atom_coor, file_atom_red_coor, &atomic_coor, and &atomic_red_coor.'
       end if
       call end_parallel
       stop
    end if

    if( (.not.if_cartesian) .and. iperiodic == 0)then
       if (comm_is_root(nproc_id_global))then
         write(*,"(I4,2x,A)")'Error in input: Reduced coordinate is invalid for isolated systems.'
       end if
       call end_parallel
       stop
    end if

    allocate(Rion(3,natom), Kion(natom), flag_geo_opt_atom(natom))
    Rion = 0d0
    Kion = 0
    flag_geo_opt_atom = 'n'

    if (comm_is_root(nproc_id_global))then
      open(fh_atomic_coor, file=filename_tmp, status='old')
      do i=1, natom
        if(use_geometry_opt == 'y')then
          read(fh_atomic_coor, *) char_atom, Rion(:,i), Kion(i), flag_geo_opt_atom(i)
        else
          read(fh_atomic_coor, *) char_atom, Rion(:,i), Kion(i)
        end if
      end do
      close(fh_atomic_coor)

    if(if_cartesian)Rion = Rion*ulength_to_au
    if(if_cartesian .and. iperiodic == 3)then
      Rion(1,:) = Rion(1,:)/al(1)
      Rion(2,:) = Rion(2,:)/al(2)
      Rion(3,:) = Rion(3,:)/al(3)
    end if
    end if

    call comm_bcast(Rion,nproc_group_global)
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
    case('fs','femtosecond')
      utime_to_au   = 1d0/au_time_fs
      utime_from_au = au_time_fs
    case default
      stop "Invalid unit for time."
    end select

! Unit for length
    select case(unit_length)
    case('au','a.u.')
      ulength_to_au   = 1d0
      ulength_from_au = 1d0
    case('AA','angstrom','Angstrom')
      ulength_to_au   = 1d0/au_length_aa
      ulength_from_au = au_length_aa
    case default
      stop "Invalid unit for length."
    end select

! Unit for energy
    select case(unit_energy)
    case('au','a.u.')
      uenergy_to_au   = 1d0
      uenergy_from_au = 1d0
    case('ev','eV')
      uenergy_to_au   = 1d0/au_energy_ev
      uenergy_from_au = au_energy_ev
    case default
      stop "Invalid unit for energy."
    end select

! Unit for charge
    select case(unit_charge)
    case('au','a.u.')
      ucharge_to_au   = 1d0
      ucharge_from_au = 1d0
    case default
      stop "Invalid unit for charge."
    end select

  end subroutine initialize_inputoutput_units

  subroutine dump_input_common
    use salmon_parallel
    use salmon_communication
    implicit none
    integer :: i,ierr_nml
    ierr_nml = 0

    if (comm_is_root(nproc_id_global)) then

      if(inml_calculation >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'calculation', inml_calculation
      print '("#",4X,A,"=",A)', 'calc_mode', calc_mode
      print '("#",4X,A,"=",A)', 'use_ehrenfest_md', use_ehrenfest_md
      print '("#",4X,A,"=",A)', 'use_ms_maxwell', use_ms_maxwell
      print '("#",4X,A,"=",A)', 'use_force', use_force
      print '("#",4X,A,"=",A)', 'use_geometry_opt', use_geometry_opt

      if(inml_control >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'control', inml_control
      print '("#",4X,A,"=",A)', 'restart_option', restart_option
      print '("#",4X,A,"=",I5)', 'backup_frequency', backup_frequency
      print '("#",4X,A,"=",ES12.5)', 'time_shutdown', time_shutdown
      print '("#",4X,A,"=",A)', 'sysname', sysname
      print '("#",4X,A,"=",A)', 'directory', directory

      if(inml_units >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'units', inml_units
      print '("#",4X,A,"=",A)', 'unit_time', unit_time
      print '("#",4X,A,"=",A)', 'unit_length', unit_length
      print '("#",4X,A,"=",A)', 'unit_energy', unit_energy
      print '("#",4X,A,"=",A)', 'unit_charge', unit_charge

      if(inml_parallel >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'parallel', inml_parallel
      print '("#",4X,A,"=",A)', 'domain_parallel', domain_parallel
      print '("#",4X,A,"=",I5)', 'nproc_ob', nproc_ob
      print '("#",4X,A,"=",I5)', 'nproc_domain(1)', nproc_domain(1)
      print '("#",4X,A,"=",I5)', 'nproc_domain(2)', nproc_domain(2)
      print '("#",4X,A,"=",I5)', 'nproc_domain(3)', nproc_domain(3)
      print '("#",4X,A,"=",I5)', 'nproc_domain_s(1)', nproc_domain_s(1)
      print '("#",4X,A,"=",I5)', 'nproc_domain_s(2)', nproc_domain_s(2)
      print '("#",4X,A,"=",I5)', 'nproc_domain_s(3)', nproc_domain_s(3)
      print '("#",4X,A,"=",I5)', 'num_datafiles_in', num_datafiles_in
      print '("#",4X,A,"=",I5)', 'num_datafiles_out', num_datafiles_out

      if(inml_system >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'system', inml_system
      print '("#",4X,A,"=",I1)', 'iperiodic', iperiodic
      print '("#",4X,A,"=",I1)', 'ispin', ispin
      print '("#",4X,A,"=",ES12.5)', 'al(1)', al(1)
      print '("#",4X,A,"=",ES12.5)', 'al(2)', al(2)
      print '("#",4X,A,"=",ES12.5)', 'al(3)', al(3)
      print '("#",4X,A,"=",I1)', 'isym', isym
      print '("#",4X,A,"=",A)', 'crystal_structure', crystal_structure
      print '("#",4X,A,"=",I4)', 'nstate', nstate
      print '("#",4X,A,"=",I4)', 'nelec', nelec
      print '("#",4X,A,"=",ES12.5)', 'temperature', temperature
      print '("#",4X,A,"=",I4)', 'nelem', nelem
      print '("#",4X,A,"=",I4)', 'natom', natom
      print '("#",4X,A,"=",A)', 'file_atom_coor', trim(file_atom_coor)
      print '("#",4X,A,"=",A)', 'file_atom_red_coor', trim(file_atom_red_coor)

      if(inml_pseudo >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'pseudo', inml_pseudo

      do i = 1,nelem
        print '("#",4X,A,I2,A,"=",A)', 'pseudo_file(',i,')', trim(pseudo_file(i))
        print '("#",4X,A,I2,A,"=",I4)', 'Lmax_ps(',i,')', Lmax_ps(i)
        print '("#",4X,A,I2,A"=",I4)', 'Lloc_ps(',i,')', Lloc_ps(i)
        print '("#",4X,A,I2,A"=",I4)', 'iZatom(',i,')', iZatom(i)
      end do
      print '("#",4X,A,"=",A)', 'psmask_option', psmask_option
      print '("#",4X,A,"=",ES12.5)', 'alpha_mask', alpha_mask
      print '("#",4X,A,"=",ES12.5)', 'gamma_mask', gamma_mask
      print '("#",4X,A,"=",ES12.5)', 'eta_mask', eta_mask

      if(inml_functional >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'functional', inml_functional
      print '("#",4X,A,"=",A)', 'xc', xc
      print '("#",4X,A,"=",ES12.5)', 'cval', cval

      if(inml_rgrid >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'rgrid', inml_rgrid
      print '("#",4X,A,"=",ES12.5)', 'dl(1)', dl(1)
      print '("#",4X,A,"=",ES12.5)', 'dl(2)', dl(2)
      print '("#",4X,A,"=",ES12.5)', 'dl(3)', dl(3)
      print '("#",4X,A,"=",I4)', 'num_rgrid(1)', num_rgrid(1)
      print '("#",4X,A,"=",I4)', 'num_rgrid(2)', num_rgrid(2)
      print '("#",4X,A,"=",I4)', 'num_rgrid(3)', num_rgrid(3)

      if(inml_kgrid >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I1)', 'kgrid', inml_kgrid
      print '("#",4X,A,"=",I4)', 'num_kgrid(1)', num_kgrid(1)
      print '("#",4X,A,"=",I4)', 'num_kgrid(2)', num_kgrid(2)
      print '("#",4X,A,"=",I4)', 'num_kgrid(3)', num_kgrid(3)
      print '("#",4X,A,"=",A)', 'file_kw', file_kw

      if(inml_tgrid >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'tgrid', inml_tgrid
      print '("#",4X,A,"=",I6)', 'nt', nt
      print '("#",4X,A,"=",ES12.5)', 'dt', dt

      if(inml_scf >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'scf', inml_scf
      print '("#",4X,A,"=",I3)', 'ncg', ncg
      print '("#",4X,A,"=",I3)', 'nmemory_mb', nmemory_mb
      print '("#",4X,A,"=",ES12.5)', 'alpha_mb', alpha_mb
      print '("#",4X,A,"=",A)', 'fsset_option', fsset_option
      print '("#",4X,A,"=",I3)', 'nfsset_start', nfsset_start
      print '("#",4X,A,"=",I3)', 'nfsset_every', nfsset_every
      print '("#",4X,A,"=",I3)', 'nscf', nscf
      print '("#",4X,A,"=",I3)', 'ngeometry_opt', ngeometry_opt
      print '("#",4X,A,"=",A)', 'subspace_diagonalization', subspace_diagonalization
      print '("#",4X,A,"=",A)', 'cmixing', cmixing
      print '("#",4X,A,"=",ES12.5)', 'rmixrate', rmixrate
      print '("#",4X,A,"=",A)', 'convergence', convergence
      print '("#",4X,A,"=",ES12.5)', 'threshold', threshold

      if(inml_emfield >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'emfield', inml_emfield
      print '("#",4X,A,"=",A)', 'trans_longi', trans_longi
      print '("#",4X,A,"=",A)', 'ae_shape1', ae_shape1
      print '("#",4X,A,"=",ES12.5)', 'amplitude1', amplitude1
      print '("#",4X,A,"=",ES12.5)', 'rlaser_int1', rlaser_int1
      print '("#",4X,A,"=",ES12.5)', 'pulse_tw1', pulse_tw1
      print '("#",4X,A,"=",ES12.5)', 'omega1', omega1
      print '("#",4X,A,"=",ES12.5)', 'epdir_re1(1)', epdir_re1(1)
      print '("#",4X,A,"=",ES12.5)', 'epdir_re1(2)', epdir_re1(2)
      print '("#",4X,A,"=",ES12.5)', 'epdir_re1(3)', epdir_re1(3)
      print '("#",4X,A,"=",ES12.5)', 'epdir_im1(1)', epdir_im1(1)
      print '("#",4X,A,"=",ES12.5)', 'epdir_im1(2)', epdir_im1(2)
      print '("#",4X,A,"=",ES12.5)', 'epdir_im1(3)', epdir_im1(3)
      print '("#",4X,A,"=",ES12.5)', 'phi_cep1', phi_cep1
      print '("#",4X,A,"=",A)', 'ae_shape2', ae_shape2
      print '("#",4X,A,"=",ES12.5)', 'amplitude2', amplitude2
      print '("#",4X,A,"=",ES12.5)', 'rlaser_int2', rlaser_int2
      print '("#",4X,A,"=",ES12.5)', 'pulse_tw2', pulse_tw2
      print '("#",4X,A,"=",ES12.5)', 'omega2', omega2
      print '("#",4X,A,"=",ES12.5)', 'epdir_re2(1)', epdir_re2(1)
      print '("#",4X,A,"=",ES12.5)', 'epdir_re2(2)', epdir_re2(2)
      print '("#",4X,A,"=",ES12.5)', 'epdir_re2(3)', epdir_re2(3)
      print '("#",4X,A,"=",ES12.5)', 'epdir_im2(1)', epdir_im2(1)
      print '("#",4X,A,"=",ES12.5)', 'epdir_im2(2)', epdir_im2(2)
      print '("#",4X,A,"=",ES12.5)', 'epdir_im2(3)', epdir_im2(3)
      print '("#",4X,A,"=",ES12.5)', 'phi_cep2', phi_cep2
      print '("#",4X,A,"=",ES12.5)', 't1_t2', t1_t2
      print '("#",4X,A,"=",A)', 'quadrupole', quadrupole
      print '("#",4X,A,"=",A)', 'quadrupole_pot', quadrupole_pot
      print '("#",4X,A,"=",A)', 'alocal_laser', alocal_laser
      print '("#",4X,A,"=",ES12.5)', 'rlaserbound_sta(1)', rlaserbound_sta(1)
      print '("#",4X,A,"=",ES12.5)', 'rlaserbound_sta(2)', rlaserbound_sta(2)
      print '("#",4X,A,"=",ES12.5)', 'rlaserbound_sta(3)', rlaserbound_sta(3)
      print '("#",4X,A,"=",ES12.5)', 'rlaserbound_end(1)', rlaserbound_end(1)
      print '("#",4X,A,"=",ES12.5)', 'rlaserbound_end(2)', rlaserbound_end(2)
      print '("#",4X,A,"=",ES12.5)', 'rlaserbound_end(3)', rlaserbound_end(3)

      if(inml_linear_response >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'linear_response', inml_linear_response
      print '("#",4X,A,"=",ES12.5)', 'e_impulse', e_impulse

      print '("#namelist: ",A,", status=",I3)', 'multiscale', inml_multiscale
      print '("#",4X,A,"=",A)', 'fdtddim', fdtddim
      print '("#",4X,A,"=",A)', 'twod_shape', twod_shape
      print '("#",4X,A,"=",I4)', 'nx_m', nx_m
      print '("#",4X,A,"=",I4)', 'ny_m', ny_m
      print '("#",4X,A,"=",I4)', 'nz_m', nz_m
      print '("#",4X,A,"=",ES12.5)', 'hx_m', hx_m
      print '("#",4X,A,"=",ES12.5)', 'hy_m', hy_m
      print '("#",4X,A,"=",ES12.5)', 'hz_m', hz_m
      print '("#",4X,A,"=",I4)', 'nksplit', nksplit
      print '("#",4X,A,"=",I4)', 'nxysplit', nxysplit
      print '("#",4X,A,"=",I4)', 'nxvacl_m', nxvacl_m
      print '("#",4X,A,"=",I4)', 'nxvacr_m', nxvacr_m

      if(inml_analysis >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'analysis', inml_analysis
      print '("#",4X,A,"=",A)', 'projection_option', projection_option
      print '("#",4X,A,"=",I6)', 'nenergy', nenergy
      print '("#",4X,A,"=",ES12.5)', 'de', de
      print '("#",4X,A,"=",A)', 'out_psi', out_psi
      print '("#",4X,A,"=",A)', 'out_dos', out_dos
      print '("#",4X,A,"=",A)', 'out_pdos', out_pdos
      print '("#",4X,A,"=",A)', 'out_dns', out_dns
      print '("#",4X,A,"=",A)', 'out_elf', out_elf
      print '("#",4X,A,"=",A)', 'out_dns_rt', out_dns_rt
      print '("#",4X,A,"=",I6)', 'out_dns_rt_step', out_dns_rt_step
      print '("#",4X,A,"=",A)', 'out_elf_rt', out_elf_rt
      print '("#",4X,A,"=",I6)', 'out_elf_rt_step', out_elf_rt_step
      print '("#",4X,A,"=",A)', 'out_estatic_rt', out_estatic_rt
      print '("#",4X,A,"=",I6)', 'out_estatic_rt_step', out_estatic_rt_step
      print '("#",4X,A,"=",A)', 'format3d', format3d

      if(inml_hartree >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'hartree', inml_hartree
      print '("#",4X,A,"=",I4)', 'meo', meo
      print '("#",4X,A,"=",I4)', 'num_pole_xyz(1)', num_pole_xyz(1)
      print '("#",4X,A,"=",I4)', 'num_pole_xyz(2)', num_pole_xyz(2)
      print '("#",4X,A,"=",I4)', 'num_pole_xyz(3)', num_pole_xyz(3)

      if(inml_ewald >0)ierr_nml = ierr_nml +1
      print '("#namelist: ",A,", status=",I3)', 'ewald', inml_ewald
      print '("#",4X,A,"=",I3)', 'newald', newald
      print '("#",4X,A,"=",ES12.5)', 'aewald', aewald

    end if

    call comm_bcast(ierr_nml,nproc_group_global)
    if(ierr_nml > 0)then
       if (comm_is_root(nproc_id_global)) write(*,"(I4,2x,A)")ierr_nml,'error(s) in input.'
       call end_parallel
       stop
    end if


  end subroutine dump_input_common
    
end module inputoutput



