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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
module inputoutput
  use salmon_global
  implicit none
!Physical constant
  real(8),parameter :: au_time_fs = 0.02418884326505d0
  real(8),parameter :: au_energy_ev = 27.21138505d0
  real(8),parameter :: au_length_aa = 0.52917721067d0


  integer, parameter :: fh_namelist = 901
  integer, parameter :: fh_atomic_spiecies = 902
  integer, parameter :: fh_atomic_positions = 903
  integer, parameter :: fh_reentrance = 904


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


  subroutine read_stdin(myrank)
    implicit none
    include 'mpif.h'
    integer,intent(in) :: myrank
    integer :: ierr

    integer :: cur = fh_namelist
    integer :: ret = 0
    character(100) :: buff, text


    
    if (myrank == 0) then
      open(fh_namelist, file='.namelist.tmp', status='replace')
!      open(fh_atomic_spiecies, file='.atomic_spiecies.tmp', status='replace')
      open(fh_atomic_positions, file='.atomic_positions.tmp', status='replace')
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
          if (text == '&atomic_positions') then
            cur = fh_atomic_positions
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
      close(fh_atomic_positions)
!      close(fh_atomic_spiecies)
      close(fh_reentrance)
    end if

!    call comm_sync_all()


    return
  end subroutine read_stdin

  subroutine read_input_common(myrank)
    implicit none
    include 'mpif.h'
    integer,intent(in) :: myrank
    integer :: ierr

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
      & file_atom

    namelist/pseudo/ &
      & pseudodir, &
      & Lmax_ps, &
      & Lloc_ps, &
      & iZatom, &
      & ps_format, &
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
      & quadrupole_pot

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
      & out_dns_rt, &
      & out_dns_rt_step, &
      & out_elf_rt, &
      & out_elf_rt_step, &
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

    if (myrank == 0)then
      open(fh_namelist, file='.namelist.tmp', status='old')
      read(fh_namelist, nml=units, iostat=inml_units)
      rewind(fh_namelist)
      close(fh_namelist)
    end if

    call mpi_bcast(unit_time,16,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(unit_length,16,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(unit_energy,16,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(unit_charge,16,mpi_character,0,mpi_comm_world,ierr)

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
    num_datafiles_in  = 0
    num_datafiles_out = 0
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
    file_atom          = 'none'
!! == default for &pseudo
    pseudodir     = './'
    Lmax_ps       = -1
    Lloc_ps       = -1
    iZatom        = -1
    ps_format     = 'KY'
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
    ngeometry_opt = 0
    subspace_diagonalization = 'y'
    cmixing       = 'broyden'
    rmixrate      = 0.5d0
    convergence   = 'rho'
    threshold     = 1d-6
!! == default for &emfield
    trans_longi    = 'Tr'
    ae_shape1      = ''
    amplitude1     = 0d0
    rlaser_int1    = -1d0
    pulse_tw1      = 0d0
    omega1         = 0d0
    epdir_re1      = (/1d0,0d0,0d0/)
    epdir_im1      = 0d0
    phi_cep1       = 0d0
    ae_shape2      = ''
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
!! == default for &linear_response
    e_impulse = 5d-5*uenergy_to_au/ulength_to_au*utime_to_au ! a.u.
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
    projection_option = 'no'
    nenergy           = 1000
    de                = (0.01d0/au_energy_ev)*uenergy_from_au  ! eV
    out_psi           = 'n'
    out_dos           = 'n'
    out_pdos          = 'n'
    out_dns           = 'n'
    out_dns_rt        = 'n'
    out_dns_rt_step   = 50
    out_elf_rt        = 'n'
    out_elf_rt_step   = 50
    format3d          = 'avs'
!! == default for &hartree
    meo          = 3
    num_pole_xyz = 0
!! == default for &ewald
    newald = 4
    aewald = 0.5d0

    if (myrank == 0)then
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
    call mpi_bcast(calc_mode,16,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(use_ehrenfest_md,1,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(use_ms_maxwell,1,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(use_force,1,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(use_geometry_opt,1,mpi_character,0,mpi_comm_world,ierr)
!! == bcast for &control
    call mpi_bcast(restart_option,8,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(backup_frequency,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(time_shutdown,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(sysname,256,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(directory,256,mpi_character,0,mpi_comm_world,ierr)
!! == bcast for &parallel
    call mpi_bcast(domain_parallel,1,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(nproc_ob,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nproc_domain,3,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nproc_domain_s,3,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(num_datafiles_in,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(num_datafiles_out,1,mpi_integer,0,mpi_comm_world,ierr)
!! == bcast for &system
    call mpi_bcast(iperiodic,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(ispin,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(al,3,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(isym,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(crystal_structure,32,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(nstate,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nelec,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(temperature,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(nelem,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(natom,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(file_atom,256,mpi_character,0,mpi_comm_world,ierr)
!! == bcast for &pseudo
    call mpi_bcast(pseudodir,256,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(Lmax_ps,maxMKI,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(Lloc_ps,maxMKI,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(iZatom,maxMKI,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(ps_format,16*maxMKI,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(psmask_option,1,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(alpha_mask,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(gamma_mask,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(eta_mask,1,mpi_real8,0,mpi_comm_world,ierr)
!! == bcast for &functional
    call mpi_bcast(xc,32,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(cval,1,mpi_real8,0,mpi_comm_world,ierr)
!! == bcast for &rgrid
    call mpi_bcast(dl,3,mpi_real8,0,mpi_comm_world,ierr)
    dl = dl * ulength_to_au
    call mpi_bcast(num_rgrid,3,mpi_integer,0,mpi_comm_world,ierr)
!! == bcast for &kgrid
    call mpi_bcast(num_kgrid,3,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(file_kw,256,mpi_character,0,mpi_comm_world,ierr)
!! == bcast for &kgrid
    call mpi_bcast(nt,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(dt,1,mpi_real8,0,mpi_comm_world,ierr)
    dt = dt * utime_to_au
!! == bcast for &propagation
    call mpi_bcast(n_hamil,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(propagator,16,mpi_character,0,mpi_comm_world,ierr)
!! == bcast for &scf
    call mpi_bcast(ncg,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nmemory_mb,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(alpha_mb,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(fsset_option,1,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(nfsset_start,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nfsset_every,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nscf,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(ngeometry_opt,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(subspace_diagonalization,1,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(cmixing,16,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(rmixrate,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(convergence,3,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(threshold,1,mpi_real8,0,mpi_comm_world,ierr)
!! == bcast for &emfield
    call mpi_bcast(trans_longi,2,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(ae_shape1,16,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(amplitude1,1,mpi_real8,0,mpi_comm_world,ierr)
    amplitude1 = amplitude1*(uenergy_to_au/ulength_to_au/ucharge_to_au)
    call mpi_bcast(rlaser_int1,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(pulse_tw1,1,mpi_real8,0,mpi_comm_world,ierr)
    pulse_tw1 = pulse_tw1 * utime_to_au
    call mpi_bcast(omega1,1,mpi_real8,0,mpi_comm_world,ierr)
    omega1 = omega1 * uenergy_to_au
    call mpi_bcast(epdir_re1,3,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(epdir_im1,3,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(phi_cep1,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(ae_shape2,16,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(amplitude2,1,mpi_real8,0,mpi_comm_world,ierr)
    amplitude2 = amplitude2*(uenergy_to_au/ulength_to_au/ucharge_to_au)
    call mpi_bcast(rlaser_int2,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(pulse_tw2,1,mpi_real8,0,mpi_comm_world,ierr)
    pulse_tw2 = pulse_tw2 * utime_to_au
    call mpi_bcast(omega2,1,mpi_real8,0,mpi_comm_world,ierr)
    omega2 = omega2 * uenergy_to_au
    call mpi_bcast(epdir_re2,3,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(epdir_im2,3,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(phi_cep2,1,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_bcast(t1_t2,1,mpi_real8,0,mpi_comm_world,ierr)
    t1_t2 = t1_t2 * utime_to_au
    call mpi_bcast(quadrupole,1,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(quadrupole_pot,8,mpi_character,0,mpi_comm_world,ierr)
!! == bcast for &linear_response
    call mpi_bcast(e_impulse,1,mpi_real8,0,mpi_comm_world,ierr)
    e_impulse = e_impulse *uenergy_to_au/ulength_to_au*utime_to_au
!! == bcast for &multiscale
    call mpi_bcast(fdtddim,16,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(twod_shape,16,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(nx_m,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(ny_m,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nz_m,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(hx_m,1,mpi_real8,0,mpi_comm_world,ierr)
    hx_m = hx_m * ulength_to_au
    call mpi_bcast(hy_m,1,mpi_real8,0,mpi_comm_world,ierr)
    hy_m = hy_m * ulength_to_au
    call mpi_bcast(hz_m,1,mpi_real8,0,mpi_comm_world,ierr)
    hz_m = hz_m * ulength_to_au
    call mpi_bcast(nksplit,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nxysplit,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nxvacl_m,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(nxvacr_m,1,mpi_integer,0,mpi_comm_world,ierr)
!! == bcast for &analysis
    call mpi_bcast(projection_option,2,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(nenergy,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(de,1,mpi_real8,0,mpi_comm_world,ierr)
    de = de * uenergy_to_au
    call mpi_bcast(out_psi,1,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(out_dos,1,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(out_pdos,1,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(out_dns,1,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(out_dns_rt,1,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(out_dns_rt_step,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(out_elf_rt,1,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(out_elf_rt_step,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(format3d,16,mpi_character,0,mpi_comm_world,ierr)
!! == bcast for &hartree
    call mpi_bcast(meo,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(num_pole_xyz,3,mpi_integer,0,mpi_comm_world,ierr)
!! == bcast for &ewald
    call mpi_bcast(newald,1,mpi_integer,0,mpi_comm_world,ierr)
    call mpi_bcast(aewald,1,mpi_real8,0,mpi_comm_world,ierr)

  end subroutine read_input_common

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

  subroutine dump_input_common(myrank)
    implicit none
    include 'mpif.h'
    integer,intent(in) :: myrank
    integer :: i

    if (myrank == 0) then

      print '("#namelist: ",A,", status=",I1)', 'calculation', inml_calculation
      print '("#",4X,A,"=",A)', 'calc_mode', calc_mode
      print '("#",4X,A,"=",A)', 'use_ehrenfest_md', use_ehrenfest_md
      print '("#",4X,A,"=",A)', 'use_ms_maxwell', use_ms_maxwell
      print '("#",4X,A,"=",A)', 'use_force', use_force
      print '("#",4X,A,"=",A)', 'use_geometry_opt', use_geometry_opt

      print '("#namelist: ",A,", status=",I1)', 'control', inml_control
      print '("#",4X,A,"=",A)', 'restart_option', restart_option
      print '("#",4X,A,"=",I1)', 'backup_frequency', backup_frequency
      print '("#",4X,A,"=",ES12.5)', 'time_shutdown', time_shutdown
      print '("#",4X,A,"=",A)', 'sysname', sysname
      print '("#",4X,A,"=",A)', 'directory', directory

      print '("#namelist: ",A,", status=",I1)', 'units', inml_units
      print '("#",4X,A,"=",A)', 'unit_time', unit_time
      print '("#",4X,A,"=",A)', 'unit_length', unit_length
      print '("#",4X,A,"=",A)', 'unit_energy', unit_energy
      print '("#",4X,A,"=",A)', 'unit_charge', unit_charge

      print '("#namelist: ",A,", status=",I1)', 'parallel', inml_parallel
      print '("#",4X,A,"=",A)', 'domain_parallel', domain_parallel
      print '("#",4X,A,"=",I1)', 'nproc_ob', nproc_ob
      print '("#",4X,A,"=",I1)', 'nproc_domain(1)', nproc_domain(1)
      print '("#",4X,A,"=",I1)', 'nproc_domain(2)', nproc_domain(2)
      print '("#",4X,A,"=",I1)', 'nproc_domain(3)', nproc_domain(3)
      print '("#",4X,A,"=",I1)', 'nproc_domain_s(1)', nproc_domain_s(1)
      print '("#",4X,A,"=",I1)', 'nproc_domain_s(2)', nproc_domain_s(2)
      print '("#",4X,A,"=",I1)', 'nproc_domain_s(3)', nproc_domain_s(3)
      print '("#",4X,A,"=",I1)', 'num_datafiles_in', num_datafiles_in
      print '("#",4X,A,"=",I1)', 'num_datafiles_out', num_datafiles_out

      print '("#namelist: ",A,", status=",I1)', 'system', inml_system
      print '("#",4X,A,"=",I1)', 'iperiodic', iperiodic
      print '("#",4X,A,"=",I1)', 'ispin', ispin
      print '("#",4X,A,"=",ES12.5)', 'al(1)', al(1)
      print '("#",4X,A,"=",ES12.5)', 'al(2)', al(2)
      print '("#",4X,A,"=",ES12.5)', 'al(3)', al(3)
      print '("#",4X,A,"=",I1)', 'isym', isym
      print '("#",4X,A,"=",A)', 'crystal_structure', crystal_structure
      print '("#",4X,A,"=",I1)', 'nstate', nstate
      print '("#",4X,A,"=",I1)', 'nelec', nelec
      print '("#",4X,A,"=",ES12.5)', 'temperature', temperature
      print '("#",4X,A,"=",I1)', 'nelem', nelem
      print '("#",4X,A,"=",I1)', 'natom', natom
      print '("#",4X,A,"=",A)', 'file_atom', file_atom

      print '("#namelist: ",A,", status=",I1)', 'pseudo', inml_pseudo
      print '("#",4X,A,"=",A)', 'pseudodir', pseudodir
      do i = 1,nelem
        print '("#",4X,A,"=",I1,2x,I1)', 'Lmax_ps(i)',i, Lmax_ps(i)
        print '("#",4X,A,"=",I1,2x,I1)', 'Lloc_ps(i)',i, Lloc_ps(i)
        print '("#",4X,A,"=",I1,2x,I1)', 'iZatom(i)',i, iZatom(i)
        print '("#",4X,A,"=",I1,2x,A)', 'ps_format(i)', i,ps_format(i)
      end do
      print '("#",4X,A,"=",A)', 'psmask_option', psmask_option
      print '("#",4X,A,"=",ES12.5)', 'alpha_mask', alpha_mask
      print '("#",4X,A,"=",ES12.5)', 'gamma_mask', gamma_mask
      print '("#",4X,A,"=",ES12.5)', 'eta_mask', eta_mask

      print '("#namelist: ",A,", status=",I1)', 'functional', inml_functional
      print '("#",4X,A,"=",A)', 'xc', xc
      print '("#",4X,A,"=",ES12.5)', 'cval', cval

      print '("#namelist: ",A,", status=",I1)', 'rgrid', inml_rgrid
      print '("#",4X,A,"=",ES12.5)', 'dl(1)', dl(1)
      print '("#",4X,A,"=",ES12.5)', 'dl(2)', dl(2)
      print '("#",4X,A,"=",ES12.5)', 'dl(3)', dl(3)
      print '("#",4X,A,"=",I1)', 'num_rgrid(1)', num_rgrid(1)
      print '("#",4X,A,"=",I1)', 'num_rgrid(2)', num_rgrid(2)
      print '("#",4X,A,"=",I1)', 'num_rgrid(3)', num_rgrid(3)

      print '("#namelist: ",A,", status=",I1)', 'kgrid', inml_kgrid
      print '("#",4X,A,"=",I1)', 'num_kgrid(1)', num_kgrid(1)
      print '("#",4X,A,"=",I1)', 'num_kgrid(2)', num_kgrid(2)
      print '("#",4X,A,"=",I1)', 'num_kgrid(3)', num_kgrid(3)
      print '("#",4X,A,"=",A)', 'file_kw', file_kw

      print '("#namelist: ",A,", status=",I1)', 'tgrid', inml_tgrid
      print '("#",4X,A,"=",I1)', 'nt', nt
      print '("#",4X,A,"=",ES12.5)', 'dt', dt

      print '("#namelist: ",A,", status=",I1)', 'scf', inml_scf
      print '("#",4X,A,"=",I1)', 'ncg', ncg
      print '("#",4X,A,"=",I1)', 'nmemory_mb', nmemory_mb
      print '("#",4X,A,"=",ES12.5)', 'alpha_mb', alpha_mb
      print '("#",4X,A,"=",A)', 'fsset_option', fsset_option
      print '("#",4X,A,"=",I1)', 'nfsset_start', nfsset_start
      print '("#",4X,A,"=",I1)', 'nfsset_every', nfsset_every
      print '("#",4X,A,"=",I1)', 'nscf', nscf
      print '("#",4X,A,"=",I1)', 'ngeometry_opt', ngeometry_opt
      print '("#",4X,A,"=",A)', 'subspace_diagonalization', subspace_diagonalization
      print '("#",4X,A,"=",A)', 'cmixing', cmixing
      print '("#",4X,A,"=",ES12.5)', 'rmixrate', rmixrate
      print '("#",4X,A,"=",A)', 'convergence', convergence
      print '("#",4X,A,"=",ES12.5)', 'threshold', threshold

      print '("#namelist: ",A,", status=",I1)', 'emfield', inml_emfield
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

      print '("#namelist: ",A,", status=",I1)', 'linear_response', inml_linear_response
      print '("#",4X,A,"=",ES12.5)', 'e_impulse', e_impulse

      print '("#namelist: ",A,", status=",I1)', 'multiscale', inml_multiscale
      print '("#",4X,A,"=",A)', 'fdtddim', fdtddim
      print '("#",4X,A,"=",A)', 'twod_shape', twod_shape
      print '("#",4X,A,"=",I1)', 'nx_m', nx_m
      print '("#",4X,A,"=",I1)', 'ny_m', ny_m
      print '("#",4X,A,"=",I1)', 'nz_m', nz_m
      print '("#",4X,A,"=",ES12.5)', 'hx_m', hx_m
      print '("#",4X,A,"=",ES12.5)', 'hy_m', hy_m
      print '("#",4X,A,"=",ES12.5)', 'hz_m', hz_m
      print '("#",4X,A,"=",I1)', 'nksplit', nksplit
      print '("#",4X,A,"=",I1)', 'nxysplit', nxysplit
      print '("#",4X,A,"=",I1)', 'nxvacl_m', nxvacl_m
      print '("#",4X,A,"=",I1)', 'nxvacr_m', nxvacr_m

      print '("#namelist: ",A,", status=",I1)', 'analysis', inml_analysis
      print '("#",4X,A,"=",A)', 'projection_option', projection_option
      print '("#",4X,A,"=",I1)', 'nenergy', nenergy
      print '("#",4X,A,"=",ES12.5)', 'de', de
      print '("#",4X,A,"=",A)', 'out_psi', out_psi
      print '("#",4X,A,"=",A)', 'out_dos', out_dos
      print '("#",4X,A,"=",A)', 'out_pdos', out_pdos
      print '("#",4X,A,"=",A)', 'out_dns', out_dns
      print '("#",4X,A,"=",A)', 'out_dns_rt', out_dns_rt
      print '("#",4X,A,"=",I1)', 'out_dns_rt_step', out_dns_rt_step
      print '("#",4X,A,"=",A)', 'out_elf_rt', out_elf_rt
      print '("#",4X,A,"=",I1)', 'out_elf_rt_step', out_elf_rt_step
      print '("#",4X,A,"=",A)', 'format3d', format3d

      print '("#namelist: ",A,", status=",I1)', 'hartree', inml_hartree
      print '("#",4X,A,"=",I1)', 'meo', meo
      print '("#",4X,A,"=",I1)', 'num_pole_xyz(1)', num_pole_xyz(1)
      print '("#",4X,A,"=",I1)', 'num_pole_xyz(2)', num_pole_xyz(2)
      print '("#",4X,A,"=",I1)', 'num_pole_xyz(3)', num_pole_xyz(3)

      print '("#namelist: ",A,", status=",I1)', 'ewald', inml_ewald
      print '("#",4X,A,"=",I1)', 'newald', newald
      print '("#",4X,A,"=",ES12.5)', 'aewald', aewald

    end if

  end subroutine dump_input_common
    
end module inputoutput



