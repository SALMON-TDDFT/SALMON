!
!  Copyright 2016 ARTED developers
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
module inputfile
  use inputoutput
  implicit none
  
  
contains

  
  subroutine transfer_basic_input()
    use salmon_global
    use salmon_communication, only: comm_sync_all, &
                                    comm_is_root
    use salmon_parallel, only: nproc_id_global
    use Global_Variables
    implicit none

! Be careful for backup!!
    entrance_iter = 0

    select case(calc_mode)
    case('GS_RT','gs_rt')
      iflag_calc_mode = iflag_calc_mode_gs_rt
    case('GS','gs')
      iflag_calc_mode = iflag_calc_mode_gs
    case('RT','rt')
      iflag_calc_mode = iflag_calc_mode_rt
    end select
    
!    namelist/control/ &
!    entrance_option = trim(restart_option)
!    Time_shutdown = time_shutdown
!    backup_frequency = backup_frequency
!    entrance_iter = entrance_iter
!    SYSname = sysname
!    directory = directory

!    namelist/system/ &
    functional = trim(xc)
!    cval = cval
    Sym = isym
!    crystal_structure = crystal_structure
    NB = nstate
!    Nelec = nelec
!    ext_field ! this variable is removed
!    MD_option   ! this variable is replaced by use_ehrenfest_md
!    AD_RHO ! this variable is replaced by projection_option

!    namelist/rgrid/ &
!            Nd,   ! Removed from input
    if(sum(abs(num_rgrid)) /= 0 .and. sum(abs(dl)) /= 0d0)then
      call err_finalize('Error: [num_rgrid] and [dl] are incompatible input parameters.')
    else if(sum(abs(num_rgrid)) /= 0 .and. sum(abs(dl)) == 0d0)then
      NLx = num_rgrid(1)
      NLy = num_rgrid(2)
      NLz = num_rgrid(3)
    else if(sum(abs(num_rgrid)) == 0d0 .and. sum(abs(dl)) /= 0d0)then
      NLx = nint(al(1)/dl(1))
      NLy = nint(al(2)/dl(2))
      NLz = nint(al(3)/dl(3))
      if(comm_is_root(nproc_id_global))then
        write(*,"(A)")"Grid number [num_rgrid(3)] is computed from Grid spacing [dl(3)]."
        write(*,"(A,2x,I7)")"num_rgrid(1)=",NLx
        write(*,"(A,2x,I7)")"num_rgrid(2)=",NLy
        write(*,"(A,2x,I7)")"num_rgrid(3)=",NLz
        write(*,"(A)")"Actual grid spacing becomes:"
        write(*,"(A,2x,e16.6e3,A)")"dl(1)=",al(1)/NLx," [Bohr]"
        write(*,"(A,2x,e16.6e3,A)")"dl(2)=",al(2)/NLy," [Bohr]"
        write(*,"(A,2x,e16.6e3,A)")"dl(3)=",al(3)/NLz," [Bohr]"
        write(*,"(A)")"Warning: actual grid spacing can be different from input"
        write(*,"(A)")"so that the lattice parameter is divisible by the spacing."
      end if
    else
      call err_finalize('Error: [num_rgrid] or [dl] should be specified in input.')
    end if

!! Comment: S.A.S 2017/06/04
!! conversion from dl to num_rgrid(3) should be implemented

!    namelist/kgrid/ &
     NKx = num_kgrid(1)
     NKy = num_kgrid(2)
     NKz = num_kgrid(3)
!     file_kw =  file_kw

!    namelist/tstep/ &
!     Nt = nt
!     dt = dt

!    namelist/pseudo/ &
!            & PSmask_option, &
!            & alpha_mask, &
!            & gamma_mask, &
!            & eta_mask

!    namelist/electrons/ &
!            & NEwald, &
!            & aEwald, &
     KbTev = temperature * au_energy_ev
!             Ncg
!             Nmemory_MB
!             alpha_MB
!             FSset_option
!             NFSset_start
!             NFSset_every
!             Nscf

!    namelist/incident/ &
!            & Longi_Trans, & ! This variable is replaced by trans_longi
!     dAc = e_impulse
!     AE_shape = trim(ae_shape1)
!     IWcm2_1 = rlaser_int1
!     tpulsefs_1 = pulse_tw1*au_time_fs
!     omegaev_1 = omega1*au_energy_ev
!     phi_CEP_1 = phi_cep1
!     Epdir_1 = epdir_re1
!     IWcm2_2 = rlaser_int2
!     tpulsefs_2 = pulse_tw2*au_time_fs
!     omegaev_2 = omega2*au_energy_ev
!     phi_CEP_2 = phi_cep1
!     Epdir_2 = epdir_re1
!     T1_T2fs = t1_t2*au_time_fs
!
!    namelist/propagation/ &
!            & propagator
     dAc = e_impulse
!    namelist/response/ &

     Nomega = nenergy
     domega = de

!    namelist/multiscale/ &
!     FDTDdim
!     TwoD_shape
!     NX_m
!     NY_m
!     HX_m
!     HY_m
!     NKsplit
!     NXYsplit
!     NXvacL_m
!     NXvacR_m

 !   namelist/group_atom/ &
     MI = natom; NI = natom
     MKI = nelem; NE = nelem
!     iZatom
!     ipsfileform
!     Lmax_ps
!     Lloc_ps
!     ps_format
      
    call comm_sync_all()
    return
  end subroutine transfer_basic_input


  subroutine transfer_atomic_data()
    use salmon_global
    use salmon_parallel, only: nproc_group_global, nproc_id_global
    use salmon_communication, only: comm_is_root, comm_bcast, comm_sync_all
    use Global_Variables, only: NE, Zatom, Lref
    implicit none
    integer :: i
!Note; NI = MI, NE = MKI

    allocate(Zatom(NE), Lref(NE))

    if (comm_is_root(nproc_id_global)) then
       do i=1, NE
          Zatom(i) = iZatom(i)
          Lref(i) = Lloc_ps(i)
       end do
    end if

    call comm_bcast(Zatom, nproc_group_global)
    call comm_bcast(Lref, nproc_group_global)
    call comm_bcast(ipsfileform, nproc_group_global)

!    if (comm_is_root(nproc_id_global)) then
!      open(fh_atomic_spiecies, file='.atomic_spiecies.tmp', status='old')
!      do i=1, NE
!        read(fh_atomic_spiecies, *) index, zatom_tmp, lref_tmp
!        if (i == index) then
!          Zatom(i) = Zatom_tmp
!          Lref(i) = Lref_tmp
!        else
!          call err_finalize('atomic_spiecies is not ordered')
!        end if
!      end do
!      close(fh_atomic_positions)
!    end if
!    call comm_bcast(Zatom, nproc_group_global)
!    call comm_bcast(Lref, nproc_group_global)

    call comm_sync_all()
    return
  end subroutine transfer_atomic_data


  subroutine transfer_input()
    implicit none
    
!    call extract_stdin()
    call transfer_basic_input()
    call transfer_atomic_data()

    return
  end subroutine transfer_input
  
  
  
end module inputfile
