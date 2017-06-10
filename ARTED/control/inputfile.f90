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

  
  subroutine read_namelist()
    use salmon_global
    use communication, only: comm_is_root, comm_bcast, comm_sync_all, proc_group
    use Global_Variables
    implicit none
    character(100) :: file_atoms_coo

! Be careful for backup!!
    entrance_iter = 0

    
!    namelist/control/ &
    entrance_option = trim(restart_option)
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
    NLx = num_rgrid(1)
    NLy = num_rgrid(2)
    NLz = num_rgrid(3)
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
     file_atoms_coo = trim(file_atom)
!     Lmax_ps
!     Lloc_ps
!     ps_format
      
    call comm_sync_all()
    return
  end subroutine read_namelist


  subroutine read_atomic_spiecies()
    use salmon_global
    use communication, only: comm_is_root, comm_bcast, comm_sync_all, proc_group
    use Global_Variables, only: NE, Zatom, Lref
    implicit none
    integer :: i, index, Zatom_tmp, Lref_tmp
!Note; NI = MI, NE = MKI

    allocate(Zatom(NE), Lref(NE))

    if (comm_is_root()) then
       do i=1, NE
          Zatom(i) = iZatom(i)
          Lref(i) = Lloc_ps(i)
       end do
    end if

    call comm_bcast(Zatom, proc_group(1))
    call comm_bcast(Lref, proc_group(1))
    call comm_bcast(ipsfileform, proc_group(1))

!    if (comm_is_root()) then
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
!    call comm_bcast(Zatom, proc_group(1))
!    call comm_bcast(Lref, proc_group(1))

    call comm_sync_all()
    return
  end subroutine


  subroutine read_atomic_positions()
    use communication, only: comm_is_root, comm_bcast, comm_sync_all, proc_group
    use Global_Variables, only: NI, Rion, Kion
    implicit none
    integer :: i, Kion_tmp
    real(8) :: Rion_tmp(3)
    character(32) :: char_atom
    
    allocate(Rion(3,NI), Kion(NI))
    
    if (comm_is_root()) then
      open(fh_atomic_positions, file='.atomic_positions.tmp', status='old')
      do i=1, NI
        read(fh_atomic_positions, *) char_atom, Rion_tmp, Kion_tmp
        Rion(:,i) = Rion_tmp
        Kion(i) = Kion_tmp
      end do
      close(fh_atomic_positions)
    end if
    call comm_bcast(Rion, proc_group(1))
    call comm_bcast(Kion, proc_group(1))
    call comm_sync_all()
    return
  end subroutine read_atomic_positions





  subroutine read_input()
    use global_variables, only: entrance_option
    implicit none
    
!    call extract_stdin()
    call read_namelist()
    if(entrance_option == 'reentrance')return
    call read_atomic_spiecies()
    call read_atomic_positions()
    call dump_inputdata()

    return
  end subroutine read_input
  
  
  
  subroutine dump_inputdata()
    use communication, only: comm_is_root, comm_sync_all
    use Global_Variables
    implicit none
    integer :: i
  
    if (comm_is_root()) then

      print *, "#section: atomic_positions"
      do i=1, NE
        print '("#",4X,I1,",Zatom=",I1,",Lref=",I1)', i, Zatom(i), Lref(i)
      end do
      
      print *, "#section: atomic_positions"
      do i=1, NI
        print '("#",4X,I1,",Rion=",3ES12.5,",Kion=",I1)', i, Rion(:,i), Kion(i)
      end do
    end if
    call comm_sync_all()
    return
  end subroutine dump_inputdata
  
  
  
end module inputfile
