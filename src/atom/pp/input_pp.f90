!
!  Copyright 2018 SALMON developers
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
subroutine input_pp(pp,hx,hy,hz)
  use salmon_pp,only : pp_info
  use salmon_global,only : pseudo_file
  use salmon_global,only : n_Yabana_Bertsch_psformat,n_ABINIT_psformat&
    &,n_ABINITFHI_psformat,n_FHI_psformat,ps_format,nelem,directory, &
    & psmask_option
  use salmon_parallel, only: nproc_group_global, nproc_id_global
  use salmon_communication, only: comm_bcast, comm_is_root
  implicit none
  real(8),parameter :: pi=3.141592653589793d0 ! copied from salmon_math
  type(pp_info) :: pp
  real(8),parameter :: Eps0=1d-10
  real(8),intent(in) :: hx,hy,hz
  integer :: ik,l,i
  real(8) :: rrc(0:pp%lmax0)
  real(8) :: r1
  real(8),allocatable :: rhor_nlcc(:,:)   !zero in radial index for taking derivative
  character(256) :: ps_file
  logical,allocatable :: flag_nlcc_element(:)

  allocate(rhor_nlcc(0:pp%nrmax0,0:2))

! Nonlinear core correction
  allocate(flag_nlcc_element(nelem)); flag_nlcc_element(:) = .false. ; pp%flag_nlcc = .false.

!      ps_file=trim(directory)//trim(pp%atom_symbol(ik))//trim(ps_postfix)

  if (comm_is_root(nproc_id_global)) then

    do ik=1,nelem
         
      ps_file = trim(pseudo_file(ik))

      select case (ps_format(ik))
      case('KY')
        call read_ps_ky(pp,rrc,ik,ps_file)
      case('ABINIT')
        call read_ps_abinit(pp,rrc,ik,ps_file)
      case('ABINITFHI')
        call read_ps_abinitfhi(pp,rrc,rhor_nlcc,flag_nlcc_element,ik,ps_file)
      case('FHI')
        call read_ps_fhi(pp,rrc,ik,ps_file)
!      case('ATOM')      ; call read_ps_ATOM
      case default ; stop 'Unprepared ps_format is required input_pseudopotential_YS'
      end select

! Set meaning domain in the arrays 
      pp%rps(ik)=maxval(rrc(0:pp%mlps(ik)))
      do i=1,pp%nrmax
        if(pp%rad(i,ik).gt.pp%rps(ik)) exit
      enddo
      pp%nrps(ik)=i
      if(pp%nrps(ik).ge.pp%nrmax) stop 'NRps>Nrmax at input_pseudopotential_YS'
      pp%nrloc(ik)=pp%nrps(ik)
      pp%rloc(ik)=pp%rps(ik)
      pp%radnl(:,ik)=pp%rad(:,ik)

      do l=0,pp%mlps(ik)
        pp%anorm(l,ik) = 0.d0
        do i=1,pp%mr(ik)-1
          r1 = pp%rad(i+1,ik)-pp%rad(i,ik)
          pp%anorm(l,ik) = pp%anorm(l,ik)  &
               + (pp%upp(i,l)**2*(pp%vpp(i,l)-pp%vpp(i,pp%lref(ik)))+pp%upp(i-1,l)**2*(pp%vpp(i-1,l)-pp%vpp(i-1,pp%lref(ik))))*r1
        end do
        pp%anorm(l,ik) = 0.5d0*pp%anorm(l,ik)
        pp%inorm(l,ik)=+1
        if(pp%anorm(l,ik).lt.0.d0) then
          pp%anorm(l,ik)=-pp%anorm(l,ik)
          pp%inorm(l,ik)=-1
        endif
        if(abs(pp%anorm(l,ik)).lt.Eps0)pp%inorm(l,ik)=0
        pp%anorm(l,ik)=sqrt(pp%anorm(l,ik))
      enddo

      write(*,*) '===================pseudopotential data==================='
      write(*,*) 'ik ,atom_symbol=',ik, pp%atom_symbol(ik)
      write(*,*) 'ps_format =',ps_format(ik)
      write(*,*) 'ps_file =',trim(ps_file)
      write(*,*) 'Zps(ik), Mlps(ik) =',pp%zps(ik), pp%mlps(ik)
      write(*,*) 'Rps(ik), NRps(ik) =',pp%rps(ik), pp%nrps(ik)
      write(*,*) 'Lref(ik) =',pp%lref(ik)
      write(*,*) 'anorm(ik,l) =',(real(pp%anorm(l,ik)),l=0,pp%mlps(ik))
      write(*,*) 'inorm(ik,l) =',(pp%inorm(l,ik),l=0,pp%mlps(ik))
      write(*,*) 'Mass(ik) =',pp%rmass(ik)
      write(*,*) 'flag_nlcc_element(ik) =',flag_nlcc_element(ik)
      write(*,*) '=========================================================='

      if (PSmask_option == 'y') then
        call making_ps_with_masking(pp,hx,hy,hz,ik, &
                                    rhor_nlcc,flag_nlcc_element)
        write(*,*) 'Following quantities are modified by masking procedure'
        write(*,*) 'Rps(ik), NRps(ik) =',pp%rps(ik), pp%nrps(ik)
        write(*,*) 'anorm(ik,l) =',(pp%anorm(l,ik),l=0,pp%mlps(ik))
        write(*,*) 'inorm(ik,l) =',(pp%inorm(l,ik),l=0,pp%mlps(ik))
      else if (PSmask_option == 'n') then
        call making_ps_without_masking(pp,ik,flag_nlcc_element,rhor_nlcc)

      else
        stop 'Wrong PSmask_option at input_pseudopotential_YS'
      end if

      pp%upp_f(:,:,ik)=pp%upp(:,:)
      pp%vpp_f(:,:,ik)=pp%vpp(:,:)

      open(4,file=trim(directory)//"PS_"//trim(pp%atom_symbol(ik))//"_"//trim(ps_format(ik))//"_"//trim(PSmask_option)//".dat")
      write(4,*) "# Mr=",pp%mr(ik)
      write(4,*) "# Rps(ik), NRps(ik)",pp%rps(ik), pp%nrps(ik)
      write(4,*) "# Mlps(ik), Lref(ik) =",pp%mlps(ik), pp%lref(ik)
      do i=1,pp%nrps(ik)
        write(4,'(30e21.12)') pp%rad(i,ik),(pp%udvtbl(i,l,ik),l=0,pp%mlps(ik)),(pp%dudvtbl(i,l,ik),l=0,pp%mlps(ik))
      end do
      close(4)

    enddo
  endif

  call comm_bcast(pp%mr,nproc_group_global)
  call comm_bcast(pp%zps,nproc_group_global)
  call comm_bcast(pp%mlps,nproc_group_global)
  call comm_bcast(pp%rps,nproc_group_global)
  call comm_bcast(pp%nrps,nproc_group_global)
  call comm_bcast(pp%nrloc,nproc_group_global)
  call comm_bcast(pp%rloc,nproc_group_global)
  call comm_bcast(pp%anorm,nproc_group_global)
  call comm_bcast(pp%inorm,nproc_group_global)
  call comm_bcast(pp%rad,nproc_group_global)
  call comm_bcast(pp%radnl,nproc_group_global)
  call comm_bcast(pp%vloctbl,nproc_group_global)
  call comm_bcast(pp%dvloctbl,nproc_group_global)
  call comm_bcast(pp%udvtbl,nproc_group_global)
  call comm_bcast(pp%dudvtbl,nproc_group_global)
  call comm_bcast(pp%upp_f,nproc_group_global)
  call comm_bcast(pp%vpp_f,nproc_group_global)
  call comm_bcast(pp%rho_nlcc_tbl,nproc_group_global)
  call comm_bcast(pp%tau_nlcc_tbl,nproc_group_global)
  call comm_bcast(pp%flag_nlcc,nproc_group_global)
  if(comm_is_root(nproc_id_global)) write(*,*)"flag_nlcc = ",pp%flag_nlcc

  return
end subroutine input_pp
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine read_ps_ky(pp,rrc,ik,ps_file)
  use salmon_pp,only : pp_info
  use salmon_global,only : Lmax_ps
  use inputoutput,only : au_length_aa, au_energy_ev
  implicit none
  type(pp_info),intent(inout) :: pp
!argument
  integer,intent(in) :: ik
  real(8),intent(out) :: rrc(0:pp%lmax0)
  character(256),intent(in) :: ps_file
!local variable
  integer :: l,i,irPC
  real(8) :: step,rPC,r,rhopp(0:pp%nrmax0),rzps

  open(4,file=ps_file,status='old')
  read(4,*) pp%mr(ik),step,pp%mlps(ik),rzps
  pp%zps(ik)=int(rzps+1d-10)
  if(pp%mr(ik).gt.pp%nrmax0) stop 'Mr>Nrmax0 at Read_PS_KY'
  if(pp%mlps(ik).gt.pp%lmax0) stop 'Mlps(ik)>Lmax0 at Read_PS_KY'
  if(pp%mlps(ik).gt.pp%lmax) stop 'Mlps(ik)>Lmax at Read_PS_KY'
  read(4,*) irPC,(rrc(l),l=0,pp%mlps(ik))
  rPC=real(irPC) !Radius for partial core correction: not working in this version
  do i=0,pp%mr(ik)
    read(4,*) r,rhopp(i),(pp%vpp(i,l),l=0,pp%mlps(ik))
  end do
  do i=0,pp%mr(ik)
    read(4,*) r,(pp%upp(i,l),l=0,pp%mlps(ik))
  end do
  close(4)

! change to atomic unit
  step=step/au_length_aa
  rrc(0:pp%mlps(ik))=rrc(0:pp%mlps(ik))/au_length_aa
  pp%vpp(0:pp%mr(ik),0:pp%mlps(ik))=pp%vpp(0:pp%mr(ik),0:pp%mlps(ik))/au_energy_ev
  pp%upp(0:pp%mr(ik),0:pp%mlps(ik))=pp%upp(0:pp%mr(ik),0:pp%mlps(ik))*sqrt(au_length_aa)

  do i=1,pp%nrmax
    pp%rad(i,ik)=(i-1)*step
  enddo

  if(Lmax_ps(ik) >= 0)pp%mlps(ik) = Lmax_ps(ik) ! Maximum angular momentum given by input

  return
end subroutine read_ps_KY
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine read_ps_abinit(pp,rrc,ik,ps_file)
  use salmon_pp,only : pp_info
  use salmon_global,only : Lmax_ps
!See http://www.abinit.org/downloads/psp-links/psp-links/lda_tm
  implicit none
  type(pp_info),intent(inout) :: pp
!argument
  integer,intent(in) :: ik
  real(8),intent(out) :: rrc(0:pp%lmax0)
  character(256),intent(in) :: ps_file
!local variable
  integer :: i
  real(8) :: rzps
  integer :: ll
  real(8) :: zatom, pspdat,pspcod,pspxc,lmaxabinit,lloc,mmax,r2well,l
  real(8) :: e99_0,e99_9,nproj,rcpsp,rms,ekb1,ekb2,epsatm,rchrg,fchrg,qchrg
  character(1) :: dummy_text

  open(4,file=ps_file,status='old')
  read(4,*) dummy_text
  read(4,*) zatom, pp%zion, pspdat
  rzps = pp%zion
  pp%zps(ik)=int(rzps+1d-10)
  read(4,*) pspcod,pspxc,lmaxabinit,lloc,mmax,r2well
  pp%mlps(ik)=lmaxabinit
  if(lloc .ne. pp%lref(ik)) write(*,*) "Warning! Lref(ik=",ik,") is different from intended one in ",ps_file
  pp%mr(ik) = mmax - 1
  if(pp%mr(ik).gt.pp%nrmax0) stop 'Mr>Nrmax0 at Read_PS_ABINIT'
  if(pp%mlps(ik).gt.pp%lmax0) stop 'Mlps(ik)>Lmax0 at Read_PS_ABINIT'
  if(pp%mlps(ik).gt.pp%lmax) stop 'Mlps(ik)>Lmax at Read_PS_ABINIT'
  do ll=0,pp%mlps(ik)
    read(4,*) l,e99_0,e99_9,nproj,rcpsp
    read(4,*) rms,ekb1,ekb2,epsatm
    rrc(ll) = rcpsp
  end do
  read(4,*) rchrg,fchrg,qchrg
  do ll=0,pp%mlps(ik)
    read(4,*) dummy_text
    do i=1,(pp%mr(ik)+1)/3
      read(4,*) pp%vpp(3*(i-1),ll),pp%vpp(3*(i-1)+1,ll),pp%vpp(3*(i-1)+2,ll)
    end do
  end do
  do ll=0,pp%mlps(ik)
    read(4,*) dummy_text
    do i=1,(pp%mr(ik)+1)/3
      read(4,*) pp%upp(3*(i-1),ll),pp%upp(3*(i-1)+1,ll),pp%upp(3*(i-1)+2,ll)
    end do
  end do
  close(4)

  do i=0,pp%nrmax-1
    pp%rad(i+1,ik) = 1.0d2*(dble(i)/dble(mmax-1)+1.0d-2)**5 - 1.0d-8
  end do

  if(Lmax_ps(ik) >= 0)pp%mlps(ik) = Lmax_ps(ik) ! Maximum angular momentum given by input

  return
end subroutine read_ps_abinit
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine read_ps_abinitfhi(pp,rrc,rhor_nlcc,flag_nlcc_element,ik,ps_file)
!This is for  FHI pseudopotential listed in abinit web page and not for original FHI98PP.
!See http://www.abinit.org/downloads/psp-links/lda_fhi
  use salmon_pp,only : pp_info
  use salmon_global,only : nelem, Lmax_ps
  implicit none
  real(8),parameter :: pi=3.141592653589793d0 ! copied from salmon_math
  type(pp_info),intent(inout) :: pp
!argument
  integer,intent(in) :: ik
  real(8),intent(out) :: rrc(0:pp%lmax0)
  character(256),intent(in) :: ps_file
  real(8),intent(out) :: rhor_nlcc(0:pp%nrmax0,0:2)
  logical,intent(inout) :: flag_nlcc_element(nelem)
!local variable
  character(50) :: temptext
  integer :: i,j
  real(8) :: step,rzps,dummy
  integer :: mr_l(0:pp%lmax0),l,ll
  real(8) :: step_l(0:pp%lmax),rrc_mat(0:pp%lmax,0:pp%lmax)

  rrc_mat(0:pp%lmax,0:pp%lmax)=-1.d0

  open(4,file=ps_file,status='old')
  write(*,*) '===================Header of ABINITFHI pseudo potential==================='
  do i=1,7
    read(4,'(a)') temptext
    write(*,*) temptext
  end do
  write(*,*) '===================Header of ABINITFHI pseudo potential==================='
  read(4,*) rzps,pp%mlps(ik)
  pp%zps(ik)=int(rzps+1d-10)
  pp%mlps(ik) = pp%mlps(ik)-1
  if(pp%mlps(ik).gt.pp%lmax0) stop 'Mlps(ik)>Lmax0 at Read_PS_FHI'
  if(pp%mlps(ik).gt.pp%lmax) stop 'Mlps(ik)>Lmax at Read_PS_FHI'
  do i=1,10
     read(4,*) dummy
  end do
  do l=0,pp%mlps(ik)
    read(4,*) mr_l(l),step_l(l)
    if(mr_l(l).gt.pp%nrmax0) stop 'Mr>Nrmax0 at Read_PS_FHI'
    do i=1,mr_l(l)
      read(4,*) j,pp%rad(i+1,ik),pp%upp(i,l),pp%vpp(i,l) !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate
    end do
    pp%rad(1,ik)=0.d0
    pp%upp(0,l)=0.d0
    pp%vpp(0,l)=pp%vpp(1,l)-(pp%vpp(2,l)-pp%vpp(1,l))/(pp%rad(3,ik)-pp%rad(2,ik))*(pp%rad(2,ik))
  end do
!Nonlinear core-correction
  do i=1,mr_l(0)
    read(4,*,end=940) pp%rad(i+1,ik),rhor_nlcc(i,0),rhor_nlcc(i,1),rhor_nlcc(i,2)
  end do
  rhor_nlcc(0,:)=rhor_nlcc(1,:)-(rhor_nlcc(2,:)-rhor_nlcc(1,:)) &
    /(pp%rad(3,ik)-pp%rad(2,ik))*(pp%rad(2,ik))
  rhor_nlcc = rhor_nlcc/(4d0*pi)
  flag_nlcc_element(ik) = .true.
!Nonlinear core-correction


940 close(4)
  

  if(minval(mr_l(0:pp%mlps(ik))).ne.maxval(mr_l(0:pp%mlps(ik)))) then
    stop 'Mr are diffrent at Read_PS_FHI'
  else 
    pp%mr(ik) = minval(mr_l(0:pp%mlps(ik)))
  end if
  if((maxval(step_l(0:pp%mlps(ik)))-minval(step_l(0:pp%mlps(ik)))).ge.1.d-14) then
    stop 'step are different at Read_PS_FHI'
  else 
    step = minval(step_l(0:pp%mlps(ik)))
  end if

  do i=pp%mr(ik)+1,pp%nrmax-1
    pp%rad(i+1,ik) = pp%rad(i,ik)*step
  end do

  do l=0,pp%mlps(ik)
    do ll=0,pp%mlps(ik)
      do i=pp%mr(ik),1,-1
        if(abs(pp%vpp(i,l)-pp%vpp(i,ll)).gt.1.d-10) then
          rrc_mat(l,ll) = pp%rad(i+1+1,ik)
          exit
        end if
      end do
    end do
    rrc(l)=maxval(rrc_mat(l,:))
  end do

  if(Lmax_ps(ik) >= 0)pp%mlps(ik) = Lmax_ps(ik) ! Maximum angular momentum given by input

  return
end subroutine read_ps_abinitfhi
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine read_ps_fhi(pp,rrc,ik,ps_file)
!This is for original FHI98PP and not for FHI pseudopotential listed in abinit web page
!See http://th.fhi-berlin.mpg.de/th/fhi98md/fhi98PP/
  use salmon_pp,only : pp_info
  use salmon_global,only : Lmax_ps
  implicit none
  type(pp_info),intent(inout) :: pp
!argument
  integer,intent(in) :: ik
  real(8),intent(out) :: rrc(0:pp%lmax0)
  character(256),intent(in) :: ps_file
!local variable
  integer :: i,j
  real(8) :: step,rzps,dummy
  integer :: mr_l(0:pp%lmax0),l,ll
  real(8) :: step_l(0:pp%lmax),rrc_mat(0:pp%lmax,0:pp%lmax)

  rrc_mat(0:pp%lmax,0:pp%lmax)=-1.d0

  open(4,file=ps_file,status='old')
  read(4,*) rzps,pp%mlps(ik)
  pp%zps(ik)=int(rzps+1d-10)
  pp%mlps(ik) = pp%mlps(ik)-1
  if(pp%mlps(ik).gt.pp%lmax0) stop 'Mlps(ik)>Lmax0 at Read_PS_FHI'
  if(pp%mlps(ik).gt.pp%lmax) stop 'Mlps(ik)>Lmax at Read_PS_FHI'
  do i=1,10
     read(4,*) dummy
  end do
  do l=0,pp%mlps(ik)
    read(4,*) mr_l(l),step_l(l)
    if(mr_l(l).gt.pp%nrmax0) stop 'Mr>Nrmax0 at Read_PS_FHI'
    do i=1,mr_l(l)
      read(4,*) j,pp%rad(i+1,ik),pp%upp(i,l),pp%vpp(i,l) !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate
    end do
    pp%rad(1,ik)=0.d0
    pp%upp(0,l)=0.d0
    pp%vpp(0,l)=pp%vpp(1,l)-(pp%vpp(2,l)-pp%vpp(1,l))/(pp%rad(3,ik)-pp%rad(2,ik))*(pp%rad(2,ik))
  end do
  close(4)

  if(minval(mr_l(0:pp%mlps(ik))).ne.maxval(mr_l(0:pp%mlps(ik)))) then
    stop 'Mr are diffrent at Read_PS_FHI'
  else 
    pp%mr(ik) = minval(mr_l(0:pp%mlps(ik)))
  end if
  if((maxval(step_l(0:pp%mlps(ik)))-minval(step_l(0:pp%mlps(ik)))).ge.1.d-14) then
    stop 'step are different at Read_PS_FHI'
  else 
    step = minval(step_l(0:pp%mlps(ik)))
  end if

  do i=pp%mr(ik)+1,pp%nrmax
    pp%rad(i+1,ik) = pp%rad(i,ik)*step
  end do

  do l=0,pp%mlps(ik)
    do ll=0,pp%mlps(ik)
      do i=pp%mr(ik),1,-1
        if(abs(pp%vpp(i,l)-pp%vpp(i,ll)).gt.1.d-10) then
          rrc_mat(l,ll) = pp%rad(i+1+1,ik)
          exit
        end if
      end do
    end do
    rrc(l)=maxval(rrc_mat(l,:))
  end do

  if(Lmax_ps(ik) >= 0)pp%mlps(ik) = Lmax_ps(ik) ! Maximum angular momentum given by input

  return
end subroutine read_ps_fhi
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
!    subroutine read_ps_ATOM !.psf format created by ATOM for SIESTA
!      implicit none
!      return
!    end subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine making_ps_with_masking(pp,hx,hy,hz,ik, &
                          rhor_nlcc,flag_nlcc_element)
  use salmon_pp,only : pp_info
  use salmon_global, only: ps_format, nelem, alpha_mask, eta_mask
  implicit none
  real(8),parameter :: pi=3.141592653589793d0 ! copied from salmon_math
  type(pp_info),intent(inout) :: pp
  integer,intent(in) :: ik
  real(8),intent(in) :: hx,hy,hz
  real(8),intent(in) :: rhor_nlcc(0:pp%nrmax0,0:2)
  logical,intent(in) :: flag_nlcc_element(nelem)
  real(8) :: eta
  integer :: ncounter
  real(8) :: uvpp(0:pp%nrmax0,0:pp%lmax0),duvpp(0:pp%nrmax0,0:pp%lmax0)
  real(8) :: vpploc(0:pp%nrmax0),dvpploc(0:pp%nrmax0)
  real(8) :: grid_function(0:pp%nrmax0)
  integer :: i,l
  real(8) :: r1,r2,r3,r4

  ncounter = 0
  do i=0,pp%mr(ik)
    if (pp%rad(i+1,ik) > dble(ncounter+1.d0)*max(Hx,Hy,Hz)) then
      ncounter = ncounter + 1
    end if
    if (ncounter/2*2 == ncounter) then 
      grid_function(i) = 1.d0
    else
      grid_function(i) = 0.d0
    end if
  end do

  vpploc(:) = pp%vpp(:,pp%lref(ik))
  do l=0,pp%mlps(ik)
    do i=0,pp%mr(ik)
       uvpp(i,l)=pp%upp(i,l)*(pp%vpp(i,l)-pp%vpp(i,pp%lref(ik)))
    end do
    do i=1,pp%mr(ik)-1
      r1 = pp%rad(i+1,ik)-pp%rad(i,ik)
      r2 = pp%rad(i+1,ik)-pp%rad(i+2,ik)
      r3 = pp%rad(i+2,ik)-pp%rad(i,ik)
      r4 = r1/r2
      duvpp(i,l)=(r4+1.d0)*(uvpp(i,l)-uvpp(i-1,l))/r1-(uvpp(i+1,l)-uvpp(i-1,l))/r3*r4
      dvpploc(i)=(r4+1.d0)*(vpploc(i)-vpploc(i-1))/r1-(vpploc(i+1)-vpploc(i-1))/r3*r4
    end do
    duvpp(0,l)=2.d0*duvpp(1,l)-duvpp(2,l)
    duvpp(pp%mr(ik),l)=2.d0*duvpp(pp%mr(ik)-1,l)-duvpp(pp%mr(ik)-2,l)
    dvpploc(0)=dvpploc(1)-(dvpploc(2)-dvpploc(1))/(pp%rad(3,ik)-pp%rad(2,ik))*(pp%rad(2,ik)-pp%rad(1,ik))
    dvpploc(pp%mr(ik))=dvpploc(pp%mr(ik)-1)+(dvpploc(pp%mr(ik)-1)-dvpploc(pp%mr(ik)-2))/(pp%rad(pp%mr(ik),ik)  &
                      -pp%rad(pp%mr(ik)-1,ik))*(pp%rad(pp%mr(ik)+1,ik)-pp%rad(pp%mr(ik),ik))
  end do

  open(4,file="PSbeforemask_"//trim(pp%atom_symbol(ik))//"_"//trim(ps_format(ik))//".dat")
  write(4,*) "# Mr =",pp%mr(ik)
  write(4,*) "# Rps(ik), NRps(ik)",pp%rps(ik), pp%nrps(ik)
  write(4,*) "# Mlps(ik), Lref(ik) =",pp%mlps(ik), pp%lref(ik)
  do i=0,pp%mr(ik)
    write(4,'(30e21.12)') pp%rad(i+1,ik),(uvpp(i,l),l=0,pp%mlps(ik)),(duvpp(i,l),l=0,pp%mlps(ik)),  &
                          vpploc(i),dvpploc(i),grid_function(i)
  end do
  close(4)

  call ps_masking(pp,uvpp,duvpp,ik,hx,hy,hz)

  open(4,file="PSaftermask_"//trim(pp%atom_symbol(ik))//"_"//trim(ps_format(ik))//".dat")
  write(4,*) "# Mr =",pp%mr(ik)
  write(4,*) "# Rps(ik), NRps(ik)",pp%rps(ik), pp%nrps(ik)
  write(4,*) "# Mlps(ik), Lref(ik) =",pp%mlps(ik), pp%lref(ik)
  eta = alpha_mask*Pi*pp%rps(ik)/max(Hx,Hy,Hz)
  write(4,*) "# eta_mask, eta =",eta_mask,eta
  do i=0,pp%mr(ik)
    write(4,'(30e21.12)') pp%rad(i+1,ik),(uvpp(i,l),l=0,pp%mlps(ik)),(duvpp(i,l),l=0,pp%mlps(ik)),  &
                          vpploc(i),dvpploc(i),grid_function(i)
  end do
  close(4)

! multiply sqrt((2l+1)/4pi)/r**(l+1) for radial w.f.
  do l=0,pp%mlps(ik)
    do i=1,pp%mr(ik)
      uvpp(i,l)=uvpp(i,l)*sqrt((2*l+1.d0)/(4*pi))/(pp%rad(i+1,ik))**(l+1)
      duvpp(i,l)=duvpp(i,l)*sqrt((2*l+1.d0)/(4*pi))/(pp%rad(i+1,ik))**(l+1) &
           &- (l+1.d0)*uvpp(i,l)/pp%rad(i+1,ik)
    enddo
    uvpp(0,l)=2.d0*uvpp(1,l)-uvpp(2,l)
    duvpp(0,l)=2.d0*duvpp(1,l)-duvpp(2,l)
  enddo

  do l=0,pp%mlps(ik)
    do i=1,pp%nrps(ik)
      pp%vloctbl(i,ik)=vpploc(i-1)
      pp%dvloctbl(i,ik)=dvpploc(i-1)
      pp%udvtbl(i,l,ik)=uvpp(i-1,l)
      pp%dudvtbl(i,l,ik)=duvpp(i-1,l)
    enddo
    if (pp%inorm(l,ik) == 0) cycle
    pp%udvtbl(1:pp%nrps(ik),l,ik)=pp%udvtbl(1:pp%nrps(ik),l,ik)/pp%anorm(l,ik)
    pp%dudvtbl(1:pp%nrps(ik),l,ik)=pp%dudvtbl(1:pp%nrps(ik),l,ik)/pp%anorm(l,ik)
  enddo


  pp%flag_nlcc = pp%flag_nlcc.or.flag_nlcc_element(ik)
  pp%rho_nlcc_tbl(:,ik)=0d0; pp%tau_nlcc_tbl(:,ik)=0d0
  if(.not.flag_nlcc_element(ik))return
  do i=1,pp%nrps(ik)
    if(rhor_nlcc(i-1,0)/rhor_nlcc(0,0) < 1d-7)exit
    pp%rho_nlcc_tbl(i,ik)=rhor_nlcc(i-1,0)
    pp%tau_nlcc_tbl(i,ik)=0.25d0*rhor_nlcc(i-1,1)**2/rhor_nlcc(i-1,0)
  end do

  return
end subroutine making_ps_with_masking
!====
subroutine making_ps_without_masking(pp,ik,flag_nlcc_element,rhor_nlcc)
  use salmon_pp,only : pp_info
  use salmon_global, only: nelem
  implicit none
  real(8),parameter :: pi=3.141592653589793d0 ! copied from salmon_math
  type(pp_info),intent(inout) :: pp
  integer,intent(in) :: ik
  logical,intent(in) :: flag_nlcc_element(nelem)
  real(8),intent(in) :: rhor_nlcc(0:pp%nrmax0,0:2)
  integer :: i,l
  real(8) :: r1,r2,r3,r4

! multiply sqrt((2l+1)/4pi)/r**(l+1) for radial w.f.
  do l=0,pp%mlps(ik)
    do i=1,pp%mr(ik)
      pp%upp(i,l)=pp%upp(i,l)*sqrt((2*l+1.d0)/(4*pi))/(pp%rad(i+1,ik))**(l+1)
    enddo
    pp%upp(0,l)=pp%upp(1,l)
!    pp%upp(0,l)=2.d0*pp%upp(1,l)-pp%upp(2,l)
  enddo

  do l=0,pp%mlps(ik)
    do i=1,pp%mr(ik)-1
      r1 = pp%rad(i+1,ik)-pp%rad(i,ik)
      r2 = pp%rad(i+1,ik)-pp%rad(i+2,ik)
      r3 = pp%rad(i+2,ik)-pp%rad(i,ik)
      r4 = r1/r2
      pp%dvpp(i,l)=(r4+1.d0)*(pp%vpp(i,l)-pp%vpp(i-1,l))/r1-(pp%vpp(i+1,l)-pp%vpp(i-1,l))/r3*r4
      pp%dupp(i,l)=(r4+1.d0)*(pp%upp(i,l)-pp%upp(i-1,l))/r1-(pp%upp(i+1,l)-pp%upp(i-1,l))/r3*r4
    end do
    pp%dvpp(0,l)=pp%dvpp(1,l)
    pp%dvpp(pp%mr(ik),l)=pp%dvpp(pp%mr(ik)-1,l)
    pp%dupp(0,l)=pp%dupp(1,l)
    pp%dupp(pp%mr(ik),l)=pp%dupp(pp%mr(ik)-1,l)
!    pp%dvpp(0,l)=2.d0*pp%dvpp(1,l)-pp%dvpp(2,l)
!    pp%dvpp(pp%mr(ik),l)=2.d0*pp%dvpp(pp%mr(ik)-1,l)-pp%dvpp(pp%mr(ik)-2,l)
!    pp%dupp(0,l)=2.d0*pp%dupp(1,l)-pp%dupp(2,l)
!    pp%dupp(pp%mr(ik),l)=2.d0*pp%dupp(pp%mr(ik)-1,l)-pp%dupp(pp%mr(ik)-2,l)
  end do

  do l=0,pp%mlps(ik)
    do i=1,pp%nrps(ik)
      pp%vloctbl(i,ik)=pp%vpp(i-1,pp%lref(ik))
      pp%dvloctbl(i,ik)=pp%dvpp(i-1,pp%lref(ik))
      pp%udvtbl(i,l,ik)=(pp%vpp(i-1,l)-pp%vpp(i-1,pp%lref(ik)))*pp%upp(i-1,l)
      pp%dudvtbl(i,l,ik)=(pp%dvpp(i-1,l)-pp%dvpp(i-1,pp%lref(ik)))*pp%upp(i-1,l)   &
                          + (pp%vpp(i-1,l)-pp%vpp(i-1,pp%lref(ik)))*pp%dupp(i-1,l)
    enddo
    if (pp%inorm(l,ik) == 0) cycle
    pp%udvtbl(1:pp%nrps(ik),l,ik)=pp%udvtbl(1:pp%nrps(ik),l,ik)/pp%anorm(l,ik)
    pp%dudvtbl(1:pp%nrps(ik),l,ik)=pp%dudvtbl(1:pp%nrps(ik),l,ik)/pp%anorm(l,ik)
  enddo

  pp%flag_nlcc = pp%flag_nlcc.or.flag_nlcc_element(ik)
  pp%rho_nlcc_tbl(:,ik)=0d0; pp%tau_nlcc_tbl(:,ik)=0d0
  if(.not.flag_nlcc_element(ik))return
  do i=1,pp%nrps(ik)
    if(rhor_nlcc(i-1,0)/rhor_nlcc(0,0) < 1d-7)exit
    pp%rho_nlcc_tbl(i,ik)=rhor_nlcc(i-1,0)
    pp%tau_nlcc_tbl(i,ik)=0.25d0*rhor_nlcc(i-1,1)**2/rhor_nlcc(i-1,0)
  end do

  return
end subroutine making_ps_without_masking
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine ps_masking(pp,uvpp,duvpp,ik,hx,hy,hz)
  use salmon_pp,only : pp_info
  use salmon_global,only :ps_format,alpha_mask,gamma_mask
  use salmon_math, only: xjl, dxjl
  implicit none
  type(pp_info),intent(inout) :: pp
  real(8),parameter :: pi=3.141592653589793d0 ! copied from salmon_math
!argument
  integer,intent(in) :: ik
  real(8),intent(inout) :: uvpp(0:pp%nrmax0,0:pp%lmax0)
  real(8),intent(out) :: duvpp(0:pp%nrmax0,0:pp%lmax0)
  real(8),intent(in) :: hx,hy,hz
!local variable
!Normalized mask function
  integer,parameter :: nkmax=1000
  integer :: i,j,l
  real(8) :: kmax,k1,k2,kr1,kr2,dr,dk
  real(8),allocatable :: radk(:),wk(:,:) !Fourier staffs
!Mask function
  real(8),allocatable :: mask(:),dmask(:)
!Functions

!Reconstruct radial coordinate Rps(ik) and NRps(ik)
  pp%rps(ik) = gamma_mask*pp%rps(ik)
  do i=1,pp%nrmax0
    if (pp%rad(i,ik) > pp%rps(ik)) exit
  end do
  pp%nrps(ik)=i
  pp%rps(ik) = pp%rad(pp%nrps(ik),ik)
  allocate(mask(pp%nrps(ik)),dmask(pp%nrps(ik)))

  call make_mask_function(pp,mask,dmask,ik)

!Make
  do i = 0,pp%nrps(ik)-1
    do l = 0,pp%mlps(ik)
      uvpp(i,l) = uvpp(i,l)/mask(i+1)
    end do
  end do

  allocate(radk(nkmax),wk(nkmax,0:pp%mlps(ik)))
  wk(:,:)=0.d0
!  kmax = alpha_mask*Pi*sqrt(1.d0/Hx**2+1.d0/Hy**2+1.d0/Hz**2)
  kmax = alpha_mask*pi/max(Hx,Hy,Hz)
  do i = 1,nkmax
    radk(i) = kmax*(dble(i-1)/dble(nkmax-1))
  end do

  do i=1,nkmax
    do j=1,pp%mr(ik)-1
      kr1 = radk(i)*pp%rad(j,ik)
      kr2 = radk(i)*pp%rad(j+1,ik)
      dr = pp%rad(j+1,ik) - pp%rad(j,ik)
      do l=0,pp%mlps(ik)
        wk(i,l) = wk(i,l) &
             &+ 0.5d0*(xjl(kr1,l)*uvpp(j-1,l) + xjl(kr2,l)*uvpp(j,l))*dr
      end do
    end do
  end do

  open(4,file="PSFourier_"//trim(pp%atom_symbol(ik))//"_"//trim(ps_format(ik))//".dat")
  write(4,*) "# Kmax, NKmax =",kmax,nkmax
  write(4,*) "# Mlps(ik), Lref(ik) =",pp%mlps(ik), pp%lref(ik)
  write(4,*) "#  Pi/max(Hx,Hy,Hz) =", pi/max(Hx,Hy,Hz)
  write(4,*) "#  Pi*sqrt(1.d0/Hx**2+1.d0/Hy**2+1.d0/Hz**2) =", pi*sqrt(1.d0/Hx**2+1.d0/Hy**2+1.d0/Hz**2)
  do i=1,nkmax
    if(radk(i) < (Pi/max(Hx,Hy,Hz))) then
      write(4,'(8e21.12)') radk(i),(wk(i,l),l=0,pp%mlps(ik)),1.d0
    else 
      write(4,'(8e21.12)') radk(i),(wk(i,l),l=0,pp%mlps(ik)),0.d0
    end if
  end do
  close(4)

  uvpp = 0.d0; duvpp=0.d0
  do i=1,nkmax-1
    do j=1,pp%mr(ik)
      kr1 = radk(i)*pp%rad(j,ik)
      kr2 = radk(i+1)*pp%rad(j,ik)
      k1 = radk(i)
      k2 = radk(i+1)
      dk = radk(i+1) - radk(i)
      do l=0,pp%mlps(ik)
        uvpp(j-1,l) = uvpp(j-1,l) &
             &+ 0.5d0*(xjl(kr1,l)*wk(i,l) + xjl(kr2,l)*wk(i+1,l))*dk
        duvpp(j-1,l) = duvpp(j-1,l) &
             &+ 0.5d0*(k1*dxjl(kr1,l)*wk(i,l) + k2*dxjl(kr2,l)*wk(i+1,l))*dk
      end do
    end do
  end do

  do l=0,pp%mlps(ik)
    uvpp(pp%mr(ik),l) = 2.d0*uvpp(pp%mr(ik)-1,l) - uvpp(pp%mr(ik)-2,l)
    duvpp(pp%mr(ik),l) = 2.d0*duvpp(pp%mr(ik)-1,l) - duvpp(pp%mr(ik)-2,l)
  end do
  uvpp = (2.d0/pi)*uvpp
  duvpp = (2.d0/pi)*duvpp

  do i=0,pp%nrps(ik)-1
    do l = 0,pp%mlps(ik)
!Derivative calculation before constructing uvpp to avoid overwrite
      duvpp(i,l) = duvpp(i,l)*mask(i+1) + uvpp(i,l)*dmask(i+1)
      uvpp(i,l) = uvpp(i,l)*mask(i+1)
    end do
  end do

  deallocate(radk,wk)

  return

  contains

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

  subroutine make_mask_function(pp,rmask,dmask,ik)
!Subroutine Make_mask_function
!Name of variables are taken from ***
    use salmon_pp,only : pp_info
    use salmon_global, only : eta_mask
    implicit none
    real(8),parameter :: pi=3.141592653589793d0 ! copied from salmon_math
    type(pp_info),intent(inout) :: pp
!Arguments
    integer,intent(in) :: ik
    real(8),intent(inout) :: rmask(pp%nrps(ik)),dmask(pp%nrps(ik))
!local variables
    integer,parameter :: M = 200
!  real(8),parameter :: eta = 15.d0
    integer :: i,j, i3,i2,i1
    real(8) :: xp,xm,dx,nmask0,kx1,dk
    real(8) :: x(M),nmask(M),mat(M,M),k(M),nmask_k(M)
!Lapack dsyev
    integer :: INFO,LWORK
    real(8),allocatable :: WORK(:),W(:)

!Making normalized mask function in radial coordinate
    do i = 1,M
      x(i) = dble(i)/dble(M)
    end do
    do i = 1,M
      xp = 2.d0*x(i)
      mat(i,i) = sin(xp*eta_mask)/xp + dble(M)*pi - eta_mask
      do j = i+1,M
        xp = x(i) + x(j)
        xm = x(i) - x(j)
        mat(i,j) = sin(xp*eta_mask)/xp - sin(xm*eta_mask)/xm
        mat(j,i) = mat(i,j)
      end do
    end do
  
    allocate(W(M))
    LWORK = max(1,3*M - 1)
    allocate(WORK(LWORK))
    call dsyev('V','U',M,mat,M,W,WORK,LWORK,INFO)
    deallocate(WORK,W)
    nmask0 = 3.d0*mat(1,1)/x(1) - 3.d0*mat(2,1)/x(2) + mat(3,1)/x(3)
    do i = 1,M
      nmask(i) = mat(i,1)/x(i)/nmask0
    end do
    nmask0 = nmask0/nmask0
  
    open(4,file="nmask.dat")
    write(4,*) "# M =",M
    write(4,*) 0,nmask0
    do i= 1,M
      write(4,*) x(i),nmask(i)
    end do
    close(4)
  
!Taking Fourier transformation
    do i = 1,M
      k(i)=pi*dble(i)
    end do
    dx = x(2)-x(1)
    nmask_k(:) = 0.d0
    do i = 1,M
      do j = 1,M
        kx1 = k(i)*x(j)
        nmask_k(i) = nmask_k(i) + nmask(j)*kx1*sin(kx1) 
      end do
      nmask_k(i) = nmask_k(i)*dx/k(i)**2
    end do
  
    open(4,file="nmask_k.dat")
    write(4,*) 0,  3.d0*nmask_k(1) - 3.d0*nmask_k(2) + nmask_k(3)
    do i= 1,M
      write(4,*) k(i),nmask_k(i)
    end do
    close(4)
  
!  allocate(mask(M),dmask(M))!debug
!Making normalized mask function in radial coordinate
    rmask(:) = 0.d0; dmask(:)=0.d0
    dk = k(2) - k(1)
    do i=2,pp%nrps(ik) !Avoiding divide by zero
      do j = 1,M
        kx1 = k(j)*pp%rad(i,ik)/pp%rps(ik)
        rmask(i) =  rmask(i) + nmask_k(j)*kx1*sin(kx1)
        dmask(i)= dmask(i) + nmask_k(j)*(kx1**2*cos(kx1)-kx1*sin(kx1))
      end do
      rmask(i) = (2.d0/pi)* rmask(i)*dk*pp%rps(ik)**2/pp%rad(i,ik)**2
      dmask(i)= (2.d0/pi)*dmask(i)*dk*pp%rps(ik)**2/pp%rad(i,ik)**3 
    end do
    rmask(1) =  rmask(2)-( rmask(3)- rmask(2))/(pp%rad(3,ik)-pp%rad(2,ik))*(pp%rad(2,ik)-pp%rad(1,ik))
    dmask(1)= dmask(2)-(dmask(3)-dmask(2))/(pp%rad(3,ik)-pp%rad(2,ik))*(pp%rad(2,ik)-pp%rad(1,ik))
    i1=pp%nrps(ik)-2
    i2=pp%nrps(ik)-1
    i3=pp%nrps(ik)
     rmask(i3)=  rmask(i2)+( rmask(i2)- rmask(i1))/(pp%rad(i2,ik)-pp%rad(i1,ik))*(pp%rad(i3,ik)-pp%rad(i2,ik))
    dmask(i3)= dmask(i2)+(dmask(i2)-dmask(i1))/(pp%rad(i2,ik)-pp%rad(i1,ik))*(pp%rad(i3,ik)-pp%rad(i2,ik))
  
    open(4,file="mask.dat")
    write(4,*) "# rps(ik), nrps(ik) =",pp%rps(ik), pp%nrps(ik)
    do i= 1,pp%nrps(ik)
      write(4,'(8e22.10)') pp%rad(i,ik),rmask(i),dmask(i)
    end do
    close(4)
  
    return
  end subroutine make_mask_function
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
end subroutine ps_masking
