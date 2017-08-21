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
module read_pslfile_sub

use scf_data

integer,parameter :: Nlps=3, Nlmps=16

integer,allocatable :: Mlps0(:)

integer :: maxMr
real(8) :: rmin_step
real(8) :: rmaxRps

real(8),allocatable :: Zps(:)              ! Pseudo charge
real(8),allocatable :: Rps(:)              ! Core radius
real(8),allocatable :: Mass(:)             ! Atomic weight

real(8) :: rPC,rRC(0:Nlps)

integer,allocatable :: Mr(:)
real(8),allocatable :: step(:)

real(8), allocatable :: upp_f(:,:,:)
real(8), allocatable :: rhopp_f(:,:)
real(8), allocatable :: vpp_f(:,:,:)

real(8), allocatable :: rad_f(:,:)

contains
!==================================================================================================
subroutine read_pslfile
use salmon_parallel, only: nproc_id_global
use salmon_communication, only: comm_is_root
implicit none
integer :: ak
character(256) :: ps_file(MKI)
character(256) :: ps_file_tmp
integer :: ips_type,nlen_psf

allocate( Mlps0(MKI) )
allocate( Mr(MKI) )
allocate( step(MKI) )

allocate( Zps(MKI) )
allocate( Rps(MKI) )
allocate( Mass(MKI) )

maxMr=0
rmin_step=1.d8
rmaxRps=0

do ak=1,MKI

  ps_file_tmp = trim(pseudo_file(ak))
  nlen_psf = len_trim(ps_file_tmp)

  if(ps_file_tmp(nlen_psf+1-8:nlen_psf) == '_rps.dat')then
     ips_type = n_Yabana_Bertsch_psformat
     ps_format(ak) = 'KY'
  else if(ps_file_tmp(nlen_psf+1-6:nlen_psf) == '.pspnc')then
     ips_type = n_ABINIT_psformat
     ps_format(ak) = 'ABINIT'
  else if(ps_file_tmp(nlen_psf+1-4:nlen_psf) == '.cpi')then
     ips_type = n_FHI_psformat
     ps_format(ak) = 'FHI'
  else
     stop 'Unprepared ps_format is required input_pseudopotential_YS'
  end if

  ipsfileform(ak) = ips_type

  if(ipsfileform(ak)==3)then
    Mlps0(ak)=Mlps(ak)
  end if

  if(comm_is_root(nproc_id_global))then
    write(*,*) iZatom(ak)
  end if

  select case( iZatom(ak) )
    case(1) ; Atomname(ak)='H'  ; Mass(ak) = 1.008d0
    case(2) ; Atomname(ak)='He' ; Mass(ak) = 4.003d0
    case(3) ; Atomname(ak)='Li' ; Mass(ak) = 6.941d0
    case(4) ; Atomname(ak)='Be' ; Mass(ak) = 9.012d0
    case(5) ; Atomname(ak)='B'  ; Mass(ak) = 10.811d0
    case(6) ; Atomname(ak)='C'  ; Mass(ak) = 12.011d0
    case(7) ; Atomname(ak)='N'  ; Mass(ak) = 14.007d0
    case(8) ; Atomname(ak)='O'  ; Mass(ak) = 15.999d0
    case(9) ; Atomname(ak)='F'  ; Mass(ak) = 18.998d0
    case(10); Atomname(ak)='Ne' ; Mass(ak) = 20.180d0
    case(11); Atomname(ak)='Na' ; Mass(ak) = 22.990d0
    case(12); Atomname(ak)='Mg' ; Mass(ak) = 24.305d0
    case(13); Atomname(ak)='Al' ; Mass(ak) = 26.982d0
    case(14); Atomname(ak)='Si' ; Mass(ak) = 28.086d0
    case(15); Atomname(ak)='P'  ; Mass(ak) = 30.974d0
    case(16); Atomname(ak)='S'  ; Mass(ak) = 32.065d0
    case(17); Atomname(ak)='Cl' ; Mass(ak) = 35.453d0
    case(18); Atomname(ak)='Ar' ; Mass(ak) = 39.948d0
    case(19); Atomname(ak)='K'  ; Mass(ak) = 39.098d0
    case(20); Atomname(ak)='Ca' ; Mass(ak) = 40.078d0
    case(21); Atomname(ak)='Sc' ; Mass(ak) = 44.956d0
    case(22); Atomname(ak)='Ti' ; Mass(ak) = 47.867d0
    case(23); Atomname(ak)='V'  ; Mass(ak) = 50.942d0
    case(24); Atomname(ak)='Cr' ; Mass(ak) = 51.996d0
    case(25); Atomname(ak)='Mn' ; Mass(ak) = 54.938d0
    case(26); Atomname(ak)='Fe' ; Mass(ak) = 55.845d0
    case(27); Atomname(ak)='Co' ; Mass(ak) = 58.933d0
    case(28); Atomname(ak)='Ni' ; Mass(ak) = 58.693d0
    case(29); Atomname(ak)='Cu' ; Mass(ak) = 63.546d0
    case(30); Atomname(ak)='Zn' ; Mass(ak) = 65.38d0
    case(31); Atomname(ak)='Ga' ; Mass(ak) = 69.723d0
    case(32); Atomname(ak)='Ge' ; Mass(ak) = 72.64d0
    case(33); Atomname(ak)='As' ; Mass(ak) = 74.922d0
    case(34); Atomname(ak)='Se' ; Mass(ak) = 78.96d0
    case(35); Atomname(ak)='Br' ; Mass(ak) = 79.904d0
    case(36); Atomname(ak)='Kr' ; Mass(ak) = 83.798d0
    case(47); Atomname(ak)='Ag' ; Mass(ak) = 107.868d0
    case(79); Atomname(ak)='Au' ; Mass(ak) = 196.967d0
  end select


  select case( ipsfileform(ak) )
    case(1)  
      select case( iZatom(ak) )
        case default ; ps_file(ak)=trim(ps_file_tmp) !trim(Atomname(ak))//'_rps.dat'
        case(10,18,36,47,79) ; ps_file(ak)= trim(ps_file_tmp)!trim(Atomname(ak))//'_crps.dat'
      end select
      call read_Mr_YB(ak,ps_file)
    case(2)
      ps_file(ak)=trim(ps_file_tmp) !trim(Atomname(ak))//'.pspnc'
      call read_Mr_ABINIT(ak,ps_file)
    case(3)
      ps_file(ak)=trim(ps_file_tmp) !trim(Atomname(ak))//'.cpi'
      call read_Mr_fhi(ak,ps_file)
  end select
end do

maxMr=maxval(Mr)

allocate(upp_f(0:maxMr,0:Nlps,MKI))
allocate(rhopp_f(0:maxMr,MKI))
allocate(vpp_f(0:maxMr,0:Nlps,MKI))
allocate(rad_f(0:maxMr,MKI) )

do ak=1,MKI
  select case( ipsfileform(ak) )
    case(1) 
      call read_psl_YB(ak,ps_file)
    case(2) 
      call read_psl_ABINIT(ak,ps_file)
    case(3) 
      call read_psl_fhi(ak,ps_file)
      call setRps_fhi(ak)
  end select
end do

maxlm=0
do ak=1,MKI
  if(Mlps(ak)>maxlm) maxlm=Mlps(ak)
end do
maxlm=(maxlm+1)**2

end subroutine read_pslfile

!======================================================================
subroutine read_Mr_YB(ak,ps_file)
implicit none
integer :: ak
character(256) :: ps_file(MKI)

open(4,file=ps_file(ak),status='old')
read(4,*) Mr(ak)
close(4)
return

end subroutine read_Mr_YB
!======================================================================
subroutine read_Mr_ABINIT(ak,ps_file)
implicit none
integer :: ak
character(256) :: ps_file(MKI)
real(8) :: zatom, zion, pspdat,pspcod,pspxc,lmaxabinit,lloc,mmax,r2well
character(1) :: dummy_text

open(4,file=ps_file(ak),status='old')
read(4,*) dummy_text
read(4,*) zatom, zion, pspdat
read(4,*) pspcod,pspxc,lmaxabinit,lloc,mmax,r2well
close(4)

Mr(ak)=mmax

end subroutine read_Mr_ABINIT

!======================================================================
subroutine read_Mr_fhi(ak,ps_file)
implicit none
integer :: ak
integer :: i
character(256) :: ps_file(MKI)
character(1) :: dummy_text

open(4,file=ps_file(ak),status='old')
do i=1,11
  read(4,*) dummy_text
end do
read(4,*) Mr(ak)
close(4)

end subroutine read_Mr_fhi

!======================================================================
!======================================================================
subroutine read_psl_YB(ak,ps_file)
implicit none
integer :: ak
character(256) :: ps_file(MKI)
integer :: L
integer :: i
real(8) :: r

open(4,file=ps_file(ak),status='old')

read(4,*) Mr(ak),step(ak),Mlps0(ak),Zps(ak)

step(ak)=step(ak)/a_B

if( Mr(ak) > maxMr ) maxMr=Mr(ak)
if( step(ak) < rmin_step ) rmin_step=step(ak)

if ( Mlps0(ak) < Mlps(ak) ) then
  write(*,*) "too large Mlps value for ",ak," atom"
  stop
end if

rRC=0.d0

read(4,*) rPC, ( rRC(L) , L=0,Mlps0(ak) )
do L=0,Mlps0(ak)
  rRC(L)=rRC(L)/a_B
end do

Rps(ak)=maxval(rRC)
if(maxval( rRC ) > rmaxRps) rmaxRps=maxval( rRC )

do i=0,Mr(ak)
  read(4,*) r,rhopp_f(i,ak),(vpp_f(i,L,ak),L=0,Mlps0(ak))
end do
do i=0,Mr(ak)
  read(4,*) r,(upp_f(i,L,ak),L=0,Mlps0(ak))
end do

do i=0,Mr(ak)
  rhopp_f(i,ak)=rhopp_f(i,ak)*(a_B)**3d0
  do L=0,Mlps0(ak)
    vpp_f(i,L,ak)=vpp_f(i,L,ak)/2d0/Ry
  end do
end do

close(4)

return

end subroutine read_psl_YB

!==================================================================================================
!======================================================================
subroutine read_psl_ABINIT(ak,ps_file)
!See http://www.abinit.org/downloads/psp-links/psp-links/lda_tm
  implicit none
  integer,intent(in) :: ak
  character(256),intent(in) :: ps_file(MKI)
  integer :: i
  real(8) :: rZps
  integer :: ll
  real(8) :: zatom, zion, pspdat,pspcod,pspxc,lmaxabinit,lloc,mmax,r2well,dl
  real(8) :: e99_0,e99_9,nproj,rcpsp,rms,ekb1,ekb2,epsatm,rchrg,fchrg,qchrg
  character(1) :: dummy_text

  step(ak)=0.01d0/a_B

  open(4,file=ps_file(ak),status='old')
  read(4,*) dummy_text
  read(4,*) zatom, zion, pspdat
  rZps = zion
  Zps(ak)=int(rZps+1d-10)
  read(4,*) pspcod,pspxc,lmaxabinit,lloc,mmax,r2well
  Mlps0(ak)=lmaxabinit
  Mr(ak) = mmax - 1

  rRC=0.d0
  do ll=0,Mlps0(ak)
    read(4,*) dl,e99_0,e99_9,nproj,rcpsp
    read(4,*) rms,ekb1,ekb2,epsatm
    rRC(ll) = rcpsp
  end do

  Rps(ak)=maxval(rRC)
  if(maxval( rRC ) > rmaxRps) rmaxRps=maxval( rRC )

  read(4,*) rchrg,fchrg,qchrg
  do ll=0,Mlps0(ak)
    read(4,*) dummy_text
    do i=1,(Mr(ak)+1)/3
      read(4,*) vpp_f(3*(i-1),ll,ak),vpp_f(3*(i-1)+1,ll,ak),vpp_f(3*(i-1)+2,ll,ak)
    end do
  end do
  do ll=0,Mlps0(ak)
    read(4,*) dummy_text
    do i=1,(Mr(ak)+1)/3
      read(4,*) upp_f(3*(i-1),ll,ak),upp_f(3*(i-1)+1,ll,ak),upp_f(3*(i-1)+2,ll,ak)
    end do
  end do
  close(4)

  do i=1,Mr(ak)
    rhopp_f(i,ak)=0.d0
  end do


return

End Subroutine read_psl_ABINIT

!==================================================================================================
!======================================================================
subroutine read_psl_ABINIT_PBE(ak,ps_file)
!See http://www.abinit.org/downloads/psp-links/psp-links/gga-phi
implicit none
integer,intent(in) :: ak
character(256),intent(in) :: ps_file(MKI)
integer :: i,ibox
real(8) :: rZps
integer :: ll
real(8) :: zatom, zion, pspdat,pspcod,pspxc,lmaxabinit,lloc,mmax,r2well
real(8) :: rchrg,fchrg,qchrg
character(1) :: dummy_text
real(8) :: step_tmp(0:3)

step_tmp(0:3)=1.d8

open(4,file=ps_file(ak),status='old')
read(4,*) dummy_text
read(4,*) zatom, zion, pspdat
rZps = zion
Zps(ak)=int(rZps+1d-10)
read(4,*) pspcod,pspxc,lmaxabinit,lloc,mmax,r2well
Mlps0(ak)=lmaxabinit
Mr(ak) = mmax - 1
read(4,*) rchrg,fchrg,qchrg
do i=5,18
  read(4,*)
end do

do ll=0,Mlps0(ak)
  read(4,*) ibox,step_tmp(ll)
  do i=0,Mr(ak)
    read(4,*) ibox,rad_f(i,ak),upp_f(i,ll,ak),vpp_f(i,ll,ak)
  end do
  rad_f(0,ak)=0.d0
end do
close(4)

step(ak) = minval(step_tmp(0:Mlps0(ak)))

do i=1,Mr(ak)
  rhopp_f(i,ak)=0.d0
end do

return

End Subroutine read_psl_ABINIT_PBE

!==================================================================================================
!======================================================================
subroutine read_psl_fhi(ak,ps_file)
implicit none
integer,intent(in) :: ak
character(256),intent(in) :: ps_file(MKI)
integer :: i,ibox
integer :: ll
character(1) :: dummy_text
real(8) :: step_tmp(0:3)

step_tmp(0:3)=1.d8

open(4,file=ps_file(ak),status='old')

call set_Zps(ak)

do i=1,11
  read(4,*) dummy_text
end do

do ll=0,Mlps0(ak)
  read(4,*) ibox,step_tmp(ll)
  Mr(ak)=ibox-1
  do i=0,Mr(ak)
    read(4,*) ibox,rad_f(i,ak),upp_f(i,ll,ak),vpp_f(i,ll,ak)
  end do
  rad_f(0,ak)=0.d0
end do
close(4)

step(ak) = minval(step_tmp(0:Mlps0(ak)))

do i=1,Mr(ak)
  rhopp_f(i,ak)=0.d0
end do

return

End Subroutine read_psl_fhi

!==================================================================================================
!======================================================================
Subroutine setRps_ABINIT_PBE(ak)
implicit none
integer,intent(in) :: ak
real(8) :: rcore(0:4,1:120)
integer :: ll

!rcore list in unit of bohr

rcore(0,1)=1.2763463d0
rcore(1,1)=1.2763463d0
rcore(2,1)=0.3502204d0
rcore(3,1)=1.2763463d0

rcore(0:3,2)=1.0916080d0

rcore(0,3)=2.1818685d0
rcore(1,3)=2.1818685d0
rcore(2,3)=2.4649735d0
rcore(3,3)=2.4649735d0

rcore(0,4)=2.0886093d0
rcore(1,4)=2.0886093d0
rcore(2,4)=2.3596136d0
rcore(3,4)=2.3596136d0

rcore(0:3,5)=1.6708874d0

rcore(0:3,6)=1.4981530d0

rcore(0:3,7)=1.4157818d0

rcore(0:3,8)=1.3995488d0

rcore(0:3,9)=1.3385227d0

rcore(0:3,10)=1.2961595d0

rcore(0,11)=2.7011879d0
rcore(1,11)=2.7011879d0
rcore(2,11)=2.9781169d0
rcore(3,11)=2.9781169d0

rcore(0,12)=1.7907575d0
rcore(1,12)=1.9743480d0
rcore(2,12)=2.1242906d0
rcore(3,12)=2.1242906d0

rcore(0,13)=2.0873202d0
rcore(1:3,13)=2.4760889d0

rcore(0,14)=1.7039185d0
rcore(1,14)=1.8786063d0
rcore(2,14)=2.0212776d0
rcore(3,14)=2.0212776d0

rcore(0,15)=1.6296049d0
rcore(1,15)=1.7966740d0
rcore(2,15)=1.9331230d0
rcore(3,15)=1.9331230d0

rcore(0,16)=1.5654902d0
rcore(1,16)=1.6843819d0
rcore(2,16)=1.8570667d0
rcore(3,16)=1.8570667d0

rcore(0,17)=1.5853006d0
rcore(1,17)=1.5853006d0
rcore(2,17)=1.7909988d0
rcore(3,17)=1.7909988d0

rcore(0,18)=1.4259180d0
rcore(1,18)=1.4259180d0
rcore(2,18)=1.7332789d0
rcore(3,18)=1.7332789d0

rcore(0,19)=2.8781425d0
rcore(1,19)=2.6749895d0
rcore(2,19)=3.1732130d0
rcore(3,19)=3.1732130d0

rcore(0,20)=2.8709747d0
rcore(1,20)=2.8709747d0
rcore(2,20)=2.5412400d0
rcore(3,20)=2.5412400d0

rcore(0,21)=2.7342616d0
rcore(1,21)=2.7342616d0
rcore(2,21)=2.4202286d0
rcore(3,21)=2.4202286d0

rcore(0,22)=2.6099770d0
rcore(1,22)=2.6099770d0
rcore(2,22)=2.3672806d0
rcore(3,22)=2.3672806d0

rcore(0,23)=2.5581633d0
rcore(1,23)=2.5581633d0
rcore(2,23)=2.3775959d0
rcore(3,23)=2.3775959d0

rcore(0,24)=2.5121270d0
rcore(1,24)=2.5121270d0
rcore(2,24)=2.3348091d0
rcore(3,24)=2.3348091d0

rcore(0,25)=2.4116419d0
rcore(1,25)=2.4712095d0
rcore(2,25)=2.2967797d0
rcore(3,25)=2.2967797d0

rcore(0,26)=2.3188865d0
rcore(1,26)=2.4348542d0
rcore(2,26)=2.2084421d0
rcore(3,26)=2.2084421d0

rcore(0,27)=2.2881569d0
rcore(1,27)=2.3446744d0
rcore(2,27)=2.1791761d0
rcore(3,27)=2.1791761d0

rcore(0,28)=2.6173915d0
rcore(1,28)=2.3167811d0
rcore(2,28)=2.1013484d0
rcore(3,28)=2.1013484d0

rcore(0,29)=2.6535190d0
rcore(1,29)=2.2921434d0
rcore(2,29)=2.0790016d0
rcore(3,29)=2.0790016d0

rcore(0,30)=2.3265479d0
rcore(1,30)=2.2704673d0
rcore(2,30)=2.0097016d0
rcore(3,30)=2.6933476d0

rcore(0,31)=2.0925765d0
rcore(1,31)=2.2514979d0
rcore(2,31)=2.0925765d0
rcore(3,31)=2.4823242d0

rcore(0,32)=1.9783190d0
rcore(1,32)=2.7167764d0
rcore(2,32)=2.4641489d0
rcore(3,32)=2.4641489d0

rcore(0,33)=1.9657537d0
rcore(1,33)=2.1672851d0
rcore(2,33)=2.1150435d0
rcore(3,33)=2.4484978d0

rcore(0,34)=1.8619473d0
rcore(1,34)=2.0528364d0
rcore(2,34)=1.8619473d0
rcore(3,34)=2.3764832d0

rcore(0,35)=1.7651496d0
rcore(1,35)=2.0939132d0
rcore(2,35)=1.6405570d0
rcore(3,35)=2.1986300d0

rcore(0,36)=1.7161177d0
rcore(1,36)=1.9866780d0
rcore(2,36)=1.4118000d0
rcore(3,36)=2.0860320d0

rcore(0,47)=2.4191908d0
rcore(1,47)=2.6029170d0
rcore(2,47)=2.4191908d0
rcore(3,47)=2.4191908d0

rcore(0,79)=2.5227012d0
rcore(1,79)=2.6488617d0
rcore(2,79)=2.5227012d0
rcore(3,79)=2.5227012d0

rRC=0.d0
do ll=0,Mlps0(ak)
  rRC(ll) = rcore(ll,iZatom(ak))
end do

Rps(ak)=maxval(rRC)
if(maxval( rRC ) > rmaxRps) rmaxRps=maxval( rRC )

End Subroutine setRps_ABINIT_PBE

!==================================================================================================
!======================================================================
Subroutine set_Zps(ak)
integer :: ak

select case( iZatom(ak) )
  case(1)   ; Zps(ak) = 1.d0
  case(2)   ; Zps(ak) = 2.d0
  case(3)   ; Zps(ak) = 1.d0
  case(4)   ; Zps(ak) = 2.d0
  case(5)   ; Zps(ak) = 3.d0
  case(6)   ; Zps(ak) = 4.d0
  case(7)   ; Zps(ak) = 5.d0
  case(8)   ; Zps(ak) = 6.d0
  case(9)   ; Zps(ak) = 7.d0
  case(10)  ; Zps(ak) = 8.d0
  case(11)  ; Zps(ak) = 1.d0
  case(12)  ; Zps(ak) = 2.d0
  case(13)  ; Zps(ak) = 3.d0
  case(14)  ; Zps(ak) = 4.d0
  case(15)  ; Zps(ak) = 5.d0
  case(16)  ; Zps(ak) = 6.d0
  case(17)  ; Zps(ak) = 7.d0
  case(18)  ; Zps(ak) = 8.d0
  case(19)  ; Zps(ak) = 1.d0
  case(20)  ; Zps(ak) = 2.d0
  case(21)  ; Zps(ak) = 3.d0
  case(22)  ; Zps(ak) = 4.d0
  case(23)  ; Zps(ak) = 5.d0
  case(24)  ; Zps(ak) = 6.d0
  case(25)  ; Zps(ak) = 7.d0
  case(26)  ; Zps(ak) = 8.d0
  case(27)  ; Zps(ak) = 9.d0
  case(28)  ; Zps(ak) = 10.d0
  case(29)  ; Zps(ak) = 11.d0
  case(30)  ; Zps(ak) = 12.d0
  case(31)  ; Zps(ak) = 3.d0
  case(32)  ; Zps(ak) = 4.d0
  case(33)  ; Zps(ak) = 5.d0
  case(34)  ; Zps(ak) = 6.d0
  case(35)  ; Zps(ak) = 7.d0
  case(36)  ; Zps(ak) = 8.d0
  case(37)  ; Zps(ak) = 1.d0
  case(38)  ; Zps(ak) = 2.d0
  case(39)  ; Zps(ak) = 3.d0
  case(40)  ; Zps(ak) = 4.d0
  case(41)  ; Zps(ak) = 5.d0
  case(42)  ; Zps(ak) = 6.d0
  case(43)  ; Zps(ak) = 7.d0
  case(44)  ; Zps(ak) = 8.d0
  case(45)  ; Zps(ak) = 9.d0
  case(46)  ; Zps(ak) = 10.d0
  case(47)  ; Zps(ak) = 11.d0
  case(48)  ; Zps(ak) = 12.d0
  case(49)  ; Zps(ak) = 3.d0
  case(50)  ; Zps(ak) = 4.d0
  case(51)  ; Zps(ak) = 5.d0
  case(52)  ; Zps(ak) = 6.d0
  case(53)  ; Zps(ak) = 7.d0
  case(54)  ; Zps(ak) = 8.d0
  case(55)  ; Zps(ak) = 1.d0
  case(56)  ; Zps(ak) = 2.d0
  case(57)  ; Zps(ak) = 3.d0
  case(58)  ; Zps(ak) = 4.d0
  case(59)  ; Zps(ak) = 5.d0
  case(60)  ; Zps(ak) = 6.d0
  case(61)  ; Zps(ak) = 7.d0
  case(62)  ; Zps(ak) = 8.d0
  case(63)  ; Zps(ak) = 9.d0
  case(64)  ; Zps(ak) = 10.d0
  case(65)  ; Zps(ak) = 11.d0
  case(66)  ; Zps(ak) = 12.d0
  case(67)  ; Zps(ak) = 13.d0
  case(68)  ; Zps(ak) = 14.d0
  case(69)  ; Zps(ak) = 15.d0
  case(70)  ; Zps(ak) = 16.d0
  case(71)  ; Zps(ak) = 17.d0
  case(72)  ; Zps(ak) = 4.d0
  case(73)  ; Zps(ak) = 5.d0
  case(74)  ; Zps(ak) = 6.d0
  case(75)  ; Zps(ak) = 7.d0
  case(76)  ; Zps(ak) = 8.d0
  case(77)  ; Zps(ak) = 9.d0
  case(78)  ; Zps(ak) = 10.d0
  case(79)  ; Zps(ak) = 11.d0
  case(80)  ; Zps(ak) = 12.d0
  case(81)  ; Zps(ak) = 3.d0
  case(82)  ; Zps(ak) = 4.d0
  case(83)  ; Zps(ak) = 5.d0
  case(84)  ; Zps(ak) = 6.d0
  case(85)  ; Zps(ak) = 7.d0
  case(86)  ; Zps(ak) = 8.d0
  case(87)  ; Zps(ak) = 1.d0
  case(88)  ; Zps(ak) = 2.d0
  case(89)  ; Zps(ak) = 3.d0
  case(90)  ; Zps(ak) = 4.d0
  case(91)  ; Zps(ak) = 5.d0
  case(92)  ; Zps(ak) = 6.d0
  case(93)  ; Zps(ak) = 7.d0
  case(94)  ; Zps(ak) = 8.d0
  case(95)  ; Zps(ak) = 9.d0
  case(96)  ; Zps(ak) = 10.d0
  case(97)  ; Zps(ak) = 11.d0
  case(98)  ; Zps(ak) = 12.d0
  case(99)  ; Zps(ak) = 13.d0
  case(100) ; Zps(ak) = 14.d0
  case(101) ; Zps(ak) = 15.d0
  case(102) ; Zps(ak) = 16.d0
  case(103) ; Zps(ak) = 17.d0
  case(104) ; Zps(ak) = 18.d0
  case(105) ; Zps(ak) = 19.d0
  case(106) ; Zps(ak) = 20.d0
  case(107) ; Zps(ak) = 21.d0
  case(108) ; Zps(ak) = 22.d0
  case(109) ; Zps(ak) = 23.d0
end select

End Subroutine set_Zps

!==================================================================================================
!======================================================================
Subroutine setRps_fhi(ak)
implicit none
integer,intent(in) :: ak
real(8) :: rcore(0:4,1:120)
integer :: ll

!rcore list in unit of bohr

rcore(0,1)=  1.2763463d0
rcore(1,1)=  1.2763463d0  
rcore(2,1)=  1.2763463d0  
rcore(3,1)=  1.2763463d0  

rcore(0,2)=  1.0916080d0  
rcore(1,2)=  1.0916080d0  
rcore(2,2)=  1.0916080d0  
rcore(3,2)=  1.0916080d0  

rcore(0,3)=  2.1818685d0  
rcore(1,3)=  2.1818685d0  
rcore(2,3)=  2.4649735d0  
rcore(3,3)=  2.4649735d0  

rcore(0,4)=  2.0886093d0  
rcore(1,4)=  2.0886093d0  
rcore(2,4)=  2.3596136d0  
rcore(3,4)=  2.3596136d0  

rcore(0,5)=  1.6708874d0  
rcore(1,5)=  1.6708874d0  
rcore(2,5)=  1.6708874d0  
rcore(3,5)=  1.6708874d0  

rcore(0,6)=  1.4981530d0  
rcore(1,6)=  1.4981530d0 
rcore(2,6)=  1.4981530d0 
rcore(3,6)=  1.4981530d0 

rcore(0,7)=  1.4157818d0  
rcore(1,7)=  1.4157818d0  
rcore(2,7)=  1.4157818d0  
rcore(3,7)=  1.4157818d0  

rcore(0,8)=  1.3995488d0  
rcore(1,8)=  1.3995488d0  
rcore(2,8)=  1.3995488d0  
rcore(3,8)=  1.3995488d0  

rcore(0,9)=  1.3385227d0  
rcore(1,9)=  1.3385227d0  
rcore(2,9)=  1.3385227d0  
rcore(3,9)=  1.3385227d0  

rcore(0,10)=  1.2961595d0  
rcore(1,10)=  1.2961595d0  
rcore(2,10)=  1.2961595d0  
rcore(3,10)=  1.2961595d0  

rcore(0,11)=  2.7011879d0  
rcore(1,11)=  2.7011879d0  
rcore(2,11)=  2.9781169d0  
rcore(3,11)=  2.9781169d0  

rcore(0,12)=  2.0873202d0  
rcore(1,12)=  2.4760889d0  
rcore(2,12)=  2.4760889d0  
rcore(3,12)=  2.4760889d0  

rcore(0,13)=  1.7907575d0
rcore(1,13)=  1.9743480d0
rcore(2,13)=  2.1242906d0
rcore(3,13)=  2.1242906d0

rcore(0,14)=  1.7039185d0
rcore(1,14)=  1.8786063d0
rcore(2,14)=  2.0212776d0
rcore(3,14)=  2.0212776d0

rcore(0,15)=  1.6296049d0
rcore(1,15)=  1.7966740d0
rcore(2,15)=  1.9331230d0
rcore(3,15)=  1.9331230d0

rcore(0,16)=  1.5654902d0
rcore(1,16)=  1.6843819d0
rcore(2,16)=  1.8570667d0
rcore(3,16)=  1.8570667d0

rcore(0,17)=  1.5853006d0
rcore(1,17)=  1.5853006d0
rcore(2,17)=  1.7909988d0
rcore(3,17)=  1.7909988d0

rcore(0,18)=  1.4259180d0
rcore(1,18)=  1.4259180d0
rcore(2,18)=  1.7332789d0
rcore(3,18)=  1.7332789d0

rcore(0,19)=  3.2515914d0
rcore(1,19)=  3.1732130d0
rcore(2,19)=  3.1732130d0
rcore(3,19)=  3.1732130d0

rcore(0,20)=  2.7342353d0
rcore(1,20)=  2.9418878d0
rcore(2,20)=  2.5412400d0
rcore(3,20)=  2.5412400d0

rcore(0,21)=  2.7342616d0
rcore(1,21)=  2.7342616d0
rcore(2,21)=  2.4202286d0
rcore(3,21)=  2.4202286d0

rcore(0,22)=  2.6099770d0
rcore(1,22)=  2.6099770d0
rcore(2,22)=  2.3672806d0
rcore(3,22)=  2.3672806d0

rcore(0,23)=  2.5581633d0
rcore(1,23)=  2.5581633d0
rcore(2,23)=  2.3775959d0
rcore(3,23)=  2.3775959d0

rcore(0,24)=  2.5121270d0
rcore(1,24)=  2.5121270d0
rcore(2,24)=  2.3348091d0
rcore(3,24)=  2.3348091d0

rcore(0,25)=  2.4116419d0
rcore(1,25)=  2.4712095d0
rcore(2,25)=  2.2967797d0
rcore(3,25)=  2.2967797d0

rcore(0,26)=  2.3188865d0
rcore(1,26)=  2.4348542d0
rcore(2,26)=  2.2084421d0
rcore(3,26)=  2.2084421d0

rcore(0,27)=  2.2881569d0
rcore(1,27)=  2.3446744d0
rcore(2,27)=  2.1791761d0
rcore(3,27)=  2.1791761d0

rcore(0,28)=  2.2064370d0
rcore(1,28)=  2.3167811d0
rcore(2,28)=  2.1013484d0
rcore(3,28)=  2.1013484d0

rcore(0,29)=  2.0790016d0
rcore(1,29)=  2.2921434d0
rcore(2,29)=  2.0790016d0
rcore(3,29)=  2.0790016d0

rcore(0,30)=  2.0097016d0
rcore(1,30)=  2.2704673d0
rcore(2,30)=  2.0097016d0
rcore(3,30)=  2.6933476d0

rcore(0,31)=  2.0925765d0
rcore(1,31)=  2.2514979d0
rcore(2,31)=  2.0925765d0
rcore(3,31)=  2.4823242d0

rcore(0,32)=  1.9783190d0
rcore(1,32)=  2.1811386d0
rcore(2,32)=  2.4641489d0
rcore(3,32)=  2.4641489d0

rcore(0,33)=  1.9657537d0
rcore(1,33)=  2.1672851d0
rcore(2,33)=  2.4484978d0
rcore(3,33)=  2.4484978d0

rcore(0,34)=  1.8619473d0
rcore(1,34)=  2.0528364d0
rcore(2,34)=  2.3764832d0
rcore(3,34)=  2.3764832d0

rcore(0,35)=  1.7651496d0
rcore(1,35)=  2.0939132d0
rcore(2,35)=  2.1986300d0
rcore(3,35)=  2.1986300d0

rcore(0,36)=  1.7161177d0
rcore(1,36)=  1.9866780d0
rcore(2,36)=  2.0860320d0
rcore(3,36)=  2.0860320d0

rcore(0,38)=  3.4638956d0
rcore(1,38)=  3.9133480d0
rcore(2,38)=  3.4638956d0
rcore(3,38)=  3.4638956d0

rcore(0,39)=  3.1368487d0
rcore(1,39)=  3.1368487d0
rcore(2,39)=  2.6443347d0
rcore(3,39)=  2.6443347d0

rcore(0,40)=  2.9127601d0
rcore(1,40)=  3.0584275d0
rcore(2,40)=  2.5782263d0
rcore(3,40)=  2.5782263d0

rcore(0,41)=  2.7732187d0
rcore(1,41)=  2.9838317d0
rcore(2,41)=  2.4547114d0
rcore(3,41)=  2.4547114d0

rcore(0,42)=  2.6419339d0
rcore(1,42)=  2.7740572d0
rcore(2,42)=  2.3962659d0
rcore(3,42)=  2.3962659d0

rcore(0,43)=  2.5804936d0
rcore(1,43)=  2.6442318d0
rcore(2,43)=  2.6442318d0
rcore(3,43)=  2.6442318d0

rcore(0,44)=  2.4610579d0
rcore(1,44)=  2.5841356d0
rcore(2,44)=  2.4017350d0
rcore(3,44)=  2.4017350d0

rcore(0,45)=  2.4658050d0
rcore(1,45)=  2.5891201d0
rcore(2,45)=  2.3483631d0
rcore(3,45)=  2.3483631d0

rcore(0,46)=  2.4717819d0
rcore(1,46)=  2.5953959d0
rcore(2,46)=  2.3540553d0
rcore(3,46)=  2.3540553d0

rcore(0,47)=  2.4191908d0
rcore(1,47)=  2.6029170d0
rcore(2,47)=  2.4191908d0
rcore(3,47)=  2.4191908d0

rcore(0,48)=  2.3687909d0
rcore(1,48)=  2.4872544d0
rcore(2,48)=  2.3116922d0
rcore(3,48)=  2.3116922d0

rcore(0,49)=  2.4364941d0
rcore(1,49)=  2.4966755d0
rcore(2,49)=  2.3204483d0
rcore(3,49)=  2.3204483d0

rcore(0,50)=  2.3302081d0
rcore(1,50)=  2.3877642d0
rcore(2,50)=  2.6975849d0
rcore(3,50)=  2.6975849d0

rcore(0,51)=  2.2294503d0
rcore(1,51)=  2.3409453d0
rcore(2,51)=  2.6446911d0
rcore(3,51)=  2.6446911d0

rcore(0,52)=  2.2405847d0
rcore(1,52)=  2.2959271d0
rcore(2,52)=  2.6578993d0
rcore(3,52)=  2.6578993d0

rcore(0,53)=  2.1453201d0
rcore(1,53)=  2.1983095d0
rcore(2,53)=  2.6721617d0
rcore(3,53)=  2.6721617d0

rcore(0,54)=  2.0548375d0
rcore(1,54)=  2.0548375d0
rcore(2,54)=  2.4977638d0
rcore(3,54)=  2.4977638d0

rcore(0,55)=  3.9950035d0
rcore(1,55)=  4.4045760d0
rcore(2,55)=  3.9950035d0
rcore(3,55)=  3.9950035d0

rcore(0,56)=  4.2216481d0
rcore(1,56)=  4.8872263d0
rcore(2,56)=  3.9236641d0
rcore(3,56)=  3.9236641d0

rcore(0,58)=  3.4360930d0
rcore(1,58)=  3.9778218d0
rcore(2,58)=  2.9681408d0
rcore(3,58)=  3.1935572d0

rcore(0,59)=  2.4597067d0
rcore(1,59)=  2.4597067d0
rcore(2,59)=  2.4597067d0
rcore(3,59)=  2.4597067d0

rcore(0,61)=  2.4980376d0
rcore(1,61)=  2.4980376d0
rcore(2,61)=  2.4980376d0
rcore(3,61)=  2.4980376d0

rcore(0,62)=  2.4577467d0
rcore(1,62)=  2.4577467d0
rcore(2,62)=  2.4577467d0
rcore(3,62)=  2.4577467d0

rcore(0,63)=  2.4784776d0
rcore(1,63)=  2.4784776d0
rcore(2,63)=  2.4784776d0
rcore(3,63)=  2.4784776d0

rcore(0,64)=  2.4397514d0
rcore(1,64)=  2.4397514d0
rcore(2,64)=  2.4397514d0
rcore(3,64)=  2.4397514d0

rcore(0,65)=  2.4615515d0
rcore(1,65)=  2.4615515d0
rcore(2,65)=  2.4615515d0
rcore(3,65)=  2.4615515d0

rcore(0,66)=  2.4841343d0
rcore(1,66)=  2.4841343d0
rcore(2,66)=  2.4841343d0
rcore(3,66)=  2.4841343d0

rcore(0,68)=  2.4706250d0
rcore(1,68)=  2.4706250d0
rcore(2,68)=  2.4706250d0
rcore(3,68)=  2.4706250d0

rcore(0,72)=  2.7679638d0
rcore(1,72)=  2.8363325d0
rcore(2,72)=  2.4500600d0
rcore(3,72)=  2.4500600d0

rcore(0,73)=  2.6642398d0
rcore(1,73)=  2.7974787d0
rcore(2,73)=  2.4164976d0
rcore(3,73)=  2.4164976d0

rcore(0,74)=  2.6282365d0
rcore(1,74)=  2.6931540d0
rcore(2,74)=  2.4427231d0
rcore(3,74)=  2.4427231d0

rcore(0,75)=  2.5931934d0
rcore(1,75)=  2.6572453d0
rcore(2,75)=  2.4101535d0
rcore(3,75)=  2.4101535d0

rcore(0,76)=  2.5590724d0
rcore(1,76)=  2.6870519d0
rcore(2,76)=  2.4973870d0
rcore(3,76)=  2.4973870d0

rcore(0,77)=  2.5258377d0
rcore(1,77)=  2.6521551d0
rcore(2,77)=  2.5258377d0
rcore(3,77)=  2.5258377d0

rcore(0,78)=  2.4934552d0
rcore(1,78)=  2.6181531d0
rcore(2,78)=  2.4934552d0
rcore(3,78)=  2.4934552d0

rcore(0,79)=  2.5227012d0
rcore(1,79)=  2.6488617d0
rcore(2,79)=  2.5227012d0
rcore(3,79)=  2.5227012d0

rcore(0,80)=  2.4911674d0
rcore(1,80)=  2.4911674d0
rcore(2,80)=  2.3725176d0
rcore(3,80)=  2.3725176d0

rcore(0,81)=  2.4604123d0
rcore(1,81)=  2.4604123d0
rcore(2,81)=  2.4604123d0
rcore(3,81)=  2.4604123d0

rcore(0,82)=  2.3718232d0
rcore(1,82)=  2.3718232d0
rcore(2,82)=  2.6149854d0
rcore(3,82)=  2.6149854d0

rcore(0,83)=  2.2867640d0
rcore(1,83)=  2.2867640d0
rcore(2,83)=  2.4604330d0
rcore(3,83)=  2.4604330d0

rcore(0,84)=  2.4911914d0
rcore(1,84)=  2.4911914d0
rcore(2,84)=  2.4911914d0
rcore(3,84)=  2.4911914d0

rcore(0,85)=  2.4618832d0
rcore(1,85)=  2.4618832d0
rcore(2,85)=  2.4618832d0
rcore(3,85)=  2.4618832d0

rcore(0,86)=  2.4933581d0
rcore(1,86)=  2.4933581d0
rcore(2,86)=  2.4933581d0
rcore(3,86)=  2.4933581d0

rcore(0,96)=  2.4626281d0
rcore(1,96)=  2.4626281d0
rcore(2,96)=  2.4626281d0
rcore(3,96)=  2.4626281d0

rcore(0,97)=  2.4974400d0
rcore(1,97)=  2.4974400d0
rcore(2,97)=  2.4974400d0
rcore(3,97)=  2.4974400d0

rcore(0,99)=  2.4469867d0
rcore(1,99)=  2.4469867d0
rcore(2,99)=  2.4469867d0
rcore(3,99)=  2.4469867d0

rcore(0,100)=  2.4823529d0
rcore(1,100)=  2.4823529d0
rcore(2,100)=  2.4823529d0
rcore(3,100)=  2.4823529d0

rcore(0,103)=  2.4695797d0
rcore(1,103)=  2.4695797d0
rcore(2,103)=  2.4695797d0
rcore(3,103)=  2.4695797d0

rcore(0,104)=  2.4458337d0
rcore(1,104)=  2.4458337d0
rcore(2,104)=  2.4458337d0
rcore(3,104)=  2.4458337d0

rcore(0,105)=  2.4823768d0
rcore(1,105)=  2.4823768d0
rcore(2,105)=  2.4823768d0
rcore(3,105)=  2.4823768d0

rcore(0,106)=  2.4589582d0
rcore(1,106)=  2.4589582d0
rcore(2,106)=  2.4589582d0
rcore(3,106)=  2.4589582d0

rcore(0,107)=  2.4961459d0
rcore(1,107)=  2.4961459d0
rcore(2,107)=  2.4961459d0
rcore(3,107)=  2.4961459d0

rcore(0,108)=  2.4730334d0
rcore(1,108)=  2.4730334d0
rcore(2,108)=  2.4730334d0
rcore(3,108)=  2.4730334d0

rcore(0,109)=  2.4503450d0
rcore(1,109)=  2.4503450d0
rcore(2,109)=  2.4503450d0
rcore(3,109)=  2.4503450d0

select case( iZatom(ak) )
  case(37,57,60,67,69-71,87-95,98,101,102) ; write(*,*) "fhi file is not prepared yet." ; stop
end select

rRC=0.d0
do ll=0,Mlps(ak)
  rRC(ll) = rcore(ll,iZatom(ak))
end do

Rps(ak)=maxval(rRC)
if(maxval( rRC ) > rmaxRps) rmaxRps=maxval( rRC )

End Subroutine setRps_fhi

end module read_pslfile_sub
