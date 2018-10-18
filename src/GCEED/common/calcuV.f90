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
SUBROUTINE calcuV
use salmon_parallel, only: nproc_id_global
use salmon_communication, only: comm_is_root
use scf_data
use allocate_psl_sub
implicit none
integer :: ikoa,LL,ii,iatom,jj,lm
real(8) :: xx,yy,zz,ratio1,ratio2,url,uVrl,rr
real(8) :: Ylm
integer :: intr

  integer :: a,i,ik,ix,iy,iz,l,m
  integer :: nl

  real(8) :: alx,aly,alz
  real(8) :: hx,hy,hz
  integer :: lx(lg_num(1)*lg_num(2)*lg_num(3))
  integer :: ly(lg_num(1)*lg_num(2)*lg_num(3))
  integer :: lz(lg_num(1)*lg_num(2)*lg_num(3))

  integer :: nlma
  integer :: lma
  character(17) :: property
  
  real(8),allocatable :: save_udVtbl_a(:,:,:)
  real(8),allocatable :: save_udVtbl_b(:,:,:)
  real(8),allocatable :: save_udVtbl_c(:,:,:)
  real(8),allocatable :: save_udVtbl_d(:,:,:)

  integer,allocatable :: lma_tbl(:,:)

  real(8),allocatable :: uv_tmp(:,:),duv_tmp(:,:,:)

  logical :: flag_use_grad_wf_on_force
  

  property='initial'
  flag_use_grad_wf_on_force=.false. 

  nl=lg_num(1)*lg_num(2)*lg_num(3)

  hx=Hgs(1)
  hy=Hgs(2)
  hz=Hgs(3)
  alx=Hgs(1)*dble(lg_num(1))
  aly=Hgs(2)*dble(lg_num(2))
  alz=Hgs(3)*dble(lg_num(3))

  do iz=1,lg_num(3)
  do iy=1,lg_num(2)
  do ix=1,lg_num(1)
    i=(iz-1)*lg_num(1)*lg_num(2)+(iy-1)*lg_num(1)+ix
    lx(i)=ix-1
    ly(i)=iy-1
    lz(i)=iz-1
  end do
  end do
  end do

  lma=0
  do a=1,MI
    ik=Kion(a)
    do l=0,Mlps(ik)
      if(pp%inorm(l,ik)==0) cycle
      do m=-l,l
        lma=lma+1
      enddo
    enddo
  enddo
  nlma=lma

  allocate(lma_tbl((pp%lmax+1)**2,MI))

  lma=0
  do a=1,MI
    ik=Kion(a)
    lm=0
    do l=0,Mlps(ik)
      if(pp%inorm(l,ik)==0) cycle
      do m=-l,l
        lm=lm+1
        lma=lma+1
        lma_tbl(lm,a)=lma
      enddo
    enddo
  enddo


  allocate( save_udVtbl_a(pp%nrmax,0:pp%lmax,natom) )
  allocate( save_udVtbl_b(pp%nrmax,0:pp%lmax,natom) )
  allocate( save_udVtbl_c(pp%nrmax,0:pp%lmax,natom) )
  allocate( save_udVtbl_d(pp%nrmax,0:pp%lmax,natom) )
     
!  allocate(a_tbl(Nlma),uV(Nps,Nlma),iuV(Nlma),duV(Nps,Nlma,3))
  allocate(uV_tmp(ppg%nps,nlma),duV_tmp(ppg%nps,nlma,3))

call calc_uv(pp,ppg,save_udvtbl_a,save_udvtbl_b,save_udvtbl_c,save_udvtbl_d,uv_tmp,duv_tmp, &
             nlma,lx,ly,lz,nl,hx,hy,hz,alx,aly,alz,  &
             lma_tbl,flag_use_grad_wf_on_force,property)

!  uV=0.d0
!  uVu=1.d0
!
!  lma = 0
!  do iatom=1,MI
!    ik=Kion(iatom)
!    do l=0,Mlps(ik)
!      if ( pp%inorm(l,ik)==0) cycle
!      do lm=l**2+1,(l+1)**2
!        lma = lma + 1
!        do jj=1,ppg%mps(iatom)
!          uV(jj,lm,iatom) = uv_tmp(jj,lma)
!        end do
!        uVu(lm,iatom)=dble(pp%inorm(l,ik))
!      end do 
!    end do
!  end do
!


do ikoa=1,MKI
  do LL=0,Mlps(ikoa)
    do ii=0,Nr
      uVnl(ii,LL,ikoa)=(vpp(ii,LL,ikoa)-vpp(ii,Lref(ikoa),ikoa))*uppr(ii,LL,ikoa)
    end do
  end do
end do

do iatom=1,MI
  ikoa=Kion(iatom)
  do LL=0,Mlps(ikoa)
  do jj=1,Mps_all(iatom)
    xx=gridcoo(Jxyz_all(1,jj,iatom),1)-Jxxyyzz_all(1,jj,iatom)*lg_num(1)*Hgs(1)-Rion(1,iatom)
    yy=gridcoo(Jxyz_all(2,jj,iatom),2)-Jxxyyzz_all(2,jj,iatom)*lg_num(2)*Hgs(2)-Rion(2,iatom)
    zz=gridcoo(Jxyz_all(3,jj,iatom),3)-Jxxyyzz_all(3,jj,iatom)*lg_num(3)*Hgs(3)-Rion(3,iatom)
    rr=sqrt(xx**2+yy**2+zz**2)+1d-50
    call bisection(rr,intr,ikoa)
    if(intr>Nr) write(*,*) "intr is larger than Nr"
    ratio1=(rr-rad_psl(intr,ikoa))/(rad_psl(intr+1,ikoa)-rad_psl(intr,ikoa)) ; ratio2=1.d0-ratio1
    url =ratio1*uppr(intr+1,LL,ikoa)+ratio2*uppr(intr,LL,ikoa)
    uVrl=ratio1*uVnl(intr+1,LL,ikoa)+ratio2*uVnl(intr,LL,ikoa)
    do lm=LL**2+1,(LL+1)**2
      call Ylm_sub(xx,yy,zz,lm,Ylm)
      ur(jj,lm)=url*Ylm
      uV_all(jj,lm,iatom)=uVrl*Ylm
    end do
  end do
  end do
  
  do lm=1,(Mlps(ikoa)+1)**2
    uVu(lm,iatom)=sum(uV_all(:Mps_all(iatom),lm,iatom)*ur(:Mps_all(iatom),lm))*Hvol
    if ( abs(uVu(lm,iatom)) < 1.d-10 ) uVu(lm,iatom)=1.d-10
  end do
end do
  
if(iSCFRT==1)then
  if(comm_is_root(nproc_id_global))then
    do iatom=1,MI
      ikoa=Kion(iatom)
      write(*,'(1x,"ion = ",i5,"  uVu integral table")') iatom
      write(*,'(1x,4f15.5)')       &
     (uVu(lm,iatom)*2d0*Ry*(a_B),lm=1,(Mlps(ikoa)+1)**2)
    end do
  end if
end if

return

END SUBROUTINE calcuV
