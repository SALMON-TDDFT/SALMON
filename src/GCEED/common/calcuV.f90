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
