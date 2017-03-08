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
subroutine current_stencil_LBLK(E, ikb_s,ikb_e)
  use global_variables, only: ik_table,ib_table,NBoccmax,NK_s,NK_e,NLx,NLy,NLz, &
  &                           nabx,naby,nabz,zI
  use opt_variables
  implicit none
  complex(8), intent(in)  :: E(0:NLz-1,0:NLy-1,0:NLx-1, NBoccmax, NK_s:NK_e)
  integer :: ikb_s,ikb_e

  real(8)    :: F,G,H
  integer    :: ix,iy,iz
  complex(8) :: v,w
  integer    :: ikb,ik,ib

#undef IDX
#undef IDY
#undef IDZ
#ifdef ARTED_DOMAIN_POWER_OF_TWO
# define IDX(dt) iz,iy,and(ix+(dt)+NLx,NLx-1),ib,ik
# define IDY(dt) iz,and(iy+(dt)+NLy,NLy-1),ix,ib,ik
# define IDZ(dt) and(iz+(dt)+NLz,NLz-1),iy,ix,ib,ik
#else
# define IDX(dt) iz,iy,modx(ix+(dt)+NLx),ib,ik
# define IDY(dt) iz,mody(iy+(dt)+NLy),ix,ib,ik
# define IDZ(dt) modz(iz+(dt)+NLz),iy,ix,ib,ik
#endif

!$acc kernels pcopy(zcx,zcy,zcz) &
#ifndef ARTED_DOMAIN_POWER_OF_TWO
!$acc pcopyin(modx,mody,modz) &
#endif
!$acc pcopyin(E,ib_table,ik_table,nabx,naby,nabz)
!$acc loop independent gang private(H,G,F)
  do ikb=ikb_s,ikb_e
    ik=ik_table(ikb)
    ib=ib_table(ikb)

    F = 0
!$acc loop collapse(3) vector(128) reduction(+:F)
    do iy=0,NLy-1
    do ix=0,NLx-1
    do iz=0,NLz-1
      w = conjg(E(iz,iy,ix, ib,ik))
      v=(nabx(1)*(E(IDX(1))) &
      & +nabx(2)*(E(IDX(2))) &
      & +nabx(3)*(E(IDX(3))) &
      & +nabx(4)*(E(IDX(4))))
      F = F + imag(w * v)
    end do
    end do
    end do
    zcx(ib,ik)=F * 2.d0

    G = 0
!$acc loop collapse(3) vector(128) reduction(+:G)
    do ix=0,NLx-1
    do iy=0,NLy-1
    do iz=0,NLz-1
      w = conjg(E(iz,iy,ix, ib,ik))
      v=(naby(1)*(E(IDY(1))) &
      & +naby(2)*(E(IDY(2))) &
      & +naby(3)*(E(IDY(3))) &
      & +naby(4)*(E(IDY(4))))
      G = G + imag(w * v)
    end do
    end do
    end do
    zcy(ib,ik)=G * 2.d0

    H = 0
!$acc loop collapse(3) vector(128) reduction(+:H)
    do ix=0,NLx-1
    do iy=0,NLy-1
    do iz=0,NLz-1
      w = conjg(E(iz,iy,ix, ib,ik))
      v=(nabz(1)*(E(IDZ(1))) &
      & +nabz(2)*(E(IDZ(2))) &
      & +nabz(3)*(E(IDZ(3))) &
      & +nabz(4)*(E(IDZ(4))))
      H = H + imag(w * v)
    end do
    end do
    end do
    zcz(ib,ik)=H * 2.d0
  end do
!$acc end kernels
end subroutine
