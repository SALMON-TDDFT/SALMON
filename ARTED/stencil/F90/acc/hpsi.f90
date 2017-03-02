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
subroutine hpsi1_RT_stencil_LBLK(A,B,C,D,E,F, ikb_s,ikb_e)
  use global_variables, only: NLx,NLy,NLz,zI
#ifndef ARTED_DOMAIN_POWER_OF_TWO
  use opt_variables, only: modx, mody, modz
#endif
  use opt_variables, only: PNLx,PNLy,PNLz
  implicit none
  integer :: ikb_s,ikb_e
  real(8),   intent(in)  :: A(ikb_s:ikb_e)
  real(8),   intent(in)  :: B(0:NLz-1,0:NLy-1,0:NLx-1)
  real(8),   intent(in)  :: C(12)
  real(8),   intent(in)  :: D(12, ikb_s:ikb_e)
  complex(8),intent(in)  :: E(0:PNLz-1,0:PNLy-1,0:PNLx-1, ikb_s:ikb_e)
  complex(8),intent(out) :: F(0:PNLz-1,0:PNLy-1,0:PNLx-1, ikb_s:ikb_e)

  integer    :: ikb, ix,iy,iz
  complex(8) :: v, w

#ifdef ARTED_DOMAIN_POWER_OF_TWO
# define IDX(dt) iz,iy,and(ix+(dt)+NLx,NLx-1),ikb
# define IDY(dt) iz,and(iy+(dt)+NLy,NLy-1),ix,ikb
# define IDZ(dt) and(iz+(dt)+NLz,NLz-1),iy,ix,ikb
#else
# define IDX(dt) iz,iy,modx(ix+(dt)+NLx),ikb
# define IDY(dt) iz,mody(iy+(dt)+NLy),ix,ikb
# define IDZ(dt) modz(iz+(dt)+NLz),iy,ix,ikb
#endif

!$acc kernels pcopy(F) &
#ifndef ARTED_DOMAIN_POWER_OF_TWO
!$acc& pcopyin(modx,mody,modz) &
#endif
!$acc& pcopyin(A,B,C,D,E)
!$acc loop gang vector(1)
  do ikb = ikb_s, ikb_e
!$acc loop collapse(2) independent gang vector(128)
    do iy=0,NLy-1
    do iz=0,NLz-1
!$acc loop seq
    do ix=0,NLx-1
      v=(C( 9)    *(E(IDZ(1))+E(IDZ(-1))) &
      & +C(10)    *(E(IDZ(2))+E(IDZ(-2))) &
      & +C(11)    *(E(IDZ(3))+E(IDZ(-3))) &
      & +C(12)    *(E(IDZ(4))+E(IDZ(-4))))
      w=(D( 9,ikb)*(E(IDZ(1))-E(IDZ(-1))) &
      & +D(10,ikb)*(E(IDZ(2))-E(IDZ(-2))) &
      & +D(11,ikb)*(E(IDZ(3))-E(IDZ(-3))) &
      & +D(12,ikb)*(E(IDZ(4))-E(IDZ(-4))))
  
      v=(C( 5)    *(E(IDY(1))+E(IDY(-1))) &
      & +C( 6)    *(E(IDY(2))+E(IDY(-2))) &
      & +C( 7)    *(E(IDY(3))+E(IDY(-3))) &
      & +C( 8)    *(E(IDY(4))+E(IDY(-4)))) + v
      w=(D( 5,ikb)*(E(IDY(1))-E(IDY(-1))) &
      & +D( 6,ikb)*(E(IDY(2))-E(IDY(-2))) &
      & +D( 7,ikb)*(E(IDY(3))-E(IDY(-3))) &
      & +D( 8,ikb)*(E(IDY(4))-E(IDY(-4)))) + w
  
      v=(C( 1)    *(E(IDX(1))+E(IDX(-1))) &
      & +C( 2)    *(E(IDX(2))+E(IDX(-2))) &
      & +C( 3)    *(E(IDX(3))+E(IDX(-3))) &
      & +C( 4)    *(E(IDX(4))+E(IDX(-4)))) + v
      w=(D( 1,ikb)*(E(IDX(1))-E(IDX(-1))) &
      & +D( 2,ikb)*(E(IDX(2))-E(IDX(-2))) &
      & +D( 3,ikb)*(E(IDX(3))-E(IDX(-3))) &
      & +D( 4,ikb)*(E(IDX(4))-E(IDX(-4)))) + w
  
      F(iz,iy,ix, ikb) = (A(ikb)+B(iz,iy,ix))*E(iz,iy,ix, ikb) &
        - 0.5d0 * v - zI * w
    end do
    end do
    end do
  end do
!$acc end kernels
end subroutine
