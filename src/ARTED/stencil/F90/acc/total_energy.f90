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

subroutine total_energy_stencil_acc(A,C,E,zu_NB,Ekin_l)
  use global_variables, only: NLx,NLy,NLz,zI,Nboccmax,NK_s,NK_e,nabx,naby,nabz, &
                           &  occ,Hxyz,kAc
  use opt_variables, only: zJxyz
#ifndef ARTED_DOMAIN_POWER_OF_TWO
  use opt_variables, only: modx, mody, modz
#endif
  implicit none
  real(8),   intent(in)  :: A
  real(8),   intent(in)  :: C(12)
  complex(8),intent(in)  :: E(0:NLz-1,0:NLy-1,0:NLx-1,zu_NB,NK_s:NK_e)
  integer,   intent(in)  :: zu_NB
  real(8),   intent(out) :: Ekin_l

  integer    :: ix,iy,iz,ib,ik
  complex(8) :: v,w,z
  complex(8) :: F
  real(8)    :: D(12)

#ifdef ARTED_DOMAIN_POWER_OF_TWO
# define IDX(dt) iz,iy,and(ix+(dt)+NLx,NLx-1),ib,ik
# define IDY(dt) iz,and(iy+(dt)+NLy,NLy-1),ix,ib,ik
# define IDZ(dt) and(iz+(dt)+NLz,NLz-1),iy,ix,ib,ik
#else
# define IDX(dt) iz,iy,modx(ix+(dt)+NLx),ib,ik
# define IDY(dt) iz,mody(iy+(dt)+NLy),ix,ib,ik
# define IDZ(dt) modz(iz+(dt)+NLz),iy,ix,ib,ik
#endif

!$acc kernels
!$acc loop gang collapse(2) private(D,F) reduction(+:Ekin_l) independent
  do ik=NK_s,NK_e
  do ib=1,NBoccmax

    D( 1: 4)=kAc(ik,1)*nabx(1:4)
    D( 5: 8)=kAc(ik,2)*naby(1:4)
    D( 9:12)=kAc(ik,3)*nabz(1:4)

    F = 0
!$acc loop vector(128) collapse(3) private(z,v,w) reduction(+:F)
    do ix=0,NLx-1
    do iy=0,NLy-1
    do iz=0,NLz-1

      z = A * E(iz,iy,ix,ib,ik)

      v=(C( 9)*(E(IDZ(1))+E(IDZ(-1))) &
      & +C(10)*(E(IDZ(2))+E(IDZ(-2))) &
      & +C(11)*(E(IDZ(3))+E(IDZ(-3))) &
      & +C(12)*(E(IDZ(4))+E(IDZ(-4))))
      w=(D( 9)*(E(IDZ(1))-E(IDZ(-1))) &
      & +D(10)*(E(IDZ(2))-E(IDZ(-2))) &
      & +D(11)*(E(IDZ(3))-E(IDZ(-3))) &
      & +D(12)*(E(IDZ(4))-E(IDZ(-4))))

      v=(C( 5)*(E(IDY(1))+E(IDY(-1))) &
      & +C( 6)*(E(IDY(2))+E(IDY(-2))) &
      & +C( 7)*(E(IDY(3))+E(IDY(-3))) &
      & +C( 8)*(E(IDY(4))+E(IDY(-4)))) + v
      w=(D( 5)*(E(IDY(1))-E(IDY(-1))) &
      & +D( 6)*(E(IDY(2))-E(IDY(-2))) &
      & +D( 7)*(E(IDY(3))-E(IDY(-3))) &
      & +D( 8)*(E(IDY(4))-E(IDY(-4)))) + w

      v=(C( 1)*(E(IDX(1))+E(IDX(-1))) &
      & +C( 2)*(E(IDX(2))+E(IDX(-2))) &
      & +C( 3)*(E(IDX(3))+E(IDX(-3))) &
      & +C( 4)*(E(IDX(4))+E(IDX(-4)))) + v
      w=(D( 1)*(E(IDX(1))-E(IDX(-1))) &
      & +D( 2)*(E(IDX(2))-E(IDX(-2))) &
      & +D( 3)*(E(IDX(3))-E(IDX(-3))) &
      & +D( 4)*(E(IDX(4))-E(IDX(-4)))) + w

      z = z - 0.5d0 * v - zI * w
      F = F + conjg(E(iz,iy,ix,ib,ik)) * z
    end do
    end do
    end do

    Ekin_l = Ekin_l + occ(ib,ik) * F * Hxyz &
    &               + occ(ib,ik) * sum(kAc(ik,:)**2) * 0.5d0

  end do
  end do
!$acc end kernels

end subroutine
