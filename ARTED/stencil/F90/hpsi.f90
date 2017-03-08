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
#define ENABLE_NONTEMPORAL_STORE

#if defined(__KNC__) || defined(__AVX512F__) || defined(__HPC_ACE2__)
# define ENABLE_OPTIMIZED_LOAD
#endif

subroutine hpsi1_RT_stencil(A,B,C,D,E,F)
  use global_variables, only: NLx,NLy,NLz,zI
#ifdef ARTED_STENCIL_ENABLE_LOOP_BLOCKING
  use opt_variables, only: LBX => STENCIL_BLOCKING_X, LBY => STENCIL_BLOCKING_Y
#endif
#ifndef ARTED_DOMAIN_POWER_OF_TWO
  use opt_variables, only: modx, mody, modz
#endif
  use opt_variables, only: PNLx,PNLy,PNLz
  implicit none
  real(8),   intent(in)  :: A
  real(8),   intent(in)  :: C(12), D(12)
  real(8),   intent(in)  :: B(0:NLz-1,0:NLy-1,0:NLx-1)
  complex(8),intent(in)  :: E(0:PNLz-1,0:PNLy-1,0:PNLx-1)
  complex(8),intent(out) :: F(0:PNLz-1,0:PNLy-1,0:PNLx-1)

#ifdef ARTED_STENCIL_ENABLE_LOOP_BLOCKING
  integer    :: bx,by
#endif
  integer    :: ix,iy,iz
  complex(8) :: v, w
#ifdef ENABLE_OPTIMIZED_LOAD
  complex(8) :: t(8)
#endif

#ifdef __INTEL_COMPILER
#if defined(__KNC__) || defined(__AVX512F__)
#   define MEM_ALIGN   64
#   define VECTOR_SIZE 4
# else
#   define MEM_ALIGN   32
#   define VECTOR_SIZE 2
# endif

!dir$ assume_aligned B:MEM_ALIGN
!dir$ assume_aligned E:MEM_ALIGN
!dir$ assume_aligned F:MEM_ALIGN
#endif

#ifdef ARTED_DOMAIN_POWER_OF_TWO
# ifdef __INTEL_COMPILER
!dir$ assume (mod(NLx, VECTOR_SIZE) == 0)
!dir$ assume (mod(NLy, VECTOR_SIZE) == 0)
!dir$ assume (mod(NLz, VECTOR_SIZE) == 0)
# endif
# define IDX(dt) iz,iy,and(ix+(dt)+NLx,NLx-1)
# define IDY(dt) iz,and(iy+(dt)+NLy,NLy-1),ix
# define IDZ(dt) and(iz+(dt)+NLz,NLz-1),iy,ix
#else
# define IDX(dt) iz,iy,modx(ix+(dt)+NLx)
# define IDY(dt) iz,mody(iy+(dt)+NLy),ix
# define IDZ(dt) modz(iz+(dt)+NLz),iy,ix
#endif


#ifdef ARTED_STENCIL_ENABLE_LOOP_BLOCKING
  do bx=0,NLx-1,LBX
  do by=0,NLy-1,LBY
  do ix=bx,min(bx+LBX-1,NLx-1)
  do iy=by,min(by+LBY-1,NLy-1)
#else
  do ix=0,NLx-1
  do iy=0,NLy-1
#endif
#ifdef __INTEL_COMPILER
!dir$ simd
#ifdef ENABLE_NONTEMPORAL_STORE
!dir$ vector nontemporal(F)
#endif
#endif
#ifdef __FUJITSU
!OCL simd
!OCL noalias
#endif
  do iz=0,NLz-1

!dir$ assume_aligned B(0,iy,ix):MEM_ALIGN
!dir$ assume_aligned E(0,iy,ix):MEM_ALIGN
!dir$ assume_aligned F(0,iy,ix):MEM_ALIGN

#ifdef ENABLE_OPTIMIZED_LOAD
    t(1) = E(IDZ( 1))
    t(2) = E(IDZ( 2))
    t(3) = E(IDZ( 3))
    t(4) = E(IDZ( 4))
    t(5) = E(IDZ(-1))
    t(6) = E(IDZ(-2))
    t(7) = E(IDZ(-3))
    t(8) = E(IDZ(-4))

    v=(C( 9)*(t(1)+t(5)) &
    & +C(10)*(t(2)+t(6)) &
    & +C(11)*(t(3)+t(7)) &
    & +C(12)*(t(4)+t(8)))
    w=(D( 9)*(t(1)-t(5)) &
    & +D(10)*(t(2)-t(6)) &
    & +D(11)*(t(3)-t(7)) &
    & +D(12)*(t(4)-t(8)))

    t(1) = E(IDY( 1))
    t(2) = E(IDY( 2))
    t(3) = E(IDY( 3))
    t(4) = E(IDY( 4))
    t(5) = E(IDY(-1))
    t(6) = E(IDY(-2))
    t(7) = E(IDY(-3))
    t(8) = E(IDY(-4))

    v=(C(5)*(t(1)+t(5)) &
    & +C(6)*(t(2)+t(6)) &
    & +C(7)*(t(3)+t(7)) &
    & +C(8)*(t(4)+t(8))) + v
    w=(D(5)*(t(1)-t(5)) &
    & +D(6)*(t(2)-t(6)) &
    & +D(7)*(t(3)-t(7)) &
    & +D(8)*(t(4)-t(8))) + w

    t(1) = E(IDX( 1))
    t(2) = E(IDX( 2))
    t(3) = E(IDX( 3))
    t(4) = E(IDX( 4))
    t(5) = E(IDX(-1))
    t(6) = E(IDX(-2))
    t(7) = E(IDX(-3))
    t(8) = E(IDX(-4))

    v=(C(1)*(t(1)+t(5)) &
    & +C(2)*(t(2)+t(6)) &
    & +C(3)*(t(3)+t(7)) &
    & +C(4)*(t(4)+t(8))) + v
    w=(D(1)*(t(1)-t(5)) &
    & +D(2)*(t(2)-t(6)) &
    & +D(3)*(t(3)-t(7)) &
    & +D(4)*(t(4)-t(8))) + w
#else
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
#endif

    F(iz,iy,ix) = B(iz,iy,ix)*E(iz,iy,ix) + A*E(iz,iy,ix) - 0.5d0 * v - zI * w
  end do
  end do
  end do
#ifdef ARTED_STENCIL_ENABLE_LOOP_BLOCKING
  end do
  end do
#endif
end subroutine
