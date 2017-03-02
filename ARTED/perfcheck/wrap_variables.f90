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
module wrap_variables
  implicit none

  public wrap_init

#ifdef ARTED_STENCIL_ORIGIN
  private init_original
#endif

contains
  subroutine wrap_init
    use global_variables
    use opt_variables
    implicit none
    integer :: ix,iy,iz,i,ik,ib

    ! wrap global_variables
    Nd = 4
    NL  = NLx*NLy*NLz
    NK_s = 1
    NK_e = NK
    dt = 0.02d0
    functional = 'PZ'
    NKB = (NK_e - NK_s + 1) * NBoccmax
    Nlma = 0

    allocate(zu(NL,NBoccmax,NK_s:NK_e))
    allocate(kAc(NK,3))
    allocate(Vloc(NL))
    allocate(lapx(-Nd:Nd),lapy(-Nd:Nd),lapz(-Nd:Nd))
    allocate(nabx(-Nd:Nd),naby(-Nd:Nd),nabz(-Nd:Nd))
    allocate(ik_table(NKB),ib_table(NKB))

    kAc(:,:) = 0.5d0
    Vloc(:) = 0.33333d0

    lapx(:) = 1.0d0
    lapy(:) = 1.0d0
    lapz(:) = 1.0d0
    lapt(:) = 1.0d0

    nabx(:) = 0.25d0
    naby(:) = 0.25d0
    nabz(:) = 0.25d0

    i=0
    do ik=NK_s,NK_e
      do ib=1,NBoccmax
        i=i+1
        ik_table(i)=ik
        ib_table(i)=ib
      end do
    end do

#ifdef ARTED_STENCIL_ORIGIN
    call init_original
#endif

    ! wrap opt_variables
    PNLx = NLx
#ifdef ARTED_STENCIL_PADDING
    PNLy = NLy + 1
#else
    PNLy = NLy
#endif
    PNLz = NLz
    PNL = PNLx*PNLy*PNLz

    allocate(zhtpsi(0:PNL-1,4,0:NUMBER_THREADS-1))
    allocate(zcx(NBoccmax,NK_s:NK_e))
    allocate(zcy(NBoccmax,NK_s:NK_e))
    allocate(zcz(NBoccmax,NK_s:NK_e))
    allocate(modx(0:NLx*2+Nd-1))
    allocate(mody(0:NLy*2+Nd-1))
    allocate(modz(0:NLz*2+Nd-1))

    do ix=0,NLx*2+Nd-1
      modx(ix) = mod(ix,NLx)
    end do
    do iy=0,NLy*2+Nd-1
      mody(iy) = mod(iy,NLy)
    end do
    do iz=0,NLz*2+Nd-1
      modz(iz) = mod(iz,NLz)
    end do

#ifdef ARTED_STENCIL_ENABLE_LOOP_BLOCKING
    call auto_blocking
#endif
  end subroutine

#ifdef ARTED_STENCIL_ORIGIN
  subroutine init_original
    use global_variables
    use opt_variables
    implicit none
    integer :: ix,iy,iz,i,j

    allocate(zifdx(-4:4,0:NL-1))
    allocate(zifdy(-4:4,0:NL-1))
    allocate(zifdz(-4:4,0:NL-1))

    do ix=0,NLx-1
    do iy=0,NLy-1
    do iz=0,NLz-1
      i = ix*NLy*NLz + iy*NLz + iz
      do j=1,4
        zifdx( j, i) = mod(ix+j+NLx,NLx)*NLy*NLz + iy*NLz + iz
        zifdx(-j, i) = mod(ix-j+NLx,NLx)*NLy*NLz + iy*NLz + iz
      end do
      do j=1,4
        zifdy( j, i) = ix*NLy*NLz + mod(iy+j+NLy,NLy)*NLz + iz
        zifdy(-j, i) = ix*NLy*NLz + mod(iy-j+NLy,NLy)*NLz + iz
      end do
      do j=1,4
        zifdz( j, i) = ix*NLy*NLz + iy*NLz + mod(iz+j+NLz,NLz)
        zifdz(-j, i) = ix*NLy*NLz + iy*NLz + mod(iz-j+NLz,NLz)
      end do
    end do
    end do
    end do
  end subroutine
#endif
end module
