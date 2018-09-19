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
module opt_variables
  implicit none

  integer :: NUMBER_THREADS_POW2

  real(8) :: lapt(12)

  integer                :: PNLx,PNLy,PNLz,PNL
  complex(8),allocatable :: zhtpsi(:,:,:),zttpsi(:,:)
  complex(8),allocatable :: ghtpsi(:,:,:)
  complex(8),allocatable :: spseudo(:,:),dpseudo(:,:)  ! (NPI, # of threads)

  real(8),allocatable :: zrhotmp(:,:)

  integer,allocatable :: zJxyz(:,:),zKxyz(:,:)

  integer,allocatable :: modx(:),mody(:),modz(:)

  integer :: STENCIL_BLOCKING_X
  integer :: STENCIL_BLOCKING_Y

  real(8),allocatable :: zcx(:,:),zcy(:,:),zcz(:,:)

  integer,allocatable    :: nprojector(:)       ! # of projector
  integer                :: NPI                 ! size of pseudo-vector (packed vector)
  integer,allocatable    :: idx_proj(:)         ! projector element index
  integer,allocatable    :: idx_lma(:)          ! start index of lma
  integer,allocatable    :: pseudo_start_idx(:) ! start index of pseudo-vector

#ifdef ARTED_STENCIL_ORIGIN
  integer,allocatable :: zifdx(:,:),zifdy(:,:),zifdz(:,:)
#endif

#ifdef ARTED_LBLK
  integer,allocatable :: t4ppt_nlma(:)    ! (PNL)
  integer,allocatable :: t4ppt_i2vi(:)    ! (PNL)
  integer,allocatable :: t4ppt_vi2i(:)    ! (PNL)
  integer,allocatable :: t4ppt_ilma(:,:)  ! (PNL?,Nlma?)
  integer,allocatable :: t4ppt_j(:,:)     ! (PNL?,Nlma?)
  integer :: t4ppt_max_vi

  integer, parameter :: at_least_parallelism = 4*1024*1024
  integer :: blk_nkb_hpsi
  integer :: blk_nkb_current

  real(8),allocatable :: t4cp_uVpsix(:,:)  ! (Nlma, NKB)
  real(8),allocatable :: t4cp_uVpsiy(:,:)
  real(8),allocatable :: t4cp_uVpsiz(:,:)
#endif

#if defined(__KNC__) || defined(__AVX512F__)
# define MEM_ALIGNED 64
#else
# define MEM_ALIGNED 32
#endif

!dir$ attributes align:MEM_ALIGNED :: lapt
!dir$ attributes align:MEM_ALIGNED :: zhtpsi,zttpsi
!dir$ attributes align:MEM_ALIGNED :: zrhotmp
!dir$ attributes align:MEM_ALIGNED :: zJxyz,zKxyz
!dir$ attributes align:MEM_ALIGNED :: zcx,zcy,zcz
!dir$ attributes align:MEM_ALIGNED :: modx,mody,modz

#ifdef ARTED_STENCIL_ORIGIN
!dir$ attributes align:MEM_ALIGNED :: zifdx,zifdy,zifdz
#endif

contains
  subroutine opt_vars_initialize_p1
    use global_variables
    use misc_routines, only: ceiling_pow2
    implicit none
    integer :: tid_range

    select case(functional)
      case('TPSS','VS98')
        call err_finalize('functional: TPSS/VS98 versions not implemented.')
    end select

    NUMBER_THREADS_POW2 = ceiling_pow2(NUMBER_THREADS)
#ifdef ARTED_REDUCE_FOR_MANYCORE
    tid_range = NUMBER_THREADS_POW2 - 1
#else
    tid_range = 0
#endif
    allocate(zrhotmp(0:NL-1,0:tid_range))
    zrhotmp(:,:) = 0.0d0
  end subroutine

  subroutine opt_vars_initialize_p2
    use global_variables
    implicit none

    integer :: ix,iy,iz

    PNLx = NLx
#ifdef ARTED_STENCIL_PADDING
    PNLy = NLy + 1
#else
    PNLy = NLy
#endif
    PNLz = NLz
    PNL  = PNLx * PNLy * PNLz

    allocate(zhtpsi(0:PNL-1,4,0:NUMBER_THREADS-1))
#ifdef ARTED_LBLK
    blk_nkb_hpsi = min(at_least_parallelism/PNL + 1, NKB)
    allocate(ghtpsi(0:PNL-1, 0:blk_nkb_hpsi-1, 4))
    !write(*,*) "blk_nkb_hpsi:", blk_nkb_hpsi

    !blk_nkb_current = min(at_least_parallelism/PNL + 1, NKB)
    blk_nkb_current = min(at_least_parallelism/(Nlma*128) + 1, NKB)
    !write(*,*) "blk_nkb_current:", blk_nkb_current
#endif
    allocate(zttpsi(0:PNL-1,0:NUMBER_THREADS-1))

    allocate(zcx(NBoccmax,NK_s:NK_e))
    allocate(zcy(NBoccmax,NK_s:NK_e))
    allocate(zcz(NBoccmax,NK_s:NK_e))

    lapt( 1: 4)=lapx(1:4)
    lapt( 5: 8)=lapy(1:4)
    lapt( 9:12)=lapz(1:4)

#ifdef ARTED_STENCIL_ORIGIN
    allocate(zifdx(-4:4,0:NL-1))
    allocate(zifdy(-4:4,0:NL-1))
    allocate(zifdz(-4:4,0:NL-1))

    zifdx(-4:4,0:NL-1) = ifdx(-4:4,1:NL) - 1
    zifdy(-4:4,0:NL-1) = ifdy(-4:4,1:NL) - 1
    zifdz(-4:4,0:NL-1) = ifdz(-4:4,1:NL) - 1
#endif

    if(.not. allocated(zJxyz)) allocate(zJxyz(Nps,NI))  !AY see subroutine prep_ps_periodic
    !allocate(zJxyz(Nps,NI))

    zJxyz(1:Nps,1:NI) = Jxyz(1:Nps,1:NI) - 1

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

#ifdef ARTED_STENCIL_PADDING
    call init_for_padding
    call init_projector(zKxyz)
#else
    call init_projector(zJxyz)
#endif

#ifdef ARTED_STENCIL_ENABLE_LOOP_BLOCKING
    call auto_blocking
#endif
  end subroutine

  subroutine init_for_padding
    use global_variables
    implicit none
    integer :: a,ik,ix,iy,iz,jx,jy,jz,i,j,k
    real(8) :: x,y,z,r

    if(allocated(zKxyz)) deallocate(zKxyz)
    allocate(zKxyz(Nps,NI))

    do a=1,NI
      ik=Kion(a)
      j=0
      do ix=-2,2
      do iy=-2,2
      do iz=-2,2
        do jx=0,NLx-1
        do jy=0,NLy-1
        do jz=0,NLz-1
          i=jx* NLy* NLz + jy* NLz + jz + 1
          k=jx*PNLy*PNLz + jy*PNLz + jz
          x=Lx(i)*Hx-(Rion(1,a)+ix*aLx)
          y=Ly(i)*Hy-(Rion(2,a)+iy*aLy)
          z=Lz(i)*Hz-(Rion(3,a)+iz*aLz)
          r=sqrt(x*x+y*y+z*z)
          if (r<Rps(ik)) then
            j=j+1
            if(j<=Nps) then
              zKxyz(j,a)=k
            end if
          end if
        end do
        end do
        end do
      end do
      end do
      end do
    end do
  end subroutine

  function count_if_integer(vec, val) result(n)
    implicit none
    integer,intent(in) :: vec(:)
    integer,intent(in) :: val
    integer :: n, i
    n = 0
    do i=1,size(vec)
      if (vec(i) == val) then
        n = n + 1
      end if
    end do
  end function

  subroutine init_projector(zJxyz)
    use global_variables
    implicit none
    integer,intent(in) :: zJxyz(Nps,NI)
    integer :: i,ioffset

    NPI = sum(Mps)

    allocate(nprojector(NI),idx_lma(NI))

    do i=1,NI
      nprojector(i) = count_if_integer(a_tbl, i)
    end do

    idx_lma(1) = 0
    do i=2,NI
      idx_lma(i) = idx_lma(i-1) + nprojector(i-1)
    end do

    allocate(idx_proj(NPI))
    allocate(pseudo_start_idx(NI))

    ioffset = 0
    do i=1,NI
      pseudo_start_idx(i) = ioffset
      idx_proj(ioffset+1:ioffset+Mps(i)) = zJxyz(1:Mps(i),i)
      ioffset = ioffset + Mps(i)
    end do

    if(allocated(spseudo)) deallocate(spseudo,dpseudo)
    allocate(spseudo(NPI,0:NUMBER_THREADS-1))
    allocate(dpseudo(NPI,0:NUMBER_THREADS-1))
  end subroutine


#ifdef ARTED_LBLK
  subroutine opt_vars_init_t4ppt
    use global_variables
    implicit none

    integer    :: ilma,ia,j,i, max_nlma,n, vi,max_vi
    ! write(*,*) "NUMBER_THREADS:", NUMBER_THREADS

    allocate(t4ppt_nlma(0:PNL-1))
    allocate(t4ppt_i2vi(0:PNL-1))
    allocate(t4ppt_vi2i(0:PNL-1))

    t4ppt_nlma(:) = 0
    do ilma=1,Nlma
       ia=a_tbl(ilma)
       do j=1,Mps(ia)
#ifdef ARTED_STENCIL_PADDING
          i=zKxyz(j,ia)
#else
          i=zJxyz(j,ia)
#endif
          t4ppt_nlma(i) = t4ppt_nlma(i) + 1
       enddo
    enddo

    max_nlma = 0
    vi = 0
    do i=0,PNL-1
       max_nlma = max(max_nlma, t4ppt_nlma(i))

       t4ppt_i2vi(i) = -1
       if (t4ppt_nlma(i) > 0) then
          t4ppt_i2vi( i) = vi
          t4ppt_vi2i(vi) =  i
          vi = vi + 1
       endif
    enddo
    max_vi = vi
    ! write(*,*) "max_nlma:", max_nlma
    ! write(*,*) "max_vi:", max_vi

    allocate(t4ppt_ilma(0:max_vi-1, max_nlma))
    allocate(t4ppt_j   (0:max_vi-1, max_nlma))
    t4ppt_max_vi = max_vi
    t4ppt_nlma(:) = 0
    do ilma=1,Nlma
       ia=a_tbl(ilma)
       do j=1,Mps(ia)
#ifdef ARTED_STENCIL_PADDING
          i=zKxyz(j,ia)
#else
          i=zJxyz(j,ia)
#endif
          vi = t4ppt_i2vi(i)

          t4ppt_nlma(vi) = t4ppt_nlma(vi) + 1
          n = t4ppt_nlma(vi)

          t4ppt_ilma(vi,n) = ilma
          t4ppt_j   (vi,n) = j
       enddo
    enddo

!$acc enter data copyin(t4ppt_nlma,t4ppt_i2vi,t4ppt_vi2i,t4ppt_ilma,t4ppt_j)

    allocate(t4cp_uVpsix(Nlma, NKB))
    allocate(t4cp_uVpsiy(Nlma, NKB))
    allocate(t4cp_uVpsiz(Nlma, NKB))

  end subroutine
#endif

  subroutine auto_blocking
    use misc_routines, only: floor_pow2
    implicit none
    integer,parameter :: L1cache_size =  8 * 1024
    integer,parameter :: value_size   = 24
    real(8) :: nyx
    integer :: sq

    nyx = dble(L1cache_size) / (PNLz * value_size)
    sq  = int(floor(sqrt(nyx)))

    STENCIL_BLOCKING_X = floor_pow2(min(sq, PNLx))
    STENCIL_BLOCKING_Y = floor_pow2(min(sq, PNLy))
  end subroutine

  subroutine opt_vars_reinitialize
    implicit none

    call opt_vars_finalize
    call opt_vars_initialize_p1
    call opt_vars_initialize_p2
  end subroutine

  subroutine opt_vars_finalize
    implicit none

#define SAFE_DEALLOCATE(var) if(allocated(var)) deallocate(var)

    SAFE_DEALLOCATE(zhtpsi)
    SAFE_DEALLOCATE(zttpsi)

    SAFE_DEALLOCATE(zrhotmp)
    SAFE_DEALLOCATE(zJxyz)
    SAFE_DEALLOCATE(zKxyz)

    SAFE_DEALLOCATE(modx)
    SAFE_DEALLOCATE(mody)
    SAFE_DEALLOCATE(modz)

    SAFE_DEALLOCATE(zcx)
    SAFE_DEALLOCATE(zcy)
    SAFE_DEALLOCATE(zcz)

    SAFE_DEALLOCATE(nprojector)
    SAFE_DEALLOCATE(idx_proj)
    SAFE_DEALLOCATE(idx_lma)
    SAFE_DEALLOCATE(pseudo_start_idx)

#ifdef ARTED_STENCIL_ORIGIN
    SAFE_DEALLOCATE(zifdx)
    SAFE_DEALLOCATE(zifdy)
    SAFE_DEALLOCATE(zifdz)
#endif

#ifdef ARTED_LBLK
    SAFE_DEALLOCATE(t4ppt_nlma)
    SAFE_DEALLOCATE(t4ppt_i2vi)
    SAFE_DEALLOCATE(t4ppt_vi2i)
    SAFE_DEALLOCATE(t4ppt_ilma)
    SAFE_DEALLOCATE(t4ppt_j)

    SAFE_DEALLOCATE(t4cp_uVpsix)
    SAFE_DEALLOCATE(t4cp_uVpsiy)
    SAFE_DEALLOCATE(t4cp_uVpsiz)
#endif
  end subroutine
end module opt_variables
