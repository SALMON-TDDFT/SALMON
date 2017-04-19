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
module opt_variables
  implicit none

  real(8) :: lapt(12)

  integer                :: PNLx,PNLy,PNLz,PNL
  complex(8),allocatable :: zhtpsi(:,:,:),zttpsi(:,:)

  real(8),allocatable :: zrhotmp(:,:)

  integer,allocatable :: zJxyz(:,:),zKxyz(:,:)

  integer,allocatable :: modx(:),mody(:),modz(:)

  integer :: STENCIL_BLOCKING_X
  integer :: STENCIL_BLOCKING_Y

  real(8),allocatable :: zcx(:,:),zcy(:,:),zcz(:,:)

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

#ifdef ARTED_REDUCE_FOR_MANYCORE
    tid_range = ceiling_pow2(NUMBER_THREADS) - 1
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

#ifndef ARTED_LBLK
    allocate(zhtpsi(0:PNL-1,4,0:NUMBER_THREADS-1))
#else
    blk_nkb_hpsi = min(at_least_parallelism/PNL + 1, NKB)
    allocate(zhtpsi(0:PNL-1, 0:blk_nkb_hpsi-1, 4))
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

    allocate(zJxyz(Nps,NI))

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

  subroutine symmetric_load_balancing(NK,NK_ave,NK_s,NK_e,NK_remainder,procid,nprocs)
    use environment
    implicit none
    integer,intent(in)    :: NK
    integer,intent(in)    :: NK_ave
    integer,intent(inout) :: NK_s
    integer,intent(inout) :: NK_e
    integer,intent(inout) :: NK_remainder
    integer,intent(in)    :: procid
    integer,intent(in)    :: nprocs

    integer :: NScpu,NSmic,NPcpu,NPmic,NPtotal
    integer :: np,npr,pos

    NPcpu   = CPU_PROCESS_PER_NODE
    NPmic   = MIC_PROCESS_PER_NODE
    NPtotal = NPcpu + NPmic

    if (procid == 0 .and. CPU_PROCESS_PER_NODE /= MIC_PROCESS_PER_NODE) then
      call err_finalize('CPU_PROCESS_PER_NODE /= MIC_PROCESS_PER_NODE')
    end if

    ! NK
    NScpu = int(NK_ave * CPU_TASK_RATIO)
    NSmic = int(NK_ave * MIC_TASK_RATIO)

    NK_remainder = NK - (NScpu * (Nprocs/2) + NSmic * (Nprocs/2))

    np  = procid / NPtotal * NPtotal
    npr = mod(procid, NPtotal)
    pos = (np / 2) * NScpu &
    &   + (np / 2) * NSmic

    if (npr < NPcpu) then
      pos = pos + npr * NScpu
    else
      pos = pos + NPcpu          * NScpu
      pos = pos + mod(npr,NPmic) * NSmic
    end if
    NK_s = pos
#ifdef __MIC__
    NK_e = pos + NSmic
#else
    NK_e = pos + NScpu
#endif
    if (procid+1 == nprocs .and. NK_remainder /= 0) then
      NK_e = NK_e + NK_remainder
    end if
    NK_s = NK_s + 1

    ! Error check
    if(procid == nprocs-1 .and. NK_e /= NK) then
      call err_finalize('prep. NK_e error')
    end if
  end subroutine

  function is_symmetric_mode()
    use global_variables
    use communication
    implicit none
    integer :: is_symmetric_mode
    logical :: arch, ret

#ifdef __MIC__
    arch = .TRUE.
#else
    arch = .FALSE.
#endif

    call comm_logical_and(arch, ret, proc_group(1))

    if(ret) then
      is_symmetric_mode = 1
    else
      is_symmetric_mode = 0
    end if
  end function
end module opt_variables
