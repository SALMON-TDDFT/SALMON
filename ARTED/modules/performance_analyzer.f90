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
module performance_analyzer
  implicit none

  public  write_performance
  public  print_stencil_size
  public  get_hamiltonian_performance

  private write_hamiltonian
  private write_loadbalance

  private summation_threads, get_gflops, get_hamiltonian_chunk_size
  private get_stencil_FLOP, get_pseudo_pt_FLOP, get_update_FLOP

contains
  subroutine print_stencil_size
    use global_variables, only: NK_s,NK_e,NBoccmax,NL,Nt,NUMBER_THREADS
    implicit none
    integer :: NK, NB

    NK = NK_e - NK_s + 1
    NB = NBoccmax
    print *, 'NK =', NK
    print *, 'NB =', NB
    print *, 'NL =', NL
    print *, 'Nt =', (Nt + 1)
    print *, 'Number of Domain/Thread =', real(NK * NB) / NUMBER_THREADS
  end subroutine

  subroutine get_hamiltonian_performance(lgflops)
    implicit none
    real(8) :: lgflops(4)
    call summation_threads(lgflops)
  end subroutine

  subroutine write_hamiltonian(iounit)
    use global_variables
    use salmon_parallel
    use salmon_communication
    use timer
    implicit none
    integer,intent(in) :: iounit

    character(*),parameter :: f = '(A,4(f15.6))'

    type(comm_maxloc_type) :: tin,tout
    real(8)                :: lgflops(4), pgflops(4), tgflops(4)
    real(8)                :: sgflops(4)

    call summation_threads(lgflops)
    pgflops = lgflops

    tin%rank = nproc_id_global
    tin%val  = lgflops(4)
    call comm_get_max(tin, tout, nproc_group_global)
    call comm_bcast(pgflops, nproc_group_global, tout%rank)

    if (calc_mode == calc_mode_ms) then
       call comm_summation(lgflops, sgflops, 4, nproc_group_tdks)
    end if

    call comm_summation(lgflops, tgflops, 4, nproc_group_global)

    if(comm_is_root(nproc_id_global)) then
      write (iounit,'(A)') 'Performance [GFLOPS]'
      write (iounit,'(A,4(A15))') 'Type           ', 'Hamiltonian', 'Stencil', 'Pseudo-Pt', 'Update'
      write (iounit,f)            'Processor      ', lgflops(4), lgflops(1), lgflops(2), lgflops(3)
      write (iounit,f)            'Processor(max) ', pgflops(4), pgflops(1), pgflops(2), pgflops(3)
      if (calc_mode == calc_mode_ms) then
        write (iounit,f)            'Macro-grid(sum)', sgflops(4), sgflops(1), sgflops(2), sgflops(3)
      endif
      write (iounit,f)            'System(sum)    ', tgflops(4), tgflops(1), tgflops(2), tgflops(3)
    end if
  end subroutine

  subroutine write_loadbalance(iounit)
    use global_variables
    use salmon_parallel
    use salmon_communication
    use timer
    implicit none
    integer,intent(in) :: iounit

    integer,parameter      :: LOG_SIZE=12
    character(*),parameter :: f = "(A,3(F12.4),F12.2)"

    real(8) :: src(LOG_SIZE), rmin(LOG_SIZE), rmax(LOG_SIZE), diff(LOG_SIZE), rel(LOG_SIZE)

    src( 1) = timer_get(LOG_DT_EVOLVE)
    src( 2) = timer_get(LOG_HPSI)
    src( 3) = timer_get(LOG_PSI_RHO)
    src( 4) = timer_get(LOG_HARTREE)
    src( 5) = timer_get(LOG_CURRENT)
    src( 6) = timer_get(LOG_TOTAL_ENERGY)
    src( 7) = timer_get(LOG_ION_FORCE)
    src( 8) = timer_get(LOG_DT_EVOLVE_AC)
    src( 9) = timer_get(LOG_K_SHIFT_WF)
    src(10) = timer_get(LOG_OTHER)
    src(11) = timer_get(LOG_ALLREDUCE)
    src(12) = timer_get(LOG_DYNAMICS)

    call comm_get_min(src,rmin,LOG_SIZE,nproc_group_global)
    call comm_get_max(src,rmax,LOG_SIZE,nproc_group_global)

    diff(:) = rmax(:) - rmin(:)
    rel(:)  = rmax(:) / rmin(:)

    if (comm_is_root(nproc_id_global)) then
      write (iounit,'(A)') 'Load balance check [sec]'
      write (iounit,'(A,4(A12))') 'Function    ','min','max','diff','rel'
      write (iounit,f)            'dt_evolve   ',rmin( 1),rmax( 1),diff( 1),rel( 1)
      write (iounit,f)            'hpsi        ',rmin( 2),rmax( 2),diff( 2),rel( 2)
      write (iounit,f)            'psi_rho     ',rmin( 3),rmax( 3),diff( 3),rel( 3)
      write (iounit,f)            'hartree     ',rmin( 4),rmax( 4),diff( 4),rel( 4)
      write (iounit,f)            'current     ',rmin( 5),rmax( 5),diff( 5),rel( 5)
      write (iounit,f)            'total_energy',rmin( 6),rmax( 6),diff( 6),rel( 6)
      write (iounit,f)            'ion_force   ',rmin( 7),rmax( 7),diff( 7),rel( 7)
      write (iounit,f)            'dt_evolve_ac',rmin( 8),rmax( 8),diff( 8),rel( 8)
      write (iounit,f)            'k_shift_wf  ',rmin( 9),rmax( 9),diff( 9),rel( 9)
      write (iounit,f)            'other       ',rmin(10),rmax(10),diff(10),rel(10)
      write (iounit,f)            'allreduce   ',rmin(11),rmax(11),diff(11),rel(11)
      write (iounit,f)            'dynamics    ',rmin(12),rmax(12),diff(12),rel(12)
    end if
  end subroutine

  subroutine write_performance(filename)
    use global_variables
    use salmon_parallel
    use salmon_communication
    use misc_routines, only: gen_logfilename
    implicit none
    character(*),intent(in) :: filename

    integer,parameter :: iounit = 999

    if(comm_is_root(nproc_id_global)) open(iounit, file=gen_logfilename(filename))
    call write_hamiltonian(iounit)
    if(comm_is_root(nproc_id_global)) write (iounit,'(A)') '==='
    call write_loadbalance(iounit)
    if(comm_is_root(nproc_id_global)) close(iounit)

    call comm_sync_all
  end subroutine

  subroutine summation_threads(lgflops)
    use global_variables, only: NUMBER_THREADS, functional, propagator
    use timer
    implicit none
    real(8), intent(out) :: lgflops(4)
    real(8) :: hflop(3), htime(4)
    integer :: i, cnt
    integer :: chunk_size(0:NUMBER_THREADS-1)
    integer :: ncalls_in_loop

    select case(functional)
      case('VS98','TPSS','TBmBJ')
        ncalls_in_loop = 3
      case default
        ncalls_in_loop = 2
    end select

    select case(propagator)
      case('middlepoint')
        ncalls_in_loop = ncalls_in_loop - 1
    end select

    call get_hamiltonian_chunk_size(chunk_size)

    lgflops = 0.0d0
#ifndef _OPENACC
    do i=0,NUMBER_THREADS-1
#else
    do i=0,0
#endif
      cnt = chunk_size(i)

      hflop(1) = get_stencil_FLOP(cnt)   * ncalls_in_loop
      hflop(2) = get_pseudo_pt_FLOP(cnt) * ncalls_in_loop
      hflop(3) = get_update_FLOP(cnt)    * ncalls_in_loop

      htime(1) = timer_thread_get(LOG_HPSI_STENCIL, i)
      htime(2) = timer_thread_get(LOG_HPSI_PSEUDO, i)
      htime(3) = timer_thread_get(LOG_HPSI_UPDATE, i)
      htime(4) = timer_thread_get(LOG_HPSI_INIT, i)

      lgflops(1) = lgflops(1) + get_gflops(hflop(1), htime(1))
      lgflops(2) = lgflops(2) + get_gflops(hflop(2), htime(2))
      lgflops(3) = lgflops(3) + get_gflops(hflop(3), htime(3))
      lgflops(4) = lgflops(4) + get_gflops(sum(hflop), sum(htime))
    end do
  end subroutine

  subroutine get_hamiltonian_chunk_size(chunk_size)
    use global_variables, only: NKB,NUMBER_THREADS
    use omp_lib, only: omp_get_thread_num
    implicit none
    integer,intent(out) :: chunk_size(0:NUMBER_THREADS-1)
    integer :: ikb, tid

#ifndef _OPENACC
    tid = 0
    chunk_size(:) = 0
!$omp parallel private(tid) shared(chunk_size)
!$ tid = omp_get_thread_num()
!$omp do private(ikb)
    do ikb=1,NKB
      chunk_size(tid) = chunk_size(tid) + 1
    end do
!$omp end do
!$omp end parallel
#else
    chunk_size(:) = 0
    chunk_size(0) = NKB
#endif
  end subroutine

  function get_stencil_FLOP(chunk_size)
    use global_variables, only: NK_s,NK_e,NXY_s,NXY_e,NBoccmax,NL,Nt
    integer,intent(in),optional :: chunk_size
    real(8),parameter           :: FLOP = 158

    real(8) :: get_stencil_FLOP
    integer :: nsize

    if(present(chunk_size)) then
      nsize = chunk_size &
            * (NXY_e - NXY_s + 1)
    else
      nsize = (NK_e  - NK_s  + 1) * NBoccmax &
            * (NXY_e - NXY_s + 1)
    end if
    get_stencil_FLOP = nsize * 4*FLOP*NL * (Nt + 1)
  end function

  function get_pseudo_pt_FLOP(chunk_size)
    use global_variables, only: NK_s,NK_e,NXY_s,NXY_e,NBoccmax,Nt,a_tbl,Mps
    implicit none
    integer,intent(in),optional :: chunk_size
    real(8),parameter           :: FLOP_reduction = (2 + 6)     + 2
    real(8),parameter           :: FLOP_scatter   = (2 + 6 + 1) + 2 ! 1 = conjg(z)
    real(8),parameter           :: FLOP_scalar    = 2 + 2

    real(8) :: get_pseudo_pt_FLOP
    real(8) :: FLOP
    integer :: nsize

    FLOP = FLOP_scalar + (FLOP_reduction + FLOP_scatter) * sum(Mps(a_tbl(:)))

    if(present(chunk_size)) then
      nsize = chunk_size &
            * (NXY_e - NXY_s + 1)
    else
      nsize = (NK_e  - NK_s  + 1) * NBoccmax &
            * (NXY_e - NXY_s + 1)
    endif
    get_pseudo_pt_FLOP = nsize * 4*FLOP * (Nt + 1)
  end function

  function get_update_FLOP(chunk_size)
    use global_variables, only: NK_s,NK_e,NXY_s,NXY_e,NBoccmax,NL,Nt
    implicit none
    integer,intent(in),optional :: chunk_size
    real(8),parameter           :: FLOP = 6 + 2

    real(8) :: get_update_FLOP
    integer :: nsize

    if(present(chunk_size)) then
      nsize = chunk_size &
            * (NXY_e - NXY_s + 1)
    else
      nsize = (NK_e  - NK_s  + 1) * NBoccmax &
            * (NXY_e - NXY_s + 1)
    endif
    get_update_FLOP = nsize * 4*FLOP*NL * (Nt + 1)
  end function

  function get_gflops(FLOP,time)
    implicit none
    real(8),intent(in) :: FLOP
    real(8),intent(in) :: time
    real(8) :: get_gflops
    get_gflops = FLOP / (time * (10**9))
  end function
end module
