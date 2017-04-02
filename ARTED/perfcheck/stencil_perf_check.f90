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
subroutine stencil_perf_check(NLx_, NLy_, NLz_, NK_, NB_, Nt_)
  use global_variables, only: NLx,NLy,NLz,NK,NBoccmax,NUMBER_THREADS,Nt
  use wrap_variables
  use performance_analyzer
  use timer
  use omp_lib
  implicit none
  integer,intent(in) :: NLx_, NLy_, NLz_
  integer,intent(in) :: NK_, NB_, Nt_

  integer :: t
  real(8) :: tbeg,tend
  real(8) :: time(4), gflops(4)

  NUMBER_THREADS = omp_get_max_threads()

  NLx      = NLx_
  NLy      = NLy_
  NLz      = NLz_
  NK       = NK_
  NBoccmax = NB_
  Nt       = Nt_

  print '("THREADS  =",i6)', NUMBER_THREADS
  print '("NK       =",i6)', NK
  print '("NBoccmax =",i6)', NBoccmax
  print '("NLx,y,z  =",i6,i6,i6)', NLx,NLy,NLz
  print '("step     =",i6)', Nt
  print *, ''
  Nt = Nt - 1

  call wrap_init
  call timer_initialize

  call dt_evolve_hpsi

  call timer_reset
  tbeg = omp_get_wtime()
  do t=0,Nt
    call dt_evolve_hpsi
  end do
  tend = omp_get_wtime()

  time(1) = timer_get(LOG_HPSI_STENCIL)
  time(2) = timer_get(LOG_HPSI_UPDATE)
  call get_hamiltonian_performance(gflops)

  print '("total   time =",f6.2," sec")', tend - tbeg
  print '("stencil time =",f6.2," sec,",f8.2," GFLOPS")', time(1), gflops(1)
  print '("update  time =",f6.2," sec,",f8.2," GFLOPS")', time(2), gflops(3)
end subroutine

subroutine err_finalize(err_message)
  use Global_Variables
  implicit none
  character(*),intent(in) :: err_message
  write(*,*) err_message
  stop
end subroutine err_finalize
