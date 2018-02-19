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
!This file is "optimization.f90"
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

!-----------------------------------------------------------------
! Optimization in the ground state by conjugate gradient method
!-----------------------------------------------------------------
module optimization
  implicit none
contains

subroutine calc_opt_ground_state
  use Global_Variables
  use timer
  use opt_variables
  use performance_analyzer
  use salmon_parallel, only: nproc_group_global, nproc_id_global
  use salmon_communication, only: comm_bcast, comm_sync_all, comm_is_root
  use misc_routines, only: get_wtime
  use salmon_global, only: format3d, out_dns, out_dns_rt, out_dns_rt_step
  use ground_state
  use io_gs_wfn_k
  use md_ground_state, only: write_xyz
  implicit none
  integer :: ia, nfree

  !check
  nfree=0
  do ia=1,NI
     if(flag_geo_opt_atom(ia)=='y') nfree=nfree+1
  enddo
  if(nfree==0)then
     if(comm_is_root(nproc_id_global)) &
     &   write(*,*) "ERROR: No free atoms in optimization. Check input file"
     call end_parallel
     stop
  endif


  if(flag_use_grad_wf_on_force) then
     call calc_opt_ground_state_useF
  else
     call calc_opt_ground_state_useE
  endif
  return

contains

!----------------------------------------------
subroutine calc_opt_ground_state_useF

  implicit none
  integer :: i,j,k,iter_perp,iter_line, Nopt_perp,Nopt_line,ia
  real(8) :: gm, tmp1,tmp2,tmp3, Rion_save(3,NI),dRion_rmsd
  real(8) :: StepLen_line0, StepLen_line(2), StepLen_line_small
  real(8) :: StepLen_line_new, StepLen_line_zero, sl_zero
  real(8) :: StepLen_line_Up, StepLen_line_Dw
  real(8) :: SearchDirection(3,NI),SearchDirection_1d(3*NI)
  real(8) :: Eall_prev,Eall_prev_line,Eall_zero,Eall_new,dEall
  real(8) :: force_prev(3,NI), force_1d(3*NI), force_prev_1d(3*NI)
  real(8) :: fmax,fave, fmax_conv, F_line_conv 
  real(8) :: F_line_2pt(2), F_line,F_line_new,F_line_zero
  logical :: flag_accept
  character(100) :: comment_line

  position_option='rewind'

  if(comm_is_root(nproc_id_global)) then
     write(*,*) "===== Grand State Optimization Start ====="
     write(*,*) "       (CG method using Force only)       "
  endif
  PrLv_scf = 0
  flag_scf_conv_ene_force=.true.
 !if(convrg_scf_ene < 0d0) then
 !   convrg_scf_ene = convrg_opt_fmax*1d-2
 !   if(convrg_scf_ene.ge.1d-8) convrg_scf_ene=1d-8
 !endif
  iflag_gs_init_wf=1   !flag to skip giving randam number initial guess
  Nopt_perp = 100
  Nopt_line = 40
  StepLen_line_small = 0.1d0   !(small step to guess initial step length)
 !StepLen_line0      = cg_alpha_ini  !(not use now)
  StepLen_line_Up    = cg_alpha_up
  StepLen_line_Dw    = cg_alpha_down
  F_line_conv        = convrg_opt_fmax
  fmax_conv          = convrg_opt_fmax

  if(comm_is_root(nproc_id_global)) then
     write(*,*) "  [Set following in optimization]"
     write(*,*) "  SCF convergence threshold(E)=",real(convrg_scf_ene)
     write(*,*) "  SCF convergence threshold(F)=",real(convrg_scf_force)
     write(*,*) "  Max optimization CG step    =",Nopt_perp
     write(*,*) "  Max line search opt step    =",Nopt_line
    !write(*,*) "  Ini. line-search step length param. =",real(StepLen_line0)
     write(*,*) "  Up rate of line-search step length  =",real(StepLen_line_Up)
     write(*,*) "  Down rate of line-search step length=",real(StepLen_line_Dw)
     write(*,*) "  Convergence threshold of F_line =",real(F_line_conv)
     write(*,*) "  Convergence threshold of Fmax   =",real(fmax_conv)
  endif

  call comm_sync_all

  !Initial Step Procedure
  call Ion_Force_omp(rion_update_on,calc_mode_gs)
  SearchDirection(:,:) = force(:,:)
  do ia=1,NI
     if(flag_geo_opt_atom(ia)=='n') SearchDirection(:,ia)=0d0  !fix atom
  enddo
  call variable_3xNto3N(NI,force,force_1d)
  call variable_3xNto3N(NI,SearchDirection,SearchDirection_1d)
  call cal_inner_product(3*NI,force_1d,SearchDirection_1d,F_line)

  write(comment_line,110) 0, 0
  call write_xyz(comment_line,"new","r  ")
  call cal_mean_max_forces(NI,force,fave,fmax)
  if(comm_is_root(nproc_id_global)) write(*,135) 0,0,fmax,fave,Eall,0d0

  !(initial check of convergence: if force is enough small, exit before opt)
  if(fmax .le. fmax_conv) then
     if(comm_is_root(nproc_id_global)) &
     &  write(*,*) " Max force is enough small: stop calculation"
     call end_parallel
     stop
  endif

  !==== Main Loop ====

  !---iteration for perpendicular direction---
  do iter_perp =1,Nopt_perp   

    if(comm_is_root(nproc_id_global)) then
       write(*,*) "==================================================="
       write(*,*) "CG Optimization Step = ", iter_perp
    endif

    !previous value
    force_prev(:,:) = force(:,:)
    Eall_prev = Eall
    Eall_prev_line = Eall

    !(store)
    Rion_save(:,:)= Rion(:,:)
    call manipulate_wfn_data('save')

    !Set initial region to be searched (=> set initial two points)
    !(calculate forces at 2 points and adjust initial step length)
    call cal_inner_product(3*NI,force_1d,SearchDirection_1d,F_line)
    StepLen_line(1) = 0d0
    F_line_2pt(1)   = F_line

    if(F_line_2pt(1).le.0d0) then
       if(comm_is_root(nproc_id_global)) &
       &  write(*,*) "WARNING: Search Direction is opposite"
       exit
    endif

    !(guess initial step length from numerical gradient of force)
    i= 2
    StepLen_line(i) = StepLen_line_small
    do
      call manipulate_wfn_data('load')
      call read_write_gs_wfn_k(iflag_read)
      Rion(:,:)   = Rion_save(:,:) + StepLen_line(i)* SearchDirection(:,:)
      Rion_eq(:,:)= Rion(:,:)
      call prep_ps_periodic('update_all       ')
      call calc_ground_state
      call Ion_Force_omp(rion_update_on,calc_mode_gs)
      call variable_3xNto3N(NI,force,force_1d)
      call cal_inner_product(3*NI,force_1d,SearchDirection_1d,F_line)
      F_line_2pt(i) = F_line
      call get_predicted_zero(StepLen_line,F_line_2pt,sl_zero)
      dRion_rmsd  = sl_zero*sqrt(sum(SearchDirection_1d(:)**2)/NI)
      if(comm_is_root(nproc_id_global)) &
      &  write(*,104) " guessed-initial-alpha=",sl_zero,"  dRion-RMSD=",dRion_rmsd

      if(abs(dRion_rmsd).lt.1d-6 .or. dRion_rmsd.gt.1d0) then
         ! just ad-hoc treatment (unusual case: I don't know why this case exists)
         StepLen_line(i) = StepLen_line_small*5d0
         if(comm_is_root(nproc_id_global)) then
            write(*,*) " Could not find good guess of initial alpha: too small or large RMSD"
            write(*,*) " --> put guessed initial alpha=",real(StepLen_line(i))
         endif
         exit
      endif

      if(sl_zero.lt.0d0) then
         StepLen_line(i) = StepLen_line(i) * StepLen_line_Dw
      else 
         StepLen_line(i) = sl_zero
         exit
      endif
    enddo

    !(adjust initial step length)
    call manipulate_wfn_data('load')
10  continue
    i= 2
   !call manipulate_wfn_data('load')
    call read_write_gs_wfn_k(iflag_read)
    Rion(:,:)   = Rion_save(:,:) + StepLen_line(i)* SearchDirection(:,:)
    Rion_eq(:,:)= Rion(:,:)
    dRion_rmsd  = StepLen_line(i)*sqrt(sum(SearchDirection_1d(:)**2)/NI)
    call prep_ps_periodic('update_all       ')
    call calc_ground_state
    call Ion_Force_omp(rion_update_on,calc_mode_gs)
    call variable_3xNto3N(NI,force,force_1d)
    call cal_inner_product(3*NI,force_1d,SearchDirection_1d,F_line)
    F_line_2pt(i) = F_line

    !(adjust initial step length so as to include zero-force between 2 points)
    if( F_line_2pt(1)*F_line_2pt(2).gt.0d0 ) then
    !Here, this algorithm can be replaced by linearly predicted point(step length)
    !if a lot of "up", modify here with prediction
       StepLen_line(1) = StepLen_line(2)
       F_line_2pt(1)   = F_line_2pt(2)
       StepLen_line(2) = StepLen_line(2) * StepLen_line_Up
       dRion_rmsd      = StepLen_line(2)*sqrt(sum(SearchDirection_1d(:)**2)/NI)
       if(comm_is_root(nproc_id_global)) &
       & write(*,104) " adjusting initial alpha(up)->",StepLen_line(2),"  dRion-RMSD=",dRion_rmsd
       goto 10
    else
       if(comm_is_root(nproc_id_global)) &
       & write(*,104) " initial alpha=",StepLen_line(2),"  dRion-RMSD=",dRion_rmsd
    endif


    !--- iteration loop for line search to find zero-force ---
    do iter_line= 1,Nopt_line

       !(narrow search range)
       call get_predicted_zero(StepLen_line,F_line_2pt,sl_zero)

       if(abs(F_line_2pt(1)/F_line_2pt(2)).gt.1d0) then
          StepLen_line_new= 0.5d0*( StepLen_line(1) + sl_zero )
       else
          StepLen_line_new= StepLen_line(1) + 0.5d0*(StepLen_line(2)-sl_zero)
       endif

       !(calculate electronic state at the new point)
       call manipulate_wfn_data('load')
       call read_write_gs_wfn_k(iflag_read)
       Rion(:,:)= Rion_save(:,:) + StepLen_line_new* SearchDirection(:,:)
       Rion_eq(:,:)= Rion(:,:)
       write(comment_line,110) iter_perp, iter_line
       call write_xyz(comment_line,"add","r  ")
       call prep_ps_periodic('update_all       ')
       call calc_ground_state
       call Ion_Force_omp(rion_update_on,calc_mode_gs)
       call variable_3xNto3N(NI,force,force_1d)
       call cal_inner_product(3*NI,force_1d,SearchDirection_1d,F_line)
       F_line_new = F_line
       Eall_new   = Eall

       flag_accept=.false.
       if(abs(F_line_2pt(1)/F_line_2pt(2)).gt.1d0) then
          if( F_line_new*F_line_2pt(1) .gt. 0d0 ) flag_accept=.true.
       else
          if( F_line_new*F_line_2pt(2) .gt. 0d0 ) flag_accept=.true.
       endif

       if(flag_accept) &
       &  call update_2pt_line_search(StepLen_line,F_line_2pt,StepLen_line_new,F_line_new)


       !(get predicted zero-crossin coordinate)
       call get_predicted_zero(StepLen_line,F_line_2pt,StepLen_line_zero)

       !(calculate electronic state at the predicted zero)
       call manipulate_wfn_data('load')
       call read_write_gs_wfn_k(iflag_read)
       Rion(:,:)   = Rion_save(:,:) + StepLen_line_zero* SearchDirection(:,:)
       Rion_eq(:,:)= Rion(:,:)
       dRion_rmsd  = StepLen_line_zero*sqrt(sum(SearchDirection_1d(:)**2)/NI)
       write(comment_line,110) iter_perp, iter_line
       call write_xyz(comment_line,"add","r  ")
       call prep_ps_periodic('update_all       ')
       call calc_ground_state
       call Ion_Force_omp(rion_update_on,calc_mode_gs)
       call variable_3xNto3N(NI,force,force_1d)
       call cal_inner_product(3*NI,force_1d,SearchDirection_1d,F_line)
       F_line_zero = F_line
       Eall_zero   = Eall

       !(log)
       if(comm_is_root(nproc_id_global)) then
          write(*,100) "   alpha-line:2pt=",(StepLen_line(i),i=1,2),&
                     & " |predict-zero=",    StepLen_line_zero
          write(*,100) "   force-line:2pt=",(F_line_2pt(i),i=1,2), &
                     & " |predict-zero=",    F_line_zero
       endif

       !Judge Convergence for line search opt
       dEall   = Eall_zero - Eall_prev_line
       if(comm_is_root(nproc_id_global)) &
       & write(*,130)iter_perp,iter_line,StepLen_line_zero,F_line_zero,dEall,dRion_RMSD

       if(abs(F_line_zero) .le. F_line_conv)then
          if(comm_is_root(nproc_id_global)) then
             write(*,*)
             write(*,*) "Converged(line search opt) in perpendicular-step",iter_perp
          endif
          call read_write_gs_wfn_k(iflag_write)
          call read_write_gs_wfn_k(iflag_read)
          call manipulate_wfn_data('save')
          exit
       else if(iter_line==Nopt_line) then
          if(comm_is_root(nproc_id_global)) then
             write(*,*) "Not Converged(line search opt) in perpendicular-step",iter_perp
             write(*,*) "==================================================="
          endif
          exit
       endif

       !(update two points among current predicted-zero and two points)
       call update_2pt_line_search(StepLen_line,F_line_2pt,StepLen_line_zero,F_line_zero)

       !(preparation for next cycle)
       call manipulate_wfn_data('load')
       Eall_prev_line = Eall_zero

    enddo

    !Force calculation
    call Ion_Force_omp(rion_update_on,calc_mode_gs)
    call variable_3xNto3N(NI,force,force_1d)
    call cal_mean_max_forces(NI,force,fave,fmax)

    !Judge Convergence for perpendicular opt (main judge)
    dEall = Eall - Eall_prev
    if(comm_is_root(nproc_id_global)) &
    &  write(*,135) iter_perp, iter_line, fmax,fave,Eall,dEall
    if(abs(Fmax) .le. fmax_conv) then
       if(comm_is_root(nproc_id_global)) then
          write(*,*) "==================================================="
          write(*,*)
          write(*,*) "Optimization Converged"
          write(*,150) iter_perp,fmax,fave,Eall
          write(*,*) "==================================================="
       endif
       exit
    else if(iter_perp==Nopt_perp) then
       if(comm_is_root(nproc_id_global)) then
          write(*,*) "==================================================="
          write(*,*)
          write(*,*) "Optimization Did Not Converged"
          write(*,150) iter_perp,fmax,fave,Eall
          write(*,*) "==================================================="
       endif
    endif


    !Update search direction vector for perpendicular step
    call variable_3xNto3N(NI,force,force_1d)
    call variable_3xNto3N(NI,force_prev,force_prev_1d)
    call cal_inner_product(3*NI,force_1d,force_1d,tmp1)
    call cal_inner_product(3*NI,force_prev_1d,force_prev_1d,tmp2)
    call cal_inner_product(3*NI,force_1d,force_prev_1d,tmp3)
    !gm = tmp1/tmp2         !(by Fletcher-Reeves)
    gm = (tmp1-tmp3)/tmp2   !(by Polak-Ribiere)--usually best, but sometimes direction is not downward.
    SearchDirection(:,:) = force(:,:) + gm * SearchDirection(:,:)
    do ia=1,NI
       if(flag_geo_opt_atom(ia)=='n') SearchDirection(:,ia)=0d0  !fix atom
    enddo
    call variable_3xNto3N(NI,SearchDirection,SearchDirection_1d)

  enddo !end of opt iteraction========================


100 format(a,2e13.5,a,e13.5)
104 format(a,e13.5,a,e13.5)
110 format("#opt   step-perp=",i4,"   step-line",i4)
120 format(a,e14.6,a,e14.6)
130 format(" step=",i3," -",i3,"  alpha=",e12.4,"  F-line=",e12.4,"  dE-line=",e12.4,"  dRion-RMSD=",e12.4)
135 format(" step=(perp=",i3,", line=",i3,")  Fmax=",e11.4,"  Fave=",e11.4,"  E=",e16.8,"  dE-perp=",e12.4)
150 format("    Iteration step=",i5,    / &
    &      "    Force(maximum)=",e18.10," [a.u.]",/ &
    &      "    Force(average)=",e18.10," [a.u.]",/ &
    &      "    Total Energy  =",e18.10," [a.u.]" )

end subroutine calc_opt_ground_state_useF

!----------------------------------------------
subroutine calc_opt_ground_state_useE

  implicit none
  integer :: i,j,k,iter_perp,iter_line,iter_line2, Nopt_perp,Nopt_line
  real(8) :: fmax,fave, gm, tmp1,tmp2,tmp3, dsl12,dsl23,ratio_too_much
  real(8) :: StepLen_LineSearch0, StepLen_LineSearch(3)
  real(8) :: StepLen_LineSearch_min, StepLen_LineSearch_new
  real(8) :: StepLen_LineSearch_Up,StepLen_LineSearch_Dw
  real(8) :: SearchDirection(3,NI), Rion_save(3,NI)
  real(8) :: Eall_prev,Eall_save,Eall_3points(3),Eall_min,Eall_new,Eall_prev_line
  real(8) :: dEall,dE_conv,dE_conv_LineSearch, fmax_conv  !,fave_conv
  real(8) :: force_prev(3,NI), force_1d(3*NI), force_prev_1d(3*NI)
  integer :: nsave
  real(8) :: alpha_save(9999), ene_save(9999)
  logical :: flag_min_found
  character(100) :: comment_line

  position_option='rewind'

  !(check)
  do i=1,NI
     if(flag_geo_opt_atom(i)/='y')then
        if(comm_is_root(nproc_id_global)) &
        &  write(*,*)'ERROR: flag of geometry opt of all atoms must be y'
        call end_parallel
        stop
     endif
  enddo

  if(comm_is_root(nproc_id_global)) &
  &  write(*,*) "===== Grand State Optimization Start ====="
  PrLv_scf = 0
  flag_scf_conv_ene_force=.true.
  if(convrg_scf_ene < 0d0) convrg_scf_ene=1d-6
  iflag_gs_init_wf=1   !flag to skip giving randam number initial guess
  Nopt_perp = 100
  Nopt_line = 40
  StepLen_LineSearch0   = cg_alpha_ini
  StepLen_LineSearch_Up = cg_alpha_up
  StepLen_LineSearch_Dw = cg_alpha_down
  dE_conv_LineSearch    = convrg_opt_ene
  dE_conv               = convrg_opt_ene
  fmax_conv             = convrg_opt_fmax
  !fave_conv             = 1d-5   !--i dont know good number

  if(comm_is_root(nproc_id_global)) then
     write(*,*) "  [Set following in optimization]"
     write(*,*) "  SCF convergence threshold(E)=",real(convrg_scf_ene)
     write(*,*) "  Max optimization CG step    =",Nopt_perp
     write(*,*) "  Max line search opt step    =",Nopt_line
     write(*,*) "  Ini. step length param. for line search    =",real(StepLen_LineSearch0)
     write(*,*) "  Up rate of step length for line search     =",real(StepLen_LineSearch_Up)
     write(*,*) "  Down rate of step length for line search   =",real(StepLen_LineSearch_Dw)
     write(*,*) "  Convergence threshold of dE: line search   =",real(dE_conv_LineSearch)
     write(*,*) "  Convergence threshold of dE: main search   =",real(dE_conv)
     write(*,*) "  Convergence threshold of Fmax: only initial=",real(fmax_conv)
     !write(*,*) "  Convergence threshold of F: only initial=",real(fave_conv)
  endif


  call comm_sync_all

  !Initial Step Procedure
  SearchDirection(:,:) = force(:,:)
  write(comment_line,110) 0, 0
  call write_xyz(comment_line,"new","r  ")
  call cal_mean_max_forces(NI,force,fave,fmax)
  if(comm_is_root(nproc_id_global)) then
     write(*,135) 0, 0, Eall, 0d0
     write(*,120) " Max-force=", fmax, "  Mean-force=", fave
  endif
  !(initial check of convergence: if force is enough small, exit before opt)
  if(fmax .le. fmax_conv) then
    if(comm_is_root(nproc_id_global)) &
    &  write(*,*) " Max force is enough small: stop calculation"
    call end_parallel
    stop
  endif
  !if(fave .le. fave_conv) then
  !  if(comm_is_root(nproc_id_global)) &
  !  &  write(*,*) " Mean force is enough small: stop calculation"
  !  call end_parallel
  !  stop
  !endif

  !--- Main Loop ----
  do iter_perp =1,Nopt_perp     !iteration for perpendicular direction

    if(comm_is_root(nproc_id_global)) then
       write(*,*) "==================================================="
       write(*,*) "CG Optimization Step = ", iter_perp
    endif

    !previous value
    force_prev(:,:) = force(:,:)
    Eall_prev = Eall

    !Line search minimization along search vector direction
    Eall_prev_line = Eall

    !(store)
    Eall_save = Eall
    Rion_save(:,:)= Rion(:,:)
    call manipulate_wfn_data('save')

    !Set initial region to be searched (=> set initial three points)
    !(calculate energy at 3 points and adjust initial step length)
    StepLen_LineSearch(1) = 0d0
    StepLen_LineSearch(2) = StepLen_LineSearch0 * 0.5d0
    StepLen_LineSearch(3) = StepLen_LineSearch0
    Eall_3points(1) = Eall
    nsave=0
    call add_alpha_save(nsave,alpha_save,ene_save,StepLen_LineSearch(1),Eall_3points(1))

1   continue
    !! stop if StepLen_LineSearch is too small or too large(add later)
    !! write(*,*) "Initial structure is bad"
    do i= 2,3
       call manipulate_wfn_data('load')
       call read_write_gs_wfn_k(iflag_read)
       Rion(:,:)= Rion_save(:,:) + StepLen_LineSearch(i)* SearchDirection(:,:)
       Rion_eq(:,:)= Rion(:,:)
       call prep_ps_periodic('update_all       ')
       call calc_ground_state
       Eall_3points(i) = Eall
       call add_alpha_save(nsave,alpha_save,ene_save,StepLen_LineSearch(i),Eall_3points(i))
    enddo

!(open following if you want to find minimum from saved points)
!    if(comm_is_root(nproc_id_global)) &
!    & write(*,200) " alpha,Eall-1,2,3:",real(StepLen_LineSearch(3)),(dble(Eall_3points(i)),i=1,3)
!200 format(a,4e16.8)
!    call find_min_from_save(nsave,alpha_save,ene_save,StepLen_LineSearch,Eall_3points,flag_min_found)
!    if(flag_min_found) goto 2


    !(adjust initial step length so that middle point has lowest energy)
    if(      Eall_3points(2).gt.Eall_3points(1) ) then
       StepLen_LineSearch(:) = StepLen_LineSearch(:) * StepLen_LineSearch_Dw
       if(comm_is_root(nproc_id_global)) &
       & write(*,*) "alpha: adjusting initial (down)--> ",real(StepLen_LineSearch(3))
       goto 1
    else if( Eall_3points(2).gt.Eall_3points(3) )then
       StepLen_LineSearch(:) = StepLen_LineSearch(:) * StepLen_LineSearch_Up
       if(comm_is_root(nproc_id_global)) &
       & write(*,*) "alpha: adjusting initial ( up )--> ",real(StepLen_LineSearch(3))
       goto 1
    endif

2   continue
    StepLen_LineSearch0 = StepLen_LineSearch(3) !for next initial step length

    iter_line2=0
    do iter_line= 1,Nopt_line      !iteration for lise search minimize

       !(narrow search range if the three points are too unbalanced)
       do
          iter_line2 = iter_line2 +1
          !(check and make new point)
          ratio_too_much = 5d0
          dsl12 = abs(StepLen_LineSearch(1)-StepLen_LineSearch(2))
          dsl23 = abs(StepLen_LineSearch(3)-StepLen_LineSearch(2))
          if( dsl12/dsl23 .gt. ratio_too_much ) then
             StepLen_LineSearch_new= 0.5d0*( StepLen_LineSearch(1) + StepLen_LineSearch(2))
          else if( dsl23/dsl12 .gt. ratio_too_much ) then
             StepLen_LineSearch_new= 0.5d0*( StepLen_LineSearch(2) + StepLen_LineSearch(3))
          else
             exit
          endif

          !(calculate electronic state at the new point)
          call manipulate_wfn_data('load')
          call read_write_gs_wfn_k(iflag_read)
          Rion(:,:)= Rion_save(:,:) + StepLen_LineSearch_new* SearchDirection(:,:)
          Rion_eq(:,:)= Rion(:,:)
          write(comment_line,110) iter_perp, iter_line2
          call write_xyz(comment_line,"add","r  ")
          call prep_ps_periodic('update_all       ')
          call calc_ground_state
          Eall_new = Eall
         !call add_alpha_save(nsave,alpha_save,ene_save,StepLen_LineSearch_new,Eall_new)

          !(update three points: current min and two closest points)
          call update_three_points_for_line_min(StepLen_LineSearch,    Eall_3points,&
                                             &  StepLen_LineSearch_new,Eall_new)
       enddo

          
       !(get minimum coordinate by 3 points interpolation)
       call get_minimum_by_interpolation(StepLen_LineSearch,Eall_3points,StepLen_LineSearch_min)

       !(calculate electronic state at the minimum)
       call manipulate_wfn_data('load')
       call read_write_gs_wfn_k(iflag_read)
       Rion(:,:)= Rion_save(:,:) + StepLen_LineSearch_min* SearchDirection(:,:)
       Rion_eq(:,:)= Rion(:,:)
       write(comment_line,110) iter_perp, iter_line2
       call write_xyz(comment_line,"add","r  ")
       call prep_ps_periodic('update_all       ')
       call calc_ground_state
       Eall_min = Eall
      !call add_alpha_save(nsave,alpha_save,ene_save,StepLen_LineSearch_min,Eall_min)

       !(log)
       if(comm_is_root(nproc_id_global)) then
          write(*,100) " alpha: 3-points=",(StepLen_LineSearch(i),i=1,3),&
                     & " | min(parabola)=", StepLen_LineSearch_min
          write(*,100) " Eall:  3-points=",(Eall_3points(i),i=1,3), &
                     & " | min(parabola)=", Eall_min
       endif

       !Judge Convergence for line search opt
       dEall = Eall_min - Eall_prev_line
       if(comm_is_root(nproc_id_global)) &
       & write(*,130) iter_perp, iter_line2, Eall, dEall

       if(abs(dEall) .le. dE_conv_LineSearch)then
          if(comm_is_root(nproc_id_global)) &
          & write(*,*) "Converged line search optimization in perpendicular-step",iter_perp
          call read_write_gs_wfn_k(iflag_write)
          call read_write_gs_wfn_k(iflag_read)
          call manipulate_wfn_data('save')
          exit
       else if(iter_line==Nopt_line) then
          if(comm_is_root(nproc_id_global)) then
             write(*,*) "Not Converged(line search opt) in perpendicular-step",iter_perp
             write(*,*) "==================================================="
          endif
          exit
       endif

       !(update three points: current min and two closest points)
       call update_three_points_for_line_min(StepLen_LineSearch,    Eall_3points, &
                                          &  StepLen_LineSearch_min,Eall_min)

       !(preparation for next cycle)
       call manipulate_wfn_data('load')
       Eall_prev_line = Eall_min

    enddo

    !Force calculation
    ! Note: force by field is zero now, i.e. E_tot=0 and FionAc=0
    call Ion_Force_omp(rion_update_on,calc_mode_gs)
    !call Total_Energy_omp(rion_update_on,calc_mode_gs)
    call cal_mean_max_forces(NI,force,fave,fmax)
    if(comm_is_root(nproc_id_global)) &
    &  write(*,120) " Max-force=", fmax, "  Mean-force=", fave

    !Judge Convergence for perpendicular opt (main judge)
    dEall = Eall - Eall_prev
    if(comm_is_root(nproc_id_global)) &
    &  write(*,135) iter_perp, iter_line2, Eall, dEall
    if(abs(dEall) .le. dE_conv) then
       if(comm_is_root(nproc_id_global)) then
          write(*,*) "Optimization Converged",iter_perp,real(Eall),real(dEall)
          write(*,*) "==================================================="
       endif
       exit
    else if(iter_perp==Nopt_perp) then
       if(comm_is_root(nproc_id_global)) then
          write(*,*) "Optimization Did Not Converged",iter_perp,real(Eall),real(dEall)
          write(*,*) "==================================================="
       endif
    endif


    !Update search direction vector for perpendicular step
    do i=1,NI
    do j=1,3
       k=3*(i-1)+j
       force_1d(k)      = force(j,i)
       force_prev_1d(k) = force_prev(j,i)
    enddo
    enddo
    call cal_inner_product(3*NI,force_1d,force_1d,tmp1)
    call cal_inner_product(3*NI,force_prev_1d,force_prev_1d,tmp2)
    call cal_inner_product(3*NI,force_1d,force_prev_1d,tmp3)
    !gm = 0d0               !just test
    !gm = tmp1/tmp2         !(by Fletcher-Reeves)
    gm = (tmp1-tmp3)/tmp2   !(by Polak-Ribiere)
    SearchDirection(:,:) = force(:,:) + gm * SearchDirection(:,:)

  enddo !end of opt iteraction========================


100 format(a17,3e16.8,a17,e16.8)
110 format("#opt   step-perp=",i4,"   step-line",i4)
120 format(a,e14.6,a,e14.6)
130 format(" step-perp=",i4,"  step-line=",i4,"  E=",e16.8,"  dE-line=",e16.8)
135 format(" step-perp=",i4,"  step-line=",i4,"  E=",e16.8,"  dE-perp=",e16.8)

end subroutine calc_opt_ground_state_useE
    
  subroutine manipulate_wfn_data(action)
    implicit none
    character(4)   :: action
    character(256) :: gs_wfn_directory, gs_wfn_file, occ_file
    character(256) :: gs_wfn_file_save, occ_file_save
    integer,parameter :: nfile_gs_wfn      =  41
    integer,parameter :: nfile_gs_wfn_save = 141
    integer,parameter :: nfile_occ      =  42
    integer,parameter :: nfile_occ_save = 142
    integer :: ik
    integer :: nproc_group_kpoint_ms
    integer :: nproc_id_kpoint_ms
    integer :: nproc_size_kpoint_ms

    write (gs_wfn_directory,'(A,A)') trim(directory),'/gs_wfn_k/'

    if(comm_is_root(nproc_id_global))then
       occ_file      = trim(gs_wfn_directory)//'occupation'
       occ_file_save = trim(gs_wfn_directory)//'occupation_save'
       open(nfile_occ,     file=trim(occ_file),     form='unformatted')
       open(nfile_occ_save,file=trim(occ_file_save),form='unformatted')
       if(action=='save') then
          read(nfile_occ)occ
          call comm_bcast(occ,nproc_group_global)
          write(nfile_occ_save)occ
       else if(action=='load') then
          read(nfile_occ_save)occ
          call comm_bcast(occ,nproc_group_global)
          write(nfile_occ)occ
       endif
       close(nfile_occ)
       close(nfile_occ_save)
    end if

    do ik=NK_s,NK_e
       write(gs_wfn_file,     '(A,A,I7.7,A)') trim(gs_wfn_directory),'/wfn_gs_k',ik,'.wfn'
       write(gs_wfn_file_save,'(A,A,I7.7,A)') trim(gs_wfn_directory),'/wfn_gs_k',ik,'.wfn_save'
       open(nfile_gs_wfn,     file=trim(gs_wfn_file),     form='unformatted')
       open(nfile_gs_wfn_save,file=trim(gs_wfn_file_save),form='unformatted')
       if(action=='save') then
          read( nfile_gs_wfn)     zu_GS(:,:,ik)
          write(nfile_gs_wfn_save)zu_GS(:,:,ik)
       else if(action=='load') then
          read( nfile_gs_wfn_save)zu_GS(:,:,ik)
          write(nfile_gs_wfn)     zu_GS(:,:,ik)
       endif
       close(nfile_gs_wfn)
       close(nfile_gs_wfn_save)
     end do


    if(action=='save')return

    zu_GS0(:,:,:)=zu_GS(:,:,:)
    zu_t(:,:,:)=zu_GS(:,1:NBoccmax,:)
    Rion_eq=Rion
    !dRion(:,:,-1)=0.d0; dRion(:,:,0)=0.d0
    call psi_rho_GS
    call Hartree
    call Exc_Cor(calc_mode_gs,NBoccmax,zu_t)

    Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
    Vloc_GS(:)=Vloc(:)
    call Total_Energy_omp(rion_update_on,calc_mode_gs)
    call Ion_Force_omp(rion_update_on,calc_mode_gs)
    Eall0=Eall
    Eall_GS0=Eall

  end subroutine

  subroutine get_predicted_zero(sl,fl,sl_zero)
    implicit none
    real(8) :: sl(2),fl(2),sl_zero,dsl,dfldsl,b

    dfldsl  = ( fl(2) - fl(1) ) / ( sl(2) - sl(1) )
    b       = dfldsl * sl(1) - fl(1)
    sl_zero = b / dfldsl

    return
  end subroutine

  subroutine update_2pt_line_search(sl,fl,sl_add,fl_add)
    implicit none
    real(8) :: sl(2),fl(2),sl_add,fl_add
    real(8) :: eps
    eps=1d-15

    !error check
    if(sl_add            .lt.-eps   .or.  &
    &  sl(1)             .lt.-eps   .or.  &
    &  sl_add            .lt. sl(1) .or.  &
    &  sl_add            .gt. sl(2) .or.  &
    &  abs(sl(1)-sl_add) .lt. eps   .or.  &
    &  abs(sl(2)-sl_add) .lt. eps   .or.  &
    &  fl(1)             .lt.-eps   .or.  &
    &  fl(1)*fl(2)       .gt. 0d0      ) then
       if(comm_is_root(nproc_id_global))then
          write(*,*) "Strange: Error1"
          write(*,*) real(sl(1)),real(sl(2)),real(sl_add)
          write(*,*) real(fl(1)),real(fl(2)),real(fl_add)
          stop
       endif
    endif

    if( sl_add.gt.sl(1) .and. sl_add.lt.sl(2) ) then
       if(fl_add.gt.0d0) then
          sl(1) = sl_add
          fl(1) = fl_add
       else
          sl(2) = sl_add
          fl(2) = fl_add
       endif
    endif

    return
  end subroutine

  subroutine get_minimum_by_interpolation(sl,ene,sl_min)
    implicit none
    real(8) :: sl(3),ene(3),sl_min,tmpA,tmpB,tmp21,tmp23,tmp21e,tmp23e

    tmp21  = sl(2)-sl(1) 
    tmp21e = ene(2) - ene(1)
    tmp23  = sl(2)-sl(3)
    tmp23e = ene(2) - ene(3)
    tmpA   = (tmp21**2)*tmp23e - (tmp23**2)*tmp21e
    tmpB   = tmp21*tmp23e - tmp23*tmp21e
    sl_min = sl(2) - 0.5d0 * tmpA / tmpB

    return
  end subroutine

  subroutine update_three_points_for_line_min(sl,ene,sl_add,ene_add)
    implicit none
    integer i,i1,i3,imin
    real(8) :: sl_old(3),ene_old(3), sl_add,ene_add, sl_new(3),ene_new(3)
    real(8) :: sl(3),ene(3),sl4(4),ene4(4), dsl,dsl_min1, dsl_min3, ene_lowest

    sl_old(:)  = sl(:)
    ene_old(:) = ene(:)
 
    !(pick up lowest energy point among four(3points+min) as new middle point)
    sl4(1:3)  = sl_old(1:3)
    sl4(4)    = sl_add
    ene4(1:3) = ene_old(1:3)
    ene4(4)   = ene_add

    ene_lowest=1d99
    do i=1,4
       if(ene4(i).le.ene_lowest)then
          imin=i
          ene_lowest=ene4(imin)
       endif
    enddo
    sl_new(2)  = sl4(imin)
    ene_new(2) = ene4(imin)

    dsl_min1=1d99
    dsl_min3=1d99
    do i=1,4
       if(i==imin) cycle
       dsl = sl4(i)-sl4(imin)
       if(dsl.lt.0d0)then
          if(abs(dsl).le.dsl_min1)then
             dsl_min1 = abs(dsl)
             i1 = i
          endif
       endif
       if(dsl.gt.0d0)then
          if(abs(dsl).le.dsl_min3)then
             dsl_min3 = abs(dsl)
             i3 = i
          endif
       endif
    enddo
    sl_new(1)  = sl4(i1)
    ene_new(1) = ene4(i1)
    sl_new(3)  = sl4(i3)
    ene_new(3) = ene4(i3)

    sl(:)  = sl_new(:)
    ene(:) = ene_new(:)

    return
  end subroutine

  subroutine variable_3xNto3N(n,val_3xn,val_3n)
    implicit none
    integer :: i,j,k,n
    real(8) :: val_3xn(3,n), val_3n(3*n)

    do i=1,n
    do j=1,3
       k=3*(i-1)+j
       val_3n(k) = val_3xN(j,i)
    enddo
    enddo

  end subroutine

  subroutine cal_inner_product(n,vec1,vec2,inner_product)
    implicit none
    integer :: i,n
    real(8) :: vec1(n),vec2(n),inner_product
    inner_product=0d0
    do i=1,n
       inner_product = inner_product + vec1(i)*vec2(i)
    enddo
    return
  end subroutine

  subroutine cal_mean_max_forces(NI,f,fave,fmax)
    implicit none
    integer ia,NI,NIfree
    real(8) :: f(3,NI),fave,fmax,fabs
    fmax   = 0d0
    fave   = 0d0
    NIfree = 0
    do ia=1,NI
       if(flag_geo_opt_atom(ia)=='n') cycle  !fix atom
       NIfree = NIfree + 1
       fabs = f(1,ia)**2 + f(2,ia)**2 + f(3,ia)**2
       fave = fave + fabs
       if(fabs .ge. fmax) fmax = fabs
    enddo
    fmax = sqrt(fmax)
    fave = sqrt(fave/NIfree)
  end subroutine

subroutine add_alpha_save(nsave,alpha_save,ene_save,alpha,ene)
! save alpha and it's Eall for line search
  integer :: i,nsave,n
  real(8) :: alpha_save(9999),ene_save(9999),alpha,ene
  real(8) :: alpha_bk(nsave),ene_bk(nsave)


  if(nsave.ne.0) then
     alpha_bk(1:nsave) = alpha_save(1:nsave)
     ene_bk(1:nsave)   = ene_save(1:nsave)
  endif

  n=nsave
  do i=1,n
     if(abs(alpha-alpha_bk(i)).le.1d-14) return  !not add

     if(i==1 .and. alpha.lt.alpha_bk(i)) then
        alpha_save(i)= alpha
        ene_save(i)  = ene
        alpha_save(i+1:nsave+1) = alpha_bk(i:nsave)
        ene_save(i+1:nsave+1)   = ene_bk(i:nsave)
        nsave=nsave+1
        return
     endif

     if(i==nsave) then
        alpha_save(nsave+1) = alpha
        ene_save(nsave+1)   = ene
        nsave=nsave+1
        return
     endif

     if(alpha.gt.alpha_bk(i) .and. alpha.lt.alpha_bk(i+1)) then
        alpha_save(i+1)= alpha
        ene_save(i+1)  = ene
        alpha_save(i+2:nsave+1) = alpha_bk(i+1:nsave)
        ene_save(i+2:nsave+1)   = ene_bk(i+1:nsave)
        nsave=nsave+1
        return
     endif
  enddo

end subroutine

subroutine find_min_from_save(nsave,alpha_save,ene_save,alpha3p,ene3p,flag_min_found)
! save alpha and it's Eall for line search
  integer    i,nsave,i_min
  real(8) :: alpha_save(9999),ene_save(9999),alpha3p(3),ene3p(3)
  real(8) :: e_min,a_min
  logical :: flag_min_found

  e_min = 1d99
  do i=1,nsave
     if(e_min.gt.ene_save(i)) then
        i_min = i
        e_min = ene_save(i)
        a_min = alpha_save(i)
     endif
  enddo

  flag_min_found=.false.
  if(i_min.ne.1 .and. i_min.ne.nsave) then
     ene3p(i_min-1)   = ene_save(i_min-1)
     ene3p(i_min)     = ene_save(i_min)
     ene3p(i_min+1)   = ene_save(i_min+1)
     alpha3p(i_min-1) = alpha_save(i_min-1)
     alpha3p(i_min)   = alpha_save(i_min)
     alpha3p(i_min+1) = alpha_save(i_min+1)
     flag_min_found=.true.
  endif

end subroutine

end subroutine calc_opt_ground_state

end module optimization
