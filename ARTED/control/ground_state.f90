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
!This file is "initialization.f90"
!This file contains initialization of solid state part
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module ground_state
  implicit none
contains
  subroutine calc_ground_state
    use Global_Variables
    use timer
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_bcast, comm_sync_all, comm_is_root
    implicit none
    integer :: iter, ik, ib, ia
    logical :: Rion_update
    character(10) :: functional_t

    Rion_update = rion_update_on

    allocate(rho_in(1:NL,1:Nscf+1),rho_out(1:NL,1:Nscf+1))
    rho_in(1:NL,1:Nscf+1)=0.d0; rho_out(1:NL,1:Nscf+1)=0.d0
    allocate(Eall_GS(0:Nscf),esp_var_ave(1:Nscf),esp_var_max(1:Nscf),dns_diff(1:Nscf))
    
    call init_wf
    call Gram_Schmidt
    rho=0.d0; Vh=0.d0
    
    call psi_rho_GS !sym
    rho_in(1:NL,1)=rho(1:NL)
    
    call Hartree
    functional_t = functional
    if(functional_t == 'TBmBJ' .or. functional_t == 'BJ_PW') functional = 'PZM'
    call Exc_Cor(calc_mode_gs,NBoccmax,zu_t)
    if(functional_t == 'TBmBJ') functional = 'TBmBJ'
    if(functional_t == 'BJ_PW') functional = 'BJ_PW'

    Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
    call Total_Energy_omp(Rion_update,calc_mode_gs)
    call Ion_Force_omp(Rion_update,calc_mode_gs)
    if (use_ehrenfest_md /= 'y') Rion_update = rion_update_off
    Eall_GS(0)=Eall

    if(comm_is_root(nproc_id_global)) then
       write(*,*) 'This is the end of preparation for ground state calculation'
       call timer_show_current_hour('elapse time=',LOG_ALL)
       write(*,*) '-----------------------------------------------------------'
    end if

    call reset_gs_timer
    call timer_begin(LOG_GROUND_STATE)
    do iter=1,Nscf
       if (comm_is_root(nproc_id_global))  write(*,*) 'iter = ',iter
       if( kbTev < 0d0 )then ! sato
          if (FSset_option == 'Y') then
             if (iter/NFSset_every*NFSset_every == iter .and. iter >= NFSset_start) then
                do ik=1,NK 
                   esp_vb_min(ik)=minval(esp(1:NBocc(ik),ik))
                   esp_vb_max(ik)=maxval(esp(1:NBocc(ik),ik))
                   esp_cb_min(ik)=minval(esp(NBocc(ik)+1:NB,ik))
                   esp_cb_max(ik)=maxval(esp(NBocc(ik)+1:NB,ik))
                end do
                if (minval(esp_cb_min(:))-maxval(esp_vb_max(:))<0.d0) then
                   call Occupation_Redistribution
                else
                   if (comm_is_root(nproc_id_global)) then
                      write(*,*) '======================================='
                      write(*,*) 'occupation redistribution is not needed'
                      write(*,*) '======================================='
                   end if
                end if
             end if
          end if
       else if( iter /= 1 )then ! sato
          call Fermi_Dirac_distribution
          if((comm_is_root(nproc_id_global)).and.(iter == Nscf))then
             open(126,file='occ.out')
             do ik=1,NK
                do ib=1,NB
                   write(126,'(2I7,e26.16E3)')ik,ib,occ(ib,ik)
                end do
             end do
             close(126)
          end if
       end if
       call Gram_Schmidt
       call diag_omp
       call Gram_Schmidt
       call CG_omp(Ncg)
       call Gram_Schmidt
       
       call psi_rho_GS
       call Density_Update(iter) 
       call Hartree

       functional_t = functional
       if(functional_t == 'TBmBJ' .and. iter < 20) functional = 'PZM'
       if(functional_t == 'BJ_PW' .and. iter < 20) functional = 'PZM'
       call Exc_Cor(calc_mode_gs,NBoccmax,zu_t)
       if(functional_t == 'TBmBJ' .and. iter < 20) functional = 'TBmBJ'
       if(functional_t == 'BJ_PW' .and. iter < 20) functional = 'BJ_PW'

       Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
       call Total_Energy_omp(Rion_update,calc_mode_gs)
       call Ion_Force_omp(Rion_update,calc_mode_gs)
       call sp_energy_omp
       call current_GS
       Eall_GS(iter)=Eall
       esp_var_ave(iter)=sum(esp_var(:,:))/(NK*Nelec/2)
       esp_var_max(iter)=maxval(esp_var(:,:))
       dns_diff(iter)=sqrt(sum((rho_out(:,iter)-rho_in(:,iter))**2))*Hxyz
       
       if (comm_is_root(nproc_id_global)) then
          write(*,*) 'Total Energy = ',Eall_GS(iter),Eall_GS(iter)-Eall_GS(iter-1)
          write(*,'(a28,3e15.6)') 'jav(1),jav(2),jav(3)= ',jav(1),jav(2),jav(3)
          write(*,'(4(i3,f12.6,2x))') (ib,esp(ib,1),ib=1,NB)
          do ia=1,NI
             write(*,'(1x,i7,3f15.6)') ia,force(1,ia),force(2,ia),force(3,ia)
          end do
          write(*,*) 'var_ave,var_max=',esp_var_ave(iter),esp_var_max(iter)
          write(*,*) 'dns. difference =',dns_diff(iter)
          if (iter/20*20 == iter) then
             write(*,*) '====='
             call timer_show_current_min('elapse time=',LOG_ALL)
          end if
          write(*,*) '-----------------------------------------------'
       end if
    end do
    call timer_end(LOG_GROUND_STATE)
    
    if(comm_is_root(nproc_id_global)) then
       call timer_show_hour('Ground State time  :', LOG_GROUND_STATE)
       call timer_show_min ('CG time            :', LOG_CG)
       call timer_show_min ('Gram Schmidt time  :', LOG_GRAM_SCHMIDT)
       call timer_show_min ('diag time          :', LOG_DIAG)
       call timer_show_min ('sp_energy time     :', LOG_SP_ENERGY)
       call timer_show_min ('hpsi time          :', LOG_HPSI)
       call timer_show_min (' - stencil time    :', LOG_HPSI_STENCIL)
       call timer_show_min (' - pseudo pt. time :', LOG_HPSI_PSEUDO)
       call timer_show_min ('psi_rho time       :', LOG_PSI_RHO)
       call timer_show_min ('Hartree time       :', LOG_HARTREE)
       call timer_show_min ('Exc_Cor time       :', LOG_EXC_COR)
       call timer_show_min ('current time       :', LOG_CURRENT)
       call timer_show_min ('Total_Energy time  :', LOG_TOTAL_ENERGY)
       call timer_show_min ('Ion_Force time     :', LOG_ION_FORCE)
    end if
    if(comm_is_root(nproc_id_global)) then
       write(*,*) 'This is the end of GS calculation'
       call timer_show_current_hour('elapse time=',LOG_ALL)
       write(*,*) '-----------------------------------------------------------'
    end if
    
    zu_GS0(:,:,:)=zu_GS(:,:,:)
    
    zu_t(:,:,:)=zu_GS(:,1:NBoccmax,:)
    Rion_eq=Rion
    dRion(:,:,-1)=0.d0; dRion(:,:,0)=0.d0
    
    call psi_rho_GS
    call Hartree
    call Exc_Cor(calc_mode_gs,NBoccmax,zu_t)

    Vloc(1:NL)=Vh(1:NL)+Vpsl(1:NL)+Vexc(1:NL)
    Vloc_GS(:)=Vloc(:)
    call Total_Energy_omp(Rion_update,calc_mode_gs)
    Eall0=Eall
    if(comm_is_root(nproc_id_global)) write(*,*) 'Eall =',Eall
    
    call timer_end(LOG_STATIC)
    if (comm_is_root(nproc_id_global)) then
       write(*,*) '-----------------------------------------------'
       call timer_show_min('static time=',LOG_STATIC)
       write(*,*) '-----------------------------------------------'
    end if
    
    if (comm_is_root(nproc_id_global)) then
       write(*,*) '-----------------------------------------------'
       write(*,*) '----some information for Band map--------------'
       do ik=1,NK 
          esp_vb_min(ik)=minval(esp(1:NBocc(ik),ik))
          esp_vb_max(ik)=maxval(esp(1:NBocc(ik),ik))
          esp_cb_min(ik)=minval(esp(NBocc(ik)+1:NB,ik))
          esp_cb_max(ik)=maxval(esp(NBocc(ik)+1:NB,ik))
       end do
       write(*,*) 'Bottom of VB',minval(esp_vb_min(:))
       write(*,*) 'Top of VB',maxval(esp_vb_max(:))
       write(*,*) 'Bottom of CB',minval(esp_cb_min(:))
       write(*,*) 'Top of CB',maxval(esp_cb_max(:))
       write(*,*) 'The Bandgap',minval(esp_cb_min(:))-maxval(esp_vb_max(:))
       write(*,*) 'BG between same k-point',minval(esp_cb_min(:)-esp_vb_max(:))
       write(*,*) 'Physicaly upper bound of CB for DOS',minval(esp_cb_max(:))
       write(*,*) 'Physicaly upper bound of CB for eps(omega)',minval(esp_cb_max(:)-esp_vb_min(:))
       write(*,*) '-----------------------------------------------'
       write(*,*) '-----------------------------------------------'
    end if
    
    call write_GS_data
    
    deallocate(rho_in,rho_out)
    deallocate(Eall_GS,esp_var_ave,esp_var_max,dns_diff)

  contains
    subroutine reset_gs_timer
      use timer
      implicit none
      integer :: i
      do i = LOG_CG,LOG_GRAM_SCHMIDT
         call timer_reset(i)
      end do
      call reset_rt_timer
    end subroutine reset_gs_timer

    subroutine reset_rt_timer
      implicit none
      integer :: i
      do i = LOG_DT_EVOLVE,LOG_ALLREDUCE
         call timer_reset(i)
      end do
    end subroutine reset_rt_timer
  end subroutine calc_ground_state


  
end module ground_state
