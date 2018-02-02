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
!This file is "k_shift_wf.f90"
!This file contain a subroutine.
!Subroutine k_shift_wf(iter,iter_GS_max)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
Subroutine write_projection_header
  use Global_Variables, only: NB
  use inputoutput, only: t_unit_time, t_unit_energy
  implicit none

 !open(404, file=file_ovlp,position = position_option) 
  write(404, '("#",1X,A)') "Projection"
  write(404, '("#",1X,A,":",1X,A)') "ik", "k-point index"
  write(404, '("#",1X,A,":",1X,A)') "ovlp_occup", "Occupation"
  write(404, '("#",1X,A,":",1X,A)') "NB", "Number of bands"
  write(404, '("#",99(1X,I0,":",A,"[",A,"]"))') &
  &           1, "ik", "none", &
  &           2, "ovlp_occup(NB)", "none"
        
 !open(408, file=file_nex, position = position_option) 
   write(408, '("#",1X,A)') "Excitation"
   write(408, '("#",1X,A,":",1X,A)') "nelec", "Number of excited electrons"
   write(408, '("#",1X,A,":",1X,A)') "nhole", "Number of excited holes"
   write(408, '("#",99(1X,I0,":",A,"[",A,"]"))')  &
   &           1, "time", trim(t_unit_time%name), &
   &           2, "nelec", "none", &
   &           3, "nhole", "none"
      
  !open(409, file=file_last_band_map,position = position_option) 
   write(409, '("#",1X,A)') "Last bandmap"
   write(409, '("#",1X,A,":",1X,A)') "ik", "k-point index"
   write(409, '("#",1X,A,":",1X,A)') "energy", "Electron energy"
   write(409, '("#",1X,A,":",1X,A)') "ovlp_occup", "Occupation"
   write(409, '("#",1X,A,":",1X,A)') "NB", "Number of bands"
   write(409, '("#",99(1X,I0,":",A,"[",A,"]"))') &
   &           1, "ik", "none", &
   &           2, "energy(NB)", trim(t_unit_energy%name), &
   &           2 + NB, "ovlp_occup(NB)", "none"      

  return
End Subroutine write_projection_header

Subroutine analysis_RT_using_GS(Rion_xyz_update,iter_GS_max,zu,it,action)
!(this subroutine was named "k_shift_wf" in past)
  use Global_Variables
  use salmon_parallel, only: nproc_group_tdks, nproc_id_global
  use salmon_communication, only: comm_summation, comm_is_root
  use inputoutput, only: t_unit_time, t_unit_energy
  implicit none
  integer,intent(in) :: iter_GS_max
  logical,intent(in) :: Rion_xyz_update
  complex(8),intent(in) :: zu(NL,NBoccmax,NK_s:NK_e)
  real(8) :: esp_all(NB,NK), Eall_prev,dEall, threshold_dEall
  real(8),allocatable :: rho_backup(:)
  integer :: iter_GS,ik,ib1,ib2,ib,ia,it
  character(10) :: action

  threshold_dEall = 1.0d-6  !now fixed, not controlled by input

  if(comm_is_root(nproc_id_global)) &
     write(*,*)'  Analysis option: ',action

  !(this condition may be changed for md option in futre)
  if(projection_option=='gs' .and.  use_ehrenfest_md=='n')then
    Vloc_t(:)=Vloc(:)
    Vloc(:)=Vloc_GS(:)
  end if
  zu_GS(:,:,:)=zu_GS0(:,:,:)

  !(GS iteration)
  Eall_prev=Eall_GS0

  call Total_Energy_omp(Rion_xyz_update,calc_mode_gs)
  dEall = Eall - Eall_prev
  if(abs(dEall).le.1d-12) then
     if(comm_is_root(nproc_id_global)) then
          write(*,*) "  No GS iteration: already converged"
          write(*,*) "  dEall=",real(dEall)
     endif
     goto 1
  endif

  if(comm_is_root(nproc_id_global)) &
     write(*,*)'  iter_GS    Eall-Eall0           dEall'

  do iter_GS=1,iter_GS_max
   !(this if-condition may be changed later)
   !if(use_ehrenfest_md=='y')then !for better conv
    if(threshold_dEall.le.1d-5)then !for better conv
       call diag_omp
       call Gram_Schmidt
    else
       if(iter_GS==1) call diag_omp
    endif
    call CG_omp(Ncg)
    call Gram_Schmidt
    call Total_Energy_omp(Rion_xyz_update,calc_mode_gs)
    call Ion_Force_omp(Rion_xyz_update,calc_mode_gs)
    dEall = Eall - Eall_prev
    if(comm_is_root(nproc_id_global)) &
       write(*,'(i6,2e20.10)') iter_GS,Eall-Eall0,dEall
    if(iter_GS.ge.3 .and. abs(dEall).le.threshold_dEall) exit
    Eall_prev = Eall
  end do
1 continue

  if(action=="proj_last ") then
     esp=0d0
     call diag_omp
     call comm_summation(esp,esp_all,NK*NB,nproc_group_tdks)
  endif

  !(analysis of projection_option)
  if(action(1:4)=="proj") then

     ovlp_occ_l=0.d0
!$omp parallel do private(ik,ib1,ib2)
     do ik=Nk_s,NK_e
     do ib1=1,NB
     do ib2=1,NBoccmax
        ovlp_occ_l(ib1,ik)=ovlp_occ_l(ib1,ik)+occ(ib2,ik)*abs(sum(conjg(zu_GS(:,ib1,ik))*zu(:,ib2,ik))*Hxyz)**2
     enddo
     enddo
     enddo
     call comm_summation(ovlp_occ_l,ovlp_occ,NK*NB,nproc_group_tdks)

     if(comm_is_root(nproc_id_global)) then

        if(action=="proj_last ") then
           do ik=1,NK
              write(409,'(i6,1000e26.16E3)')ik,(esp_all(ib,ik)*t_unit_energy%conv,ovlp_occ(ib,ik)*NKxyz,ib=1,NB)
           end do
        else &
        if(action=="projection") then
           write(404,'(i11)') it
           do ik=1,NK
              write(404,'(i6,500e16.6)')ik,(ovlp_occ(ib,ik)*NKxyz,ib=1,NB)
           enddo
           write(408,'(1x,3e16.6E3)')it*dt**t_unit_time%conv,sum(ovlp_occ(NBoccmax+1:NB,:)),sum(occ)-sum(ovlp_occ(1:NBoccmax,:))
        end if

        write(*,*) 'forces on atoms:'
        do ia=1,NI
           write(*,'(i8,3f15.6)') ia,force(1,ia),force(2,ia),force(3,ia)
        end do
        write(*,*) 'number of excited electron',sum(ovlp_occ(NBoccmax+1:NB,:)),sum(occ)-sum(ovlp_occ(1:NBoccmax,:))
        write(*,*) 'var_tot,var_max=',sum(esp_var(:,:))/(NK*Nelec/2),maxval(esp_var(:,:)) 
     end if

  endif  !action=="projection" or "proj_last "


  !(get density of ground state)
  if(action=="get_dns_gs") then
     if(.not.allocated(rho_gs_t)) allocate(rho_gs_t(NL))
     allocate(rho_backup(NL))
     rho_backup(:) = rho(:)
     call psi_rho_GS
     rho_gs_t(:) = rho(:)
     rho(:) = rho_backup(:)
     deallocate(rho_backup)
  endif

  !(this condition may be changed for md option in futre)
  if(projection_option=='gs' .and.  use_ehrenfest_md=='n')then
     Vloc(:)=Vloc_t(:)
  end if

  zu_GS0(:,:,:)=zu_GS(:,:,:) !update for the next initial guess
  Eall_GS0 = Eall

  return
End Subroutine analysis_RT_using_GS
