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
  use inputoutput, only: t_unit_time, t_unit_energy, &
  &                      projection_option,projection_decomp
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

 !open(420, file=file_nex_atom,position = position_option) 
 if(projection_option=='gs' .and. projection_decomp=='atom') then
  write(420,'("#",1X,A)') "Excitation"
  write(420,'("#",1X,A,":",1X,A)') "nelec-ia","Number of excited electrons on ia-th atom"
  write(420,'("#",99(1X,I0,":",A,"[",A,"]"))')  &
  &          1, "time", trim(t_unit_time%name), &
  &          2, "sum(nelec-i)", "none", &
  &          3, "nelec-1", "none", &
  &          4, "nelec-2", "none", &
  &          5, ".......", ""
  endif

  return
End Subroutine write_projection_header

Subroutine analysis_RT_using_GS(Rion_xyz_update,iter_GS_max,zu,it,action)
!(this subroutine was named "k_shift_wf" in past)
  use Global_Variables
  use salmon_parallel, only: nproc_group_tdks, nproc_id_global
  use salmon_communication, only: comm_summation, comm_is_root
  use inputoutput, only: t_unit_time, t_unit_energy
  use ground_state
  implicit none
  integer,intent(in) :: iter_GS_max
  logical,intent(in) :: Rion_xyz_update
  complex(8),intent(in) :: zu(NL,NBoccmax,NK_s:NK_e)
  real(8) :: esp_all(NB,NK), Eall_prev,dEall, threshold_dEall, Vloc_save(NL),rtmp
  real(8) :: time, nee, neh, var_tot, var_max, tconv
  real(8),allocatable :: rho_backup(:), nee_ia(:)
  complex(8),allocatable :: ovlp_occ_ia_l(:,:,:,:),ovlp_occ_ia(:,:,:,:)
  integer :: i, iter_GS,ik,ib1,ib2,ib,ia,it, ia1,ia2
  character(10) :: action

  threshold_dEall = 1.0d-6  !now fixed, not controlled by input

  if(comm_is_root(nproc_id_global)) &
     write(*,*)'  Analysis option: ',action

  !(this condition may be changed for md option in futre)
  Vloc_save(:) = Vloc(:)
  if(use_ehrenfest_md=='n' .and. projection_option=='gs') Vloc(:)=Vloc_GS(:)
  zu_GS(:,:,:)=zu_GS0(:,:,:)

  if(it .ne. it_last_update_zu_GS_proj) then
  if(use_ehrenfest_md=='y'.and.(projection_option=='gs'.or.action=="get_dns_gs"))then
     !(some variables are changed assuming this subroutine is called only in RT)
     iflag_gs_init_wf = 1
     PrLv_scf = 0
     flag_scf_conv_ene_force=.true.
     flag_update_only_zu_GS =.true.
     rtmp           = convrg_scf_ene
     convrg_scf_ene = threshold_dEall
     call calc_ground_state
     convrg_scf_ene = rtmp
     flag_update_only_zu_GS =.false.
  else
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
1   continue

    if(action=="proj_last ") then
       esp=0d0
       call diag_omp
       call comm_summation(esp,esp_all,NK*NB,nproc_group_tdks)
    endif

  endif
  endif
  it_last_update_zu_GS_proj=it

  !(analysis of projection_option)
  if(action(1:4)=="proj") then

     if(projection_option=='gs' .and. projection_decomp=='atom') then
        !initialize
        if(.not.allocated(assign_grid_atom)) call init_proj_decomp_atom
     endif

     !calc overlap(projection)
     ovlp_occ_l=0.d0
!$omp parallel do private(ik,ib1,ib2)
     do ik=Nk_s,NK_e
     do ib1=1,NB
     do ib2=1,NBoccmax
        ovlp_occ_l(ib1,ik) = ovlp_occ_l(ib1,ik) &
        &  + occ(ib2,ik) * abs(sum(conjg(zu_GS(:,ib1,ik))*zu(:,ib2,ik))*Hxyz)**2
     enddo
     enddo
     enddo
     call comm_summation(ovlp_occ_l,ovlp_occ,NK*NB,nproc_group_tdks)

     !print
     if(comm_is_root(nproc_id_global)) then

        if(action=="proj_last ") then
           tconv = t_unit_energy%conv
           do ik=1,NK
              write(409,10)ik,(esp_all(ib,ik)*tconv,ovlp_occ(ib,ik)*NKxyz,ib=1,NB)
           end do
10         format(i6,1000e26.16E3)

        else &
        if(action=="projection") then
           write(404,'(i11)') it
           do ik=1,NK
              write(404,20) ik,(ovlp_occ(ib,ik)*NKxyz,ib=1,NB)
           enddo
           time = it*dt**t_unit_time%conv
           nee  = sum(ovlp_occ(NBoccmax+1:NB,:))
           neh  = sum(occ)-sum(ovlp_occ(1:NBoccmax,:))
           write(408,30) time, nee, neh
20         format(i6,500e16.6)
30         format(1x,3e16.6E3)
        end if

        write(*,*) '  forces on atoms (GS):'
        do ia=1,NI
           write(*,'(i8,3f15.6)') ia,force(1,ia),force(2,ia),force(3,ia)
        end do
        nee  = sum(ovlp_occ(NBoccmax+1:NB,:))
        neh  = sum(occ)-sum(ovlp_occ(1:NBoccmax,:))
        write(*,*) 'number of excited electron', nee, neh
        var_tot = sum(esp_var(:,:))/(NK*Nelec/2)
        var_max = maxval(esp_var(:,:)) 
        write(*,*) 'var_tot,var_max=', var_tot, var_max
     end if

     if(action=="projection" .and. projection_decomp=='atom') then
        allocate(ovlp_occ_ia_l(NB,NB,NK,NI), ovlp_occ_ia(NB,NB,NK,NI),nee_ia(NI))

        ovlp_occ_ia_l=0.d0
!$omp parallel do private(ik,ib1,ib2,i,ia)
        do ik=Nk_s,NK_e
        do ib1=1,NB
        do ib2=1,NBoccmax
           do i=1,NL
              ia = assign_grid_atom(i)
              ovlp_occ_ia_l(ib1,ib2,ik,ia) = ovlp_occ_ia_l(ib1,ib2,ik,ia) &
              &  + conjg(zu_GS(i,ib1,ik)) * zu(i,ib2,ik)
           enddo
        enddo
        enddo
        enddo
        ovlp_occ_ia_l = ovlp_occ_ia_l * Hxyz
        call comm_summation(ovlp_occ_ia_l,ovlp_occ_ia,NB*NB*NK*NI,nproc_group_tdks)

        nee_ia(:) = 0d0
        do ia1=1,NI
           do ia2=1,NI
           do ib1 = NBoccmax+1,NB
           do ib2 = 1,NBoccmax
              nee_ia(ia1) = nee_ia(ia1) &
              &  + sum(occ(ib2,:)*real(conjg(ovlp_occ_ia(ib1,ib2,:,ia1))*ovlp_occ_ia(ib1,ib2,:,ia2)))
           enddo
           enddo
           enddo
        enddo

        !print
        if(comm_is_root(nproc_id_global)) then
           time = it*dt**t_unit_time%conv
           write(420,40) time, real(sum(nee_ia(:))), (nee_ia(ia),ia=1,NI)
40         format(1x,1000e16.6E3)
        endif

        deallocate(ovlp_occ_ia_l, ovlp_occ_ia, nee_ia)

     endif


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

  Vloc(:) = Vloc_save(:)
  zu_GS0(:,:,:) = zu_GS(:,:,:) !update for the next initial guess
  Eall_GS0 = Eall

  contains
    subroutine init_proj_decomp_atom
      implicit none
      integer :: i,j,ia,ia_close,ix,iy,iz, fh_ass
      real(8) :: x,y,z, r,r_close


      allocate(assign_grid_atom(NL))

      do i=1,NL

         r_close  = 1d99
         ia_close = -1
         do ix=-2,2
         do iy=-2,2
         do iz=-2,2

            do ia=1,NI
               x = Rion(1,ia)+ix*aLx - Lx(i)*Hx
               y = Rion(2,ia)+iy*aLy - Ly(i)*Hy
               z = Rion(3,ia)+iz*aLz - Lz(i)*Hz
               r = sqrt(x*x+y*y+z*z)
               if(r .le. r_close) then
                  r_close  = r
                  ia_close = ia
               endif
            enddo

         enddo
         enddo
         enddo
         assign_grid_atom(i) = ia_close

      enddo

      !for check (print assignment to cube file)
      if(comm_is_root(nproc_id_global)) then
      fh_ass = 123
      open(fh_ass,file="assign_grid_atom.cube",status="unknown")
      write(fh_ass,10)
      write(fh_ass,20) NI,  0d0, 0d0, 0d0
      write(fh_ass,20) NLx, Hx,  0d0, 0d0
      write(fh_ass,20) NLy, 0d0, Hy,  0d0
      write(fh_ass,20) NLz, 0d0, 0d0, Hz

      do i=1, NI
         write(fh_ass,30) Zatom(Kion(i)), 0d0, (Rion(j,i), j=1,3)
      end do

      j=1
      do ix=0, NLx-1
      do iy=0, NLy-1
      do iz=0, NLz-1
         i=Lxyz(ix,iy,iz)
         if(mod(j,6)==0) then
            write(fh_ass,40) dble(assign_grid_atom(i))
         else
            write(fh_ass,40,advance='no') dble(assign_grid_atom(i))
         endif
         j=j+1
      end do
      end do
      end do

10    format("# SALMON",/, &
      &      "# COMMENT" )
20    format(I5,3(F12.6))
30    format(I5,4(F12.6))
40    format(ES12.4)

      close(fh_ass)
      endif

    end subroutine init_proj_decomp_atom
!  return
End Subroutine analysis_RT_using_GS
