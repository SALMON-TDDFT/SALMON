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
!=======================================================================
!=======================================================================

SUBROUTINE time_evolution_step(shtpsi)
use salmon_parallel, only: nproc_id_global, nproc_group_global, nproc_group_grid, nproc_group_h, nproc_group_korbital
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use misc_routines, only: get_wtime
use inputoutput
use scf_data
use new_world_sub
use allocate_mat_sub
use read_pslfile_sub

implicit none
integer :: ix,iy,iz,i1,mm,jj
integer :: ii,iob,iatom,iik
real(8) :: rbox1,rbox1q,rbox1q12,rbox1q23,rbox1q31,rbox1e
complex(8),allocatable :: cmatbox1(:,:,:),cmatbox2(:,:,:)
real(8) :: absr2

integer :: idensity, idiffDensity, ielf
real(8) :: rNe

complex(8),parameter :: zi=(0.d0,1.d0)

complex(8) :: shtpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd,   &
                     1:iobnum,k_sta:k_end)

complex(8) :: cbox1,cbox2,cbox3

elp3(511)=get_wtime()

idensity=0
idiffDensity=1
ielf=2 

if(iperiodic==3) call init_k_rd(k_rd,ksquare,1)

select case(ikind_eext)
  case(0,3,9:12)
    ihpsieff=0
  case(1,2,4,6:8,15)
    ihpsieff=1
end select

if(ikind_eext==1) call calcVbox

elp3(512)=get_wtime()
elp3(532)=elp3(532)+elp3(512)-elp3(511)

if(iobnum.ge.1)then
  if(mod(itt,2)==1)then
    call taylor(zpsi_in,zpsi_out,shtpsi)
  else
    call taylor(zpsi_out,zpsi_in,shtpsi)
  end if
end if

if(iperiodic==0)then
  if(ikind_eext==0.and.itt>=2)then
    if(mod(itt,2)==1)then
      call Total_energy_groupob(zpsi_out,shtpsi,2)
    else
      call Total_energy_groupob(zpsi_in,shtpsi,2)
    end if
    call subdip(rNe,2)
  end if
end if

elp3(513)=get_wtime()
elp3(533)=elp3(533)+elp3(513)-elp3(512)

if(ilsda == 0) then

  elp3(761)=get_wtime()
  call comm_summation(rhobox,rho,mg_num(1)*mg_num(2)*mg_num(3),nproc_group_grid)
  elp3(762)=get_wtime()
  elp3(760)=elp3(760)+elp3(762)-elp3(761)

else if(ilsda==1)then

  elp3(761)=get_wtime()
  call comm_summation(rhobox_s,rho_s,mg_num(1)*mg_num(2)*mg_num(3)*2,nproc_group_grid)
!$OMP parallel do private(iz,iy,ix)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rho(ix,iy,iz)=rho_s(ix,iy,iz,1)+rho_s(ix,iy,iz,2)
  end do
  end do
  end do
  elp3(762)=get_wtime()
  elp3(760)=elp3(760)+elp3(762)-elp3(761)
end if

  elp3(515)=get_wtime()
   if(itt/=1)then
     if(mod(itt,2)==1)then
!$OMP parallel do private(iz,iy,ix)
       do iz=ng_sta(3),ng_end(3)
       do iy=ng_sta(2),ng_end(2)
       do ix=ng_sta(1),ng_end(1)
         Vh_stock2(ix,iy,iz)=2.d0*Vh_stock1(ix,iy,iz)-Vh_stock2(ix,iy,iz)
       end do
       end do
       end do
     else
!$OMP parallel do private(iz,iy,ix)
       do iz=ng_sta(3),ng_end(3)
       do iy=ng_sta(2),ng_end(2)
       do ix=ng_sta(1),ng_end(1)
         Vh_stock1(ix,iy,iz)=2.d0*Vh_stock2(ix,iy,iz)-Vh_stock1(ix,iy,iz)
       end do
       end do
       end do
     end if
   end if

  
  call Hartree_ns

  elp3(516)=get_wtime()
  elp3(536)=elp3(536)+elp3(516)-elp3(515)

  if(imesh_s_all==1.or.(imesh_s_all==0.and.nproc_id_global<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
    call exc_cor_ns
  end if

  elp3(517)=get_wtime()
  elp3(537)=elp3(537)+elp3(517)-elp3(516)

  call allgatherv_vlocal

  elp3(518)=get_wtime()
  elp3(538)=elp3(538)+elp3(518)-elp3(517)

! result

  if(iperiodic==0)then
    if(ikind_eext/=0.or.(ikind_eext==0.and.itt==itotNtime))then
      elp3(526)=get_wtime()
  
      ihpsieff=0
      if(mod(itt,2)==1)then
        call Total_Energy_groupob(zpsi_out,shtpsi,1)              ! Total energy
      else
        call Total_Energy_groupob(zpsi_in,shtpsi,1)              ! Total energy
      end if
      elp3(527)=get_wtime()
      elp3(540)=elp3(540)+elp3(527)-elp3(526)
      
      call subdip(rNe,1)
      elp3(528)=get_wtime()
      elp3(539)=elp3(539)+elp3(528)-elp3(527)
  
    end if
  end if

  if(iwrite_projection==1.and.mod(itt,itwproj)==0)then
    if(mod(itt,2)==1)then
      call projection(zpsi_out)
    else
      call projection(zpsi_in)
    end if
  end if

  if(iflag_dip2==1)then
    do jj=1,num_dip2
      do i1=1,3
        rbox1=0.d0
!$OMP parallel do reduction( + : rbox1 ) private(iz,iy,ix)
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rbox1=rbox1+vecR(i1,ix,iy,iz)*rho(ix,iy,iz)*rto_ix(ix,jj)
        end do
        end do
        end do
        rbox_array_dip2(i1,jj)=rbox1
      end do
    end do
  
    call comm_summation(rbox_array_dip2,rbox_array2_dip2,4*num_dip2,nproc_group_h)

    do jj=1,num_dip2
      Dp2(1:3,itt,jj)=rbox_array2_dip2(1:3,jj)*Hgs(1:3)*Hvol-vecDs2(1:3,jj)
    end do

!------------QUADRUPOLE-start------------

    if(quadrupole=='y')then
      rho_diff(:,:,:) = rho(:,:,:)-rho0(:,:,:)
      do jj=1,num_dip2
        vecR_tmp(:,:,:,:)=vecR(:,:,:,:)
        vecR_tmp(1,:,:,:)=vecR_tmp(1,:,:,:)-dip2center(jj)/Hgs(1)
        do i1=1,3
          rbox1q=0.d0
 !$OMP parallel do reduction( + : rbox1q ) private(absr2,iz,iy,ix)
          do iz=ng_sta(3),ng_end(3)
          do iy=ng_sta(2),ng_end(2)
          do ix=ng_sta(1),ng_end(1)
            absr2=vecR_tmp(1,ix,iy,iz)**2+vecR_tmp(2,ix,iy,iz)**2+vecR_tmp(3,ix,iy,iz)**2
            rbox1q=rbox1q+(3.d0*vecR_tmp(i1,ix,iy,iz)*vecR_tmp(i1,ix,iy,iz)-absr2)*rho_diff(ix,iy,iz)*rto_ix(ix,jj)
          end do
          end do
          end do
          rbox_array_dip2q(i1,i1,jj)=rbox1q
        end do
      end do
    
      do jj=1,num_dip2
        rbox1q12=0.d0
        rbox1q23=0.d0
        rbox1q31=0.d0
 !$OMP parallel do reduction( + : rbox1q12,rbox1q23,rbox1q31 ) private(iz,iy,ix)
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rbox1q12=rbox1q12+3.d0*vecR_tmp(1,ix,iy,iz)*vecR_tmp(2,ix,iy,iz)*rho_diff(ix,iy,iz)*rto_ix(ix,jj)
          rbox1q23=rbox1q23+3.d0*vecR_tmp(2,ix,iy,iz)*vecR_tmp(3,ix,iy,iz)*rho_diff(ix,iy,iz)*rto_ix(ix,jj)
          rbox1q31=rbox1q31+3.d0*vecR_tmp(3,ix,iy,iz)*vecR_tmp(1,ix,iy,iz)*rho_diff(ix,iy,iz)*rto_ix(ix,jj)
        end do
        end do
        end do
        rbox_array_dip2q(1,2,jj)=rbox1q12 ; rbox_array_dip2q(2,1,jj)=rbox1q12
        rbox_array_dip2q(2,3,jj)=rbox1q23 ; rbox_array_dip2q(3,2,jj)=rbox1q23
        rbox_array_dip2q(3,1,jj)=rbox1q31 ; rbox_array_dip2q(1,3,jj)=rbox1q31
      end do

      call comm_summation(rbox_array_dip2q,rbox_array2_dip2q,9*num_dip2,nproc_group_h)
 
      do jj=1,num_dip2
        do i1=1,3
          Qp2(1:3,i1,itt,jj)=rbox_array2_dip2q(1:3,i1,jj)*Hgs(1:3)**2*Hvol
        end do
      end do

    end if

!------------QUADRUPOLE-end--------------
    if(iflag_intelectron==1)then
    do jj=1,num_dip2
        rbox1e=0.d0
!$OMP parallel do reduction( + : rbox1e ) private(iz,iy,ix)
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rbox1e=rbox1e+rho(ix,iy,iz)*rto_ix(ix,jj)
        end do
        end do
        end do
        rbox_array_dip2e(jj)=rbox1e
      end do
  
      call comm_summation(rbox_array_dip2e,rbox_array2_dip2e,num_dip2,nproc_group_h)

      do jj=1,num_dip2
        rIe2(itt,jj)=rbox_array2_dip2e(jj)*Hvol
      end do
    end if
  end if


  if(iperiodic==3)then
    call subdip(rNe,1)
    if(mod(itt,2)==1)then
      elp3(526)=get_wtime()
      call calc_current(zpsi_out)
      elp3(527)=get_wtime()
      if(itt==1.or.itt==itotNtime.or.mod(itt,itcalc_ene)==0)then
        call Total_Energy_periodic(zpsi_out,shtpsi)              ! Total energy
      end if
      elp3(528)=get_wtime()
      elp3(539)=elp3(539)+elp3(527)-elp3(526)
      elp3(540)=elp3(540)+elp3(528)-elp3(527)
    else
      elp3(526)=get_wtime()
      call calc_current(zpsi_in)
      elp3(527)=get_wtime()
      if(itt==itotNtime.or.mod(itt,itcalc_ene)==0)then
        call Total_Energy_periodic(zpsi_in,shtpsi)              ! Total energy
      end if
      elp3(528)=get_wtime()
      elp3(539)=elp3(539)+elp3(527)-elp3(526)
      elp3(540)=elp3(540)+elp3(528)-elp3(527)
    end if
    rbox1=0.d0
  !$OMP parallel do private(iz,iy,ix) reduction( + : rbox1 )
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      rbox1=rbox1+rho(ix,iy,iz)
    end do
    end do
    end do
    rbox1=rbox1*Hvol
    call comm_summation(rbox1,rNe,nproc_group_h)
  !    write(*,'(1x,i7, 3e16.8, f15.8,f18.8,i5,f16.8)')       &
  !      itt, (curr(i1,itt),i1=1,3), Ne, Etot*2d0*Ry,iterVh,dble(cumnum)
    if(comm_is_root(nproc_id_global))then
      write(*,'(i8,f14.8, 3e16.8, f15.8,f18.8)')       &
        itt,dble(itt)*dt*2.41888d-2, (curr(i1,itt),i1=1,3), rNe, Etot*2d0*Ry
      write(16,'(f14.8, 3e16.8, f15.8,f18.8)')       &
        dble(itt)*dt*2.41888d-2, (curr(i1,itt),i1=1,3), rNe, Etot*2d0*Ry
      write(17,'(f14.8, 3e16.8)')       &
        dble(itt)*dt*2.41888d-2, (E_tot(i1,itt),i1=1,3)
      write(18,'(f14.8, 3e16.8)')       &
        dble(itt)*dt*2.41888d-2, (E_ext(i1,itt),i1=1,3)
      write(19,'(f14.8, 3e16.8)')       &
        dble(itt)*dt*2.41888d-2, (E_ind(i1,itt),i1=1,3)
    end if
  end if

  elp3(520)=get_wtime()

  if(icalcforce==1)then
    if(mod(itt,2)==1)then
      call calc_force_c(zpsi_out)
    else
      call calc_force_c(zpsi_in)
    end if
    if(comm_is_root(nproc_id_global))then
      do iatom=1,MI
        dRion(:,iatom,1)=2*dRion(:,iatom,0)-dRion(:,iatom,-1)+    &
                             rforce(:,iatom)*dt**2/(umass*Mass(Kion(iatom)))
        Rion(:,iatom)=Rion_eq(:,iatom)+dRion(:,iatom,1)
      enddo
    end if
    call comm_bcast(Rion,nproc_group_global)
    call init_ps
    if(comm_is_root(nproc_id_global))then
      dRion(:,:,-1)=dRion(:,:,0)
      dRion(:,:,0)=dRion(:,:,1)
    end if
    call comm_bcast(dRion,nproc_group_global)
  end if


  if(circular=='y')then
    allocate(cmatbox1(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
    allocate(cmatbox2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
    
!$OMP parallel do private(iz,iy,ix)
    do iz=lg_sta(3),lg_end(3)
    do iy=lg_sta(2),lg_end(2)
    do ix=lg_sta(1),lg_end(1)
      cmatbox1(ix,iy,iz)=0.d0
    end do
    end do
    end do
    cbox1=0.d0

    do iik=k_sta,k_end
    do iob=1,iobnum
      if(mod(itt,2)==1)then
!$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cmatbox1(ix,iy,iz)=zpsi_out(ix,iy,iz,iob,iik)
        end do
        end do
        end do
      else
!$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cmatbox1(ix,iy,iz)=zpsi_in(ix,iy,iz,iob,iik)
        end do
        end do
        end do
      end if

      call comm_summation(cmatbox1,cmatbox2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_korbital)
      cbox3=0.d0
      do iz=lg_sta(3),lg_end(3)
      do ix=1,lg_end(1)
        cbox3=cbox3+(conjg(cmatbox2(ix,0,iz))*(cmatbox2(ix,1,iz)-cmatbox2(ix,-1,iz))/(2.d0*Hgs(2)) &
                   -(conjg(cmatbox2(ix,1,iz))-conjg(cmatbox2(ix,-1,iz)))/(2.d0*Hgs(2))*cmatbox2(ix,0,iz))/Hgs(1)/Hgs(2)
      end do
      end do
      cbox1=cbox1+cbox3

      cbox3=0.d0
      do iz=lg_sta(3),lg_end(3)
      do iy=1,lg_end(2)
        cbox3=cbox3-(conjg(cmatbox2(0,iy,iz))*(cmatbox2(1,iy,iz)-cmatbox2(-1,iy,iz))/(2.d0*Hgs(1)) &
                   +(conjg(cmatbox2(1,iy,iz))-conjg(cmatbox2(-1,iy,iz)))/(2.d0*Hgs(1))*cmatbox2(0,iy,iz))/Hgs(1)/Hgs(2)
      end do
      end do
      cbox1=cbox1+cbox3

    end do
    end do

    call comm_summation(cbox1,cbox2,nproc_group_global)

    cumnum=cumnum+cbox2/zi*dt

    deallocate(cmatbox1,cmatbox2) 
  end if

  if(comm_is_root(nproc_id_global))then
    if(iflag_md==1)then
      write(15,'(3f16.8)') dble(itt)*dt*0.0241889d0,  &
                  sqrt((Rion(1,idisnum(1))-Rion(1,idisnum(2)))**2   &
                      +(Rion(2,idisnum(1))-Rion(2,idisnum(2)))**2   &
                      +(Rion(3,idisnum(1))-Rion(3,idisnum(2)))**2)*a_B, Etot*2.d0*Ry
      do ii=1,wmaxMI
        write(20+ii,'(4f16.8)') dble(itt)*dt*0.0241889d0, (Rion(jj,ii)*a_B,jj=1,3)
        write(30+ii,'(4f16.8)') dble(itt)*dt*0.0241889d0, (rforce(jj,ii)*2.d0*Ry/a_B,jj=1,3)
      end do
    end if
  end if

  if(iflag_fourier_omega==1)then
    do mm=1,num_fourier_omega
!$OMP parallel do  private(iz,iy,ix)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        zalpha2(ix,iy,iz,mm)=zalpha2(ix,iy,iz,mm)   &
                             +exp(zi*fourier_omega(mm)*(itt*dt))*(rho(ix,iy,iz)-rho0(ix,iy,iz)) & 
                             *(1-3*(itt/itotNtime2)**2+2*(itt/itotNtime2)**3)
      end do
      end do
      end do
    end do
  end if

  if(out_dns_rt=='y')then
    if(mod(itt,out_dns_rt_step)==0)then
      call writedns
    end if
  end if
  if(out_elf_rt=='y')then
    if(mod(itt,out_elf_rt_step)==0)then
      call calcELF
      call writeelf
    end if
  end if
  if(out_estatic_rt=='y')then
    if(mod(itt,out_estatic_rt_step)==0)then
      call calcEstatic
      call writeestatic
    end if
  end if

  elp3(521)=get_wtime()
  elp3(541)=elp3(541)+elp3(521)-elp3(520)
  elp3(542)=elp3(542)+elp3(521)-elp3(511)



return

END SUBROUTINE time_evolution_step
