! Copyright 2017 Katsuyuki Nobusada, Masashi Noda, Kazuya Ishimura, Kenji Iida, Maiku Yamaguchi
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!=======================================================================
!=======================================================================

SUBROUTINE time_evolution_step(shtpsi)
!$ use omp_lib
use scf_data
use new_world_sub
use allocate_mat_sub
use read_pslfile_sub

implicit none
integer :: ix,iy,iz,i1,i2,mm,jj
integer :: ii,iob,is,iatom
integer :: ista_box(3),iend_box(3)
real(8) :: box,rbox1,rbox2,rbox3,rbox1q,rbox1q12,rbox1q23,rbox1q31,rbox1e
real(8) :: rbox_array(4),rbox_array2(4)
real(8),allocatable :: rbox_array3(:,:),rbox_array4(:,:)
complex(8),allocatable :: cmatbox1(:,:,:),cmatbox2(:,:,:)
real(8) :: rho_region_center
real(8) :: absr2

integer :: idensity, idiffDensity, ielf
real(8) :: rNe

complex(8),parameter :: zi=(0.d0,1.d0)

complex(8) :: shtpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd,   &
                     1:iobnum,1)

complex(8) :: cbox1,cbox2,cbox3

integer :: ix_sta_Vbox(3),ix_end_Vbox(3)
integer :: icount
integer :: jsta,jend

!$ call omp_set_num_threads(inumthreads)

elp3(511)=MPI_Wtime()

idensity=0
idiffDensity=1
ielf=2 

select case(ikind_eext)
  case(0,3,9:12)
    ihpsieff=0
  case(1,2,4,6:8,15)
    ihpsieff=1
end select

call calcVbox

elp3(512)=MPI_Wtime()
elp3(532)=elp3(532)+elp3(512)-elp3(511)

if(iobnum.ge.1)then
  if(mod(itt,2)==1)then
    call taylor(zpsi_in,zpsi_out,shtpsi)
  else
    call taylor(zpsi_out,zpsi_in,shtpsi)
  end if
end if

if(ikind_eext==0.and.itt>=2)then
  if(mod(itt,2)==1)then
    call Total_energy_groupob(zpsi_out,shtpsi,2)
  else
    call Total_energy_groupob(zpsi_in,shtpsi,2)
  end if
  call subdip(rNe,2)
end if

elp3(513)=MPI_Wtime()
elp3(533)=elp3(533)+elp3(513)-elp3(512)

if(ilsda == 0) then

  elp3(761)=MPI_Wtime()
  call MPI_allreduce(rhobox,rho,      &
           mg_num(1)*mg_num(2)*mg_num(3),      &
           MPI_DOUBLE_PRECISION,MPI_SUM,      &
           newworld_comm_grid,ierr)
  elp3(762)=MPI_Wtime()
  elp3(760)=elp3(760)+elp3(762)-elp3(761)

else if(ilsda==1)then

  elp3(761)=MPI_Wtime()
  call MPI_allreduce(rhobox_s,rho_s,      &
           mg_num(1)*mg_num(2)*mg_num(3)*2,      &
           MPI_DOUBLE_PRECISION,MPI_SUM,      &
           newworld_comm_grid,ierr)
!$OMP parallel do
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rho(ix,iy,iz)=rho_s(ix,iy,iz,1)+rho_s(ix,iy,iz,2)
  end do
  end do
  end do
  elp3(762)=MPI_Wtime()
  elp3(760)=elp3(760)+elp3(762)-elp3(761)
end if

  elp3(515)=MPI_Wtime()
   if(itt/=1)then
     if(mod(itt,2)==1)then
!$OMP parallel do
       do iz=ng_sta(3),ng_end(3)
       do iy=ng_sta(2),ng_end(2)
       do ix=ng_sta(1),ng_end(1)
         Vh_stock2(ix,iy,iz)=2.d0*Vh_stock1(ix,iy,iz)-Vh_stock2(ix,iy,iz)
       end do
       end do
       end do
     else
!$OMP parallel do
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

  elp3(516)=MPI_Wtime()
  elp3(536)=elp3(536)+elp3(516)-elp3(515)

  if(imesh_s_all==1.or.(imesh_s_all==0.and.myrank<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
    if(ilsda==0)then
      call conv_core_exc_cor
    else if(ilsda==1)then
      call Exc_cor_ns
    end if
  end if

  elp3(517)=MPI_Wtime()
  elp3(537)=elp3(537)+elp3(517)-elp3(516)

  call mpi_allgatherv_vlocal

  elp3(518)=MPI_Wtime()
  elp3(538)=elp3(538)+elp3(518)-elp3(517)

! result

  if(ikind_eext/=0.or.(ikind_eext==0.and.itt==itotNtime))then
    elp3(526)=MPI_Wtime()

    ihpsieff=0
    if(mod(itt,2)==1)then
      call Total_Energy_groupob(zpsi_out,shtpsi,1)              ! Total energy
    else
      call Total_Energy_groupob(zpsi_in,shtpsi,1)              ! Total energy
    end if
    elp3(527)=MPI_Wtime()
    elp3(540)=elp3(540)+elp3(527)-elp3(526)
    
    call subdip(rNe,1)
    elp3(528)=MPI_Wtime()
    elp3(539)=elp3(539)+elp3(528)-elp3(527)

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
!$OMP parallel do reduction( + : rbox1 )
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
  
    call MPI_allreduce(rbox_array_dip2,rbox_array2_dip2,4*num_dip2,MPI_DOUBLE_PRECISION,MPI_SUM,      &
          newworld_comm_h,ierr)

    do jj=1,num_dip2
      Dp2(1:3,itt,jj)=rbox_array2_dip2(1:3,jj)*Hgs(1:3)*Hvol-vecDs2(1:3,jj)
    end do

!------------QUADRUPOLE-start------------

    if(iflag_quadrupole==1)then
      rho_diff(:,:,:) = rho(:,:,:)-rho0(:,:,:)
      do jj=1,num_dip2
        vecR_tmp(:,:,:,:)=vecR(:,:,:,:)
        vecR_tmp(1,:,:,:)=vecR_tmp(1,:,:,:)-dip2center(jj)/Hgs(1)
        do i1=1,3
          rbox1q=0.d0
 !$OMP parallel do reduction( + : rbox1q ) private(absr2)
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
 !$OMP parallel do reduction( + : rbox1q12,rbox1q23,rbox1q31 )
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

      call MPI_allreduce(rbox_array_dip2q,rbox_array2_dip2q,9*num_dip2,MPI_DOUBLE_PRECISION,MPI_SUM,      &
            newworld_comm_h,ierr)
 
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
!$OMP parallel do reduction( + : rbox1e )
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rbox1e=rbox1e+rho(ix,iy,iz)*rto_ix(ix,jj)
        end do
        end do
        end do
        rbox_array_dip2e(jj)=rbox1e
      end do
  
      call MPI_allreduce(rbox_array_dip2e,rbox_array2_dip2e,num_dip2,MPI_DOUBLE_PRECISION,MPI_SUM,      &
            newworld_comm_h,ierr)

      do jj=1,num_dip2
        rIe2(itt,jj)=rbox_array2_dip2e(jj)*Hvol
      end do
    end if
  end if


  elp3(520)=MPI_Wtime()

  if(icalcforce==1)then
    if(mod(itt,2)==1)then
      call calc_force_c(zpsi_out)
    else
      call calc_force_c(zpsi_in)
    end if
    if(myrank==0)then
      do iatom=1,MI
        dRion(:,iatom,1)=2*dRion(:,iatom,0)-dRion(:,iatom,-1)+    &
                             rforce(:,iatom)*dt**2/(umass*Mass(Kion(iatom)))
        Rion(:,iatom)=Rion_eq(:,iatom)+dRion(:,iatom,1)
      enddo
    end if
    call MPI_Bcast(Rion,3*MI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call init_ps
    if(myrank==0)then
      dRion(:,:,-1)=dRion(:,:,0)
      dRion(:,:,0)=dRion(:,:,1)
    end if
    call MPI_Bcast(dRion(1,1,-1),3*MI*3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  end if


  if(ikind_eext==4.or.ikind_eext==14)then
    allocate(cmatbox1(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
    allocate(cmatbox2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
    
!$OMP parallel do
    do iz=lg_sta(3),lg_end(3)
    do iy=lg_sta(2),lg_end(2)
    do ix=lg_sta(1),lg_end(1)
      cmatbox1(ix,iy,iz)=0.d0
    end do
    end do
    end do
    cbox1=0.d0

    do iob=1,iobnum
      if(mod(itt,2)==1)then
!$OMP parallel do
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cmatbox1(ix,iy,iz)=zpsi_out(ix,iy,iz,iob,1)
        end do
        end do
        end do
      else
!$OMP parallel do
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cmatbox1(ix,iy,iz)=zpsi_in(ix,iy,iz,iob,1)
        end do
        end do
        end do
      end if

      call MPI_allreduce(cmatbox1,cmatbox2,      &
                         lg_num(1)*lg_num(2)*lg_num(3),      &
                         MPI_DOUBLE_COMPLEX,MPI_SUM,      &
                         newworld_comm_orbital,ierr)
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

    call MPI_allreduce(cbox1,cbox2,1,    &
                       MPI_DOUBLE_COMPLEX,MPI_SUM,      &
                       MPI_COMM_WORLD,ierr)

    cumnum=cumnum+cbox2/zi*dt

    deallocate(cmatbox1,cmatbox2) 
  end if

  if(myrank.eq.0)then
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
!$OMP parallel do 
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
  
! calculate Estatic
  if(iflag_Estatic==1)then
    call calcEstatic
  end if

! write DFT data
! WriteDensity writes rho(xyz) if he recieves "idensity", elseif he
! recieves "idiffDensity", he writes rho-rho0(xyz) to, forexample, diff30
  if (iwdenstep /= 0) then
     if (mod(itt,iwdenstep).eq.0)then
       iSCFRT=2      
!       call calcELF 
       write(fileNumber, '(i8)') itt
       rtOutFile = trim(fileTmp)//adjustl(fileNumber)
       rtDiffOutFile = trim(fileTmp2)//adjustl(fileNumber)
       call WriteDensity(rtOutFile,iwdenoption,iwdenstep,      &
                    denplane,idensum,posplane,idensity)
       call WriteDensity(rtDiffOutFile,iwdenoption,iwdenstep,      &
                    denplane,idensum,posplane,idiffDensity)
       if(iflag_Estatic==1)then
         fileTmp3="Exsta"
         rtOutFile = trim(fileTmp3)//adjustl(fileNumber)
         call WriteDensity(rtOutFile,iwdenoption,iwdenstep,      &
                    denplane,idensum,posplane,10)
         fileTmp3="Eysta"
         rtOutFile = trim(fileTmp3)//adjustl(fileNumber)
         call WriteDensity(rtOutFile,iwdenoption,iwdenstep,      &
                  denplane,idensum,posplane,11)
         fileTmp3="Ezsta"
         rtOutFile = trim(fileTmp3)//adjustl(fileNumber)
         call WriteDensity(rtOutFile,iwdenoption,iwdenstep,      &
                  denplane,idensum,posplane,12)
       end if
     endif
  end if

  elp3(521)=MPI_Wtime()
  elp3(541)=elp3(541)+elp3(521)-elp3(520)
  elp3(542)=elp3(542)+elp3(521)-elp3(511)



return

END SUBROUTINE time_evolution_step
