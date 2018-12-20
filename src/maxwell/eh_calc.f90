!
!  Copyright 2018 SALMON developers
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
!-----------------------------------------------------------------------------------------
subroutine eh_calc(grid,tmp)
  use inputoutput,          only: iobs_num_em,iobs_samp_em,directory,utime_from_au,t1_t2,t1_delay,&
                                  amplitude1,pulse_tw1,omega1,phi_cep1,epdir_re1,epdir_im1,ae_shape1,&
                                  amplitude2,pulse_tw2,omega2,phi_cep2,epdir_re2,epdir_im2,ae_shape2
  use salmon_parallel,      only: nproc_id_global,nproc_size_global,nproc_group_global
  use salmon_communication, only: comm_is_root,comm_summation
  use salmon_maxwell,       only: fdtd_grid,fdtd_tmp
  implicit none
  type(fdtd_grid)  :: grid
  type(fdtd_tmp)   :: tmp
  integer          :: iter,ii,ix,iy,iz
  real(8),parameter :: pi=3.141592653589793d0
  character(128)    :: save_name
  
  !time-iteration
  do iter=tmp%iter_sta,tmp%iter_end
    !update iter_now
    grid%iter_now=iter
    if(comm_is_root(nproc_id_global))then
      write(*,*) grid%iter_now
    end if
    
    !update drude
    if(tmp%inum_d>0) then
      call eh_update_drude
    end if
    
    !calculate linear response
    if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
      call eh_calc_lr
    end if
    
    !update e
    call eh_fd(tmp%iex_y_sta,tmp%iex_y_end,      grid%ng_sta,grid%ng_end,tmp%Nd,&
               tmp%c1_ex_y,tmp%c2_ex_y,tmp%ex_y,tmp%hz_x,tmp%hz_y,      'e','y') !ex_y
    call eh_fd(tmp%iex_z_sta,tmp%iex_z_end,      grid%ng_sta,grid%ng_end,tmp%Nd,&
               tmp%c1_ex_z,tmp%c2_ex_z,tmp%ex_z,tmp%hy_z,tmp%hy_x,      'e','z') !ex_z
    call eh_fd(tmp%iey_z_sta,tmp%iey_z_end,      grid%ng_sta,grid%ng_end,tmp%Nd,&
               tmp%c1_ey_z,tmp%c2_ey_z,tmp%ey_z,tmp%hx_y,tmp%hx_z,      'e','z') !ey_z
    call eh_fd(tmp%iey_x_sta,tmp%iey_x_end,      grid%ng_sta,grid%ng_end,tmp%Nd,&
               tmp%c1_ey_x,tmp%c2_ey_x,tmp%ey_x,tmp%hz_x,tmp%hz_y,      'e','x') !ey_x
    call eh_fd(tmp%iez_x_sta,tmp%iez_x_end,      grid%ng_sta,grid%ng_end,tmp%Nd,&
               tmp%c1_ez_x,tmp%c2_ez_x,tmp%ez_x,tmp%hy_z,tmp%hy_x,      'e','x') !ez_x
    call eh_fd(tmp%iez_y_sta,tmp%iez_y_end,      grid%ng_sta,grid%ng_end,tmp%Nd,&
               tmp%c1_ez_y,tmp%c2_ez_y,tmp%ez_y,tmp%hx_y,tmp%hx_z,      'e','y') !ez_y
    if(tmp%inc_num>0) then                                !add incident current source
      if(tmp%inc_dist1/='none') call eh_add_inc(1,amplitude1,pulse_tw1,omega1,phi_cep1,&
                                                  epdir_re1,epdir_im1,ae_shape1,tmp%inc_dist1)
      if(tmp%inc_dist2/='none') call eh_add_inc(2,amplitude2,pulse_tw2,omega2,phi_cep2,&
                                                  epdir_re2,epdir_im2,ae_shape2,tmp%inc_dist2)
    end if
    if(tmp%inum_d>0) then
      call eh_add_curr(tmp%rjx_sum_d(:,:,:),tmp%rjy_sum_d(:,:,:),tmp%rjz_sum_d(:,:,:))
    end if
    call eh_sendrecv(grid,tmp,'e')
    
    !store old h
    if( (iobs_num_em>0).and.(mod(iter,iobs_samp_em)==0) )then
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=(grid%ng_sta(3)-tmp%Nd),(grid%ng_end(3)+tmp%Nd)
      do iy=(grid%ng_sta(2)-tmp%Nd),(grid%ng_end(2)+tmp%Nd)
      do ix=(grid%ng_sta(1)-tmp%Nd),(grid%ng_end(1)+tmp%Nd)
        tmp%hx_s(ix,iy,iz)=tmp%hx_y(ix,iy,iz)+tmp%hx_z(ix,iy,iz)
        tmp%hy_s(ix,iy,iz)=tmp%hy_z(ix,iy,iz)+tmp%hy_x(ix,iy,iz)
        tmp%hz_s(ix,iy,iz)=tmp%hz_x(ix,iy,iz)+tmp%hz_y(ix,iy,iz)
      end do
      end do
      end do
!$omp end do
!$omp end parallel
    end if
    
    !update h
    call eh_fd(tmp%ihx_y_sta,tmp%ihx_y_end,      grid%ng_sta,grid%ng_end,tmp%Nd,&
               tmp%c1_hx_y,tmp%c2_hx_y,tmp%hx_y,tmp%ez_x,tmp%ez_y,      'h','y') !hx_y
    call eh_fd(tmp%ihx_z_sta,tmp%ihx_z_end,      grid%ng_sta,grid%ng_end,tmp%Nd,&
               tmp%c1_hx_z,tmp%c2_hx_z,tmp%hx_z,tmp%ey_z,tmp%ey_x,      'h','z') !hx_z
    call eh_fd(tmp%ihy_z_sta,tmp%ihy_z_end,      grid%ng_sta,grid%ng_end,tmp%Nd,&
               tmp%c1_hy_z,tmp%c2_hy_z,tmp%hy_z,tmp%ex_y,tmp%ex_z,      'h','z') !hy_z
    call eh_fd(tmp%ihy_x_sta,tmp%ihy_x_end,      grid%ng_sta,grid%ng_end,tmp%Nd,&
               tmp%c1_hy_x,tmp%c2_hy_x,tmp%hy_x,tmp%ez_x,tmp%ez_y,      'h','x') !hy_x
    call eh_fd(tmp%ihz_x_sta,tmp%ihz_x_end,      grid%ng_sta,grid%ng_end,tmp%Nd,&
               tmp%c1_hz_x,tmp%c2_hz_x,tmp%hz_x,tmp%ey_z,tmp%ey_x,      'h','x') !hz_x
    call eh_fd(tmp%ihz_y_sta,tmp%ihz_y_end,      grid%ng_sta,grid%ng_end,tmp%Nd,&
               tmp%c1_hz_y,tmp%c2_hz_y,tmp%hz_y,tmp%ex_y,tmp%ex_z,      'h','y') !hz_y
    call eh_sendrecv(grid,tmp,'h')
    
    !observation
    if( (iobs_num_em>0).and.(mod(iter,iobs_samp_em)==0) )then
      !prepare e and h for save
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=(grid%ng_sta(3)-tmp%Nd),(grid%ng_end(3)+tmp%Nd)
      do iy=(grid%ng_sta(2)-tmp%Nd),(grid%ng_end(2)+tmp%Nd)
      do ix=(grid%ng_sta(1)-tmp%Nd),(grid%ng_end(1)+tmp%Nd)
        tmp%ex_s(ix,iy,iz)=tmp%ex_y(ix,iy,iz)+tmp%ex_z(ix,iy,iz)
        tmp%ey_s(ix,iy,iz)=tmp%ey_z(ix,iy,iz)+tmp%ey_x(ix,iy,iz)
        tmp%ez_s(ix,iy,iz)=tmp%ez_x(ix,iy,iz)+tmp%ez_y(ix,iy,iz)
        tmp%hx_s(ix,iy,iz)=( tmp%hx_s(ix,iy,iz)+(tmp%hx_y(ix,iy,iz)+tmp%hx_z(ix,iy,iz)) )/2.0d0
        tmp%hy_s(ix,iy,iz)=( tmp%hy_s(ix,iy,iz)+(tmp%hy_z(ix,iy,iz)+tmp%hy_x(ix,iy,iz)) )/2.0d0
        tmp%hz_s(ix,iy,iz)=( tmp%hz_s(ix,iy,iz)+(tmp%hz_x(ix,iy,iz)+tmp%hz_y(ix,iy,iz)) )/2.0d0
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      call eh_sendrecv(grid,tmp,'s')
      
      !save data
      do ii=1,iobs_num_em
        !point
        if(tmp%iobs_po_pe(ii)==1) then
          write(save_name,*) ii
          save_name=trim(adjustl(directory))//'/obs'//trim(adjustl(save_name))//'_at_point.data'
          open(tmp%ifn,file=save_name,status='old',position='append')
          write(tmp%ifn, '(E13.5)',advance="no") dble(iter)*grid%dt*utime_from_au
          write(tmp%ifn,'(E16.6e3)',advance="no") &
                tmp%ex_s(tmp%iobs_po_id(ii,1),tmp%iobs_po_id(ii,2),tmp%iobs_po_id(ii,3))*tmp%uVperm_from_au
          write(tmp%ifn,'(E16.6e3)',advance="no") &
                tmp%ey_s(tmp%iobs_po_id(ii,1),tmp%iobs_po_id(ii,2),tmp%iobs_po_id(ii,3))*tmp%uVperm_from_au
          write(tmp%ifn,'(E16.6e3)',advance="no") &
                tmp%ez_s(tmp%iobs_po_id(ii,1),tmp%iobs_po_id(ii,2),tmp%iobs_po_id(ii,3))*tmp%uVperm_from_au
          write(tmp%ifn,'(E16.6e3)',advance="no") &
                tmp%hx_s(tmp%iobs_po_id(ii,1),tmp%iobs_po_id(ii,2),tmp%iobs_po_id(ii,3))*tmp%uAperm_from_au
          write(tmp%ifn,'(E16.6e3)',advance="no") &
                tmp%hy_s(tmp%iobs_po_id(ii,1),tmp%iobs_po_id(ii,2),tmp%iobs_po_id(ii,3))*tmp%uAperm_from_au
          write(tmp%ifn,'(E16.6e3)',advance="no") &
                tmp%hz_s(tmp%iobs_po_id(ii,1),tmp%iobs_po_id(ii,2),tmp%iobs_po_id(ii,3))*tmp%uAperm_from_au
          close(tmp%ifn)
        end if
        
        !plane
        call eh_save_plane(tmp%iobs_po_id(ii,:),tmp%iobs_pl_pe(ii,:),tmp%uVperm_from_au,&
                           grid%ng_sta,grid%ng_end,grid%lg_sta,grid%lg_end,tmp%Nd,tmp%ifn,ii,iter,tmp%ex_s,'ex')
        call eh_save_plane(tmp%iobs_po_id(ii,:),tmp%iobs_pl_pe(ii,:),tmp%uVperm_from_au,&
                           grid%ng_sta,grid%ng_end,grid%lg_sta,grid%lg_end,tmp%Nd,tmp%ifn,ii,iter,tmp%ey_s,'ey')
        call eh_save_plane(tmp%iobs_po_id(ii,:),tmp%iobs_pl_pe(ii,:),tmp%uVperm_from_au,&
                           grid%ng_sta,grid%ng_end,grid%lg_sta,grid%lg_end,tmp%Nd,tmp%ifn,ii,iter,tmp%ez_s,'ez')
        call eh_save_plane(tmp%iobs_po_id(ii,:),tmp%iobs_pl_pe(ii,:),tmp%uAperm_from_au,&
                           grid%ng_sta,grid%ng_end,grid%lg_sta,grid%lg_end,tmp%Nd,tmp%ifn,ii,iter,tmp%hx_s,'hx')
        call eh_save_plane(tmp%iobs_po_id(ii,:),tmp%iobs_pl_pe(ii,:),tmp%uAperm_from_au,&
                           grid%ng_sta,grid%ng_end,grid%lg_sta,grid%lg_end,tmp%Nd,tmp%ifn,ii,iter,tmp%hy_s,'hy')
        call eh_save_plane(tmp%iobs_po_id(ii,:),tmp%iobs_pl_pe(ii,:),tmp%uAperm_from_au,&
                           grid%ng_sta,grid%ng_end,grid%lg_sta,grid%lg_end,tmp%Nd,tmp%ifn,ii,iter,tmp%hz_s,'hz')
      end do
      
      !check maximum
      call eh_update_max
    end if
    
  end do
  
contains
  
  !=========================================================================================
  != update drude ==========================================================================
  subroutine eh_update_drude
    implicit none
    
    !initialize
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=grid%ng_sta(3),grid%ng_end(3)
    do iy=grid%ng_sta(2),grid%ng_end(2)
    do ix=grid%ng_sta(1),grid%ng_end(1)
      tmp%rjx_sum_d(ix,iy,iz)=0.0d0; tmp%rjy_sum_d(ix,iy,iz)=0.0d0; tmp%rjz_sum_d(ix,iy,iz)=0.0d0;
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    
    !update drude current
    do ii=1,tmp%inum_d
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=grid%ng_sta(3),grid%ng_end(3)
      do iy=grid%ng_sta(2),grid%ng_end(2)
      do ix=grid%ng_sta(1),grid%ng_end(1)
        tmp%rjx_d(ix,iy,iz,ii)= tmp%c1_j_d(ii)*tmp%rjx_d(ix,iy,iz,ii) &
                               +tmp%c2_j_d(ii)*( tmp%ex_y(ix,iy,iz)+tmp%ex_z(ix,iy,iz) )&
                               *dble(tmp%idx_d(ix,iy,iz,ii))
        tmp%rjy_d(ix,iy,iz,ii)= tmp%c1_j_d(ii)*tmp%rjy_d(ix,iy,iz,ii) &
                               +tmp%c2_j_d(ii)*( tmp%ey_z(ix,iy,iz)+tmp%ey_x(ix,iy,iz) )&
                               *dble(tmp%idy_d(ix,iy,iz,ii))
        tmp%rjz_d(ix,iy,iz,ii)= tmp%c1_j_d(ii)*tmp%rjz_d(ix,iy,iz,ii) &
                               +tmp%c2_j_d(ii)*( tmp%ez_x(ix,iy,iz)+tmp%ez_y(ix,iy,iz) )&
                               *dble(tmp%idz_d(ix,iy,iz,ii))
        tmp%rjx_sum_d(ix,iy,iz)=tmp%rjx_sum_d(ix,iy,iz)+tmp%wex_d(ix,iy,iz,ii)*tmp%rjx_d(ix,iy,iz,ii)
        tmp%rjy_sum_d(ix,iy,iz)=tmp%rjy_sum_d(ix,iy,iz)+tmp%wey_d(ix,iy,iz,ii)*tmp%rjy_d(ix,iy,iz,ii)
        tmp%rjz_sum_d(ix,iy,iz)=tmp%rjz_sum_d(ix,iy,iz)+tmp%wez_d(ix,iy,iz,ii)*tmp%rjz_d(ix,iy,iz,ii)
      end do
      end do
      end do
!$omp end do
!$omp end parallel
    end do
    
  end subroutine eh_update_drude
  
  !=========================================================================================
  != calculate linear response =============================================================
  subroutine eh_calc_lr
    use inputoutput,          only: iperiodic
    use salmon_parallel,      only: nproc_group_global
    use salmon_communication, only: comm_summation
    implicit none
    real(8) :: sum_lr_x,sum_lr_y,sum_lr_z
    real(8) :: sum_lr(3),sum_lr2(3)
    
    !update time
    tmp%time_lr(tmp%iter_lr)=dble(tmp%iter_lr)*grid%dt
    
    !initialize current density
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=grid%ng_sta(3),grid%ng_end(3)
    do iy=grid%ng_sta(2),grid%ng_end(2)
    do ix=grid%ng_sta(1),grid%ng_end(1)
      tmp%rjx_lr(ix,iy,iz)=0.0d0; tmp%rjy_lr(ix,iy,iz)=0.0d0; tmp%rjz_lr(ix,iy,iz)=0.0d0;
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    
    !add all current density
    if(tmp%inum_d>0) then
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=grid%ng_sta(3),grid%ng_end(3)
      do iy=grid%ng_sta(2),grid%ng_end(2)
      do ix=grid%ng_sta(1),grid%ng_end(1)
        tmp%rjx_lr(ix,iy,iz)=tmp%rjx_lr(ix,iy,iz)+tmp%rjx_sum_d(ix,iy,iz)
        tmp%rjy_lr(ix,iy,iz)=tmp%rjy_lr(ix,iy,iz)+tmp%rjy_sum_d(ix,iy,iz)
        tmp%rjz_lr(ix,iy,iz)=tmp%rjz_lr(ix,iy,iz)+tmp%rjz_sum_d(ix,iy,iz)
      end do
      end do
      end do
!$omp end do
!$omp end parallel
    end if
    
    !calculate dip or curr
    if(iperiodic==0) then
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=grid%ng_sta(3),grid%ng_end(3)
      do iy=grid%ng_sta(2),grid%ng_end(2)
      do ix=grid%ng_sta(1),grid%ng_end(1)
        tmp%px_lr(ix,iy,iz)=tmp%px_lr(ix,iy,iz)+tmp%rjx_lr(ix,iy,iz)*grid%dt
        tmp%py_lr(ix,iy,iz)=tmp%py_lr(ix,iy,iz)+tmp%rjy_lr(ix,iy,iz)*grid%dt
        tmp%pz_lr(ix,iy,iz)=tmp%pz_lr(ix,iy,iz)+tmp%rjz_lr(ix,iy,iz)*grid%dt
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      sum_lr_x=0.0d0;  sum_lr_y=0.0d0;  sum_lr_z=0.0d0;
      sum_lr(:)=0.0d0; sum_lr2(:)=0.0d0;
!$omp parallel
!$omp do private(ix,iy,iz) reduction( + : sum_lr_x,sum_lr_y,sum_lr_z )
      do iz=grid%ng_sta(3),grid%ng_end(3)
      do iy=grid%ng_sta(2),grid%ng_end(2)
      do ix=grid%ng_sta(1),grid%ng_end(1)
        sum_lr_x=sum_lr_x+tmp%px_lr(ix,iy,iz)
        sum_lr_y=sum_lr_y+tmp%py_lr(ix,iy,iz)
        sum_lr_z=sum_lr_z+tmp%pz_lr(ix,iy,iz)
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      sum_lr(1)=sum_lr_x; sum_lr(2)=sum_lr_y; sum_lr(3)=sum_lr_z;
      call comm_summation(sum_lr,sum_lr2,3,nproc_group_global)
      tmp%dip_lr(tmp%iter_lr,:)=sum_lr2(:)*grid%hgs(1)*grid%hgs(2)*grid%hgs(3)
    elseif(iperiodic==3) then
      sum_lr_x=0.0d0;  sum_lr_y=0.0d0;  sum_lr_z=0.0d0;
      sum_lr(:)=0.0d0; sum_lr2(:)=0.0d0;
!$omp parallel
!$omp do private(ix,iy,iz) reduction( + : sum_lr_x,sum_lr_y,sum_lr_z )
      do iz=grid%ng_sta(3),grid%ng_end(3)
      do iy=grid%ng_sta(2),grid%ng_end(2)
      do ix=grid%ng_sta(1),grid%ng_end(1)
        sum_lr_x=sum_lr_x+tmp%rjx_lr(ix,iy,iz)
        sum_lr_y=sum_lr_y+tmp%rjy_lr(ix,iy,iz)
        sum_lr_z=sum_lr_z+tmp%rjz_lr(ix,iy,iz)
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      sum_lr(1)=sum_lr_x; sum_lr(2)=sum_lr_y; sum_lr(3)=sum_lr_z;
      call comm_summation(sum_lr,sum_lr2,3,nproc_group_global)
      tmp%curr_lr(tmp%iter_lr,:)=sum_lr2(:)*grid%hgs(1)*grid%hgs(2)*grid%hgs(3) &
                                 /(grid%rlsize(1)*grid%rlsize(2)*grid%rlsize(3))
    end if
    
    !update time iteration
    tmp%iter_lr=tmp%iter_lr+1
    
  end subroutine eh_calc_lr
  
  !=========================================================================================
  != add incident current source ===========================================================
  subroutine eh_add_inc(iord,amp,tw,omega,cep,ep_r,ep_i,aes,typ)
    implicit none
    integer,intent(in)       :: iord
    real(8),intent(in)       :: amp,tw,omega,cep
    real(8),intent(in)       :: ep_r(3),ep_i(3)
    character(16),intent(in) :: aes,typ
    real(8)                  :: t_sta,t,theta1,theta2_r,theta2_i,alpha,beta,gamma,tf_r,tf_i
    real(8)                  :: add_inc(3)
    
    !calculate time factor and adding current
    if(iord==1) then
      t_sta=t1_delay
    elseif(iord==2) then
      t_sta=t1_delay+t1_t2
    end if
    t=(dble(iter)-0.5d0)*grid%dt-t_sta
    theta1=pi/tw*(t-0.5d0*tw)                         !for cos(theta1)**2
    alpha =pi/tw                                      !for cos(theta1)**2
    theta2_r=omega*(t-0.5d0*tw)+cep*2d0*pi            !for cos(theta2)
    theta2_i=omega*(t-0.5d0*tw)+cep*2d0*pi+3d0/2d0*pi !for cos(theta2), where this is translated to sin.
    beta=omega                                        !for cos(theta2)
    if(t>=0.0d0.and.t<=tw) then
      gamma=1.0d0
    else
      gamma=0.0d0
    end if
    if(aes=='Ecos2')then
      tf_r=cos(theta1)**2*cos(theta2_r)*gamma
      tf_i=cos(theta1)**2*cos(theta2_i)*gamma
    else if(aes=='Acos2')then
      tf_r=-(-alpha*sin(2.d0*theta1)*cos(theta2_r)   &
             -beta*cos(theta1)**2*sin(theta2_r))*gamma
      tf_i=-(-alpha*sin(2.d0*theta1)*cos(theta2_i)   &
             -beta*cos(theta1)**2*sin(theta2_i))*gamma
    end if
!    tf_r=exp(-0.5d0*(( ((dble(iter)-0.5d0)*grid%dt-10.0d0*pulse_tw1)/pulse_tw1 )**2.0d0) ) !test time factor
    add_inc(:)=amp*(tf_r*ep_r(:)+tf_i*ep_i(:))
    
    if(typ=='point') then
      if(tmp%inc_po_pe(iord)==1) then
        ix=tmp%inc_po_id(iord,1); iy=tmp%inc_po_id(iord,2); iz=tmp%inc_po_id(iord,3);
        tmp%ex_y(ix,iy,iz)=add_inc(1)/2.0d0
        tmp%ex_z(ix,iy,iz)=add_inc(1)/2.0d0
        tmp%ey_z(ix,iy,iz)=add_inc(2)/2.0d0
        tmp%ey_x(ix,iy,iz)=add_inc(2)/2.0d0
        tmp%ez_x(ix,iy,iz)=add_inc(3)/2.0d0
        tmp%ez_y(ix,iy,iz)=add_inc(3)/2.0d0
      end if
    elseif(typ=='x-line') then
      if(tmp%inc_li_pe(iord,1)==1) then
        iy=tmp%inc_po_id(iord,2); iz=tmp%inc_po_id(iord,3);
        tmp%ex_y(tmp%iex_y_sta(1):tmp%iex_y_end(1),iy,iz)=add_inc(1)/2.0d0
        tmp%ex_z(tmp%iex_z_sta(1):tmp%iex_z_end(1),iy,iz)=add_inc(1)/2.0d0
        tmp%ey_z(tmp%iey_z_sta(1):tmp%iey_z_end(1),iy,iz)=add_inc(2)/2.0d0
        tmp%ey_x(tmp%iey_x_sta(1):tmp%iey_x_end(1),iy,iz)=add_inc(2)/2.0d0
        tmp%ez_x(tmp%iez_x_sta(1):tmp%iez_x_end(1),iy,iz)=add_inc(3)/2.0d0
        tmp%ez_y(tmp%iez_y_sta(1):tmp%iez_y_end(1),iy,iz)=add_inc(3)/2.0d0
      end if
    elseif(typ=='y-line') then
      if(tmp%inc_li_pe(iord,2)==1) then
        ix=tmp%inc_po_id(iord,1); iz=tmp%inc_po_id(iord,3);
        tmp%ex_y(ix,tmp%iex_y_sta(2):tmp%iex_y_end(2),iz)=add_inc(1)/2.0d0
        tmp%ex_z(ix,tmp%iex_z_sta(2):tmp%iex_z_end(2),iz)=add_inc(1)/2.0d0
        tmp%ey_z(ix,tmp%iey_z_sta(2):tmp%iey_z_end(2),iz)=add_inc(2)/2.0d0
        tmp%ey_x(ix,tmp%iey_x_sta(2):tmp%iey_x_end(2),iz)=add_inc(2)/2.0d0
        tmp%ez_x(ix,tmp%iez_x_sta(2):tmp%iez_x_end(2),iz)=add_inc(3)/2.0d0
        tmp%ez_y(ix,tmp%iez_y_sta(2):tmp%iez_y_end(2),iz)=add_inc(3)/2.0d0
      end if
    elseif(typ=='z-line') then
      if(tmp%inc_li_pe(iord,3)==1) then
        ix=tmp%inc_po_id(iord,1); iy=tmp%inc_po_id(iord,2);
        tmp%ex_y(ix,iy,tmp%iex_y_sta(3):tmp%iex_y_end(3))=add_inc(1)/2.0d0
        tmp%ex_z(ix,iy,tmp%iex_z_sta(3):tmp%iex_z_end(3))=add_inc(1)/2.0d0
        tmp%ey_z(ix,iy,tmp%iey_z_sta(3):tmp%iey_z_end(3))=add_inc(2)/2.0d0
        tmp%ey_x(ix,iy,tmp%iey_x_sta(3):tmp%iey_x_end(3))=add_inc(2)/2.0d0
        tmp%ez_x(ix,iy,tmp%iez_x_sta(3):tmp%iez_x_end(3))=add_inc(3)/2.0d0
        tmp%ez_y(ix,iy,tmp%iez_y_sta(3):tmp%iez_y_end(3))=add_inc(3)/2.0d0
      end if
    elseif(typ=='xy-plane') then !z propagation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(tmp%inc_pl_pe(iord,1)==1) then
        iz=tmp%inc_po_id(iord,3)
!$omp parallel
!$omp do private(ix,iy)
        do iy=tmp%iex_z_sta(2),tmp%iex_z_end(2)
        do ix=tmp%iex_z_sta(1),tmp%iex_z_end(1)
          tmp%ex_z(ix,iy,iz)=tmp%ex_z(ix,iy,iz)+tmp%c2_inc_xyz(3)*add_inc(1)
        end do
        end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iy)
        do iy=tmp%iey_z_sta(2),tmp%iey_z_end(2)
        do ix=tmp%iey_z_sta(1),tmp%iey_z_end(1)
          tmp%ey_z(ix,iy,iz)=tmp%ey_z(ix,iy,iz)+tmp%c2_inc_xyz(3)*add_inc(2)
        end do
        end do
!$omp end do
!$omp end parallel
      end if
    elseif(typ=='yz-plane') then !x propagation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(tmp%inc_pl_pe(iord,2)==1) then
        ix=tmp%inc_po_id(iord,1)
!$omp parallel
!$omp do private(iy,iz)
        do iz=tmp%iey_x_sta(3),tmp%iey_x_end(3)
        do iy=tmp%iey_x_sta(2),tmp%iey_x_end(2)
          tmp%ey_x(ix,iy,iz)=tmp%ey_x(ix,iy,iz)+tmp%c2_inc_xyz(1)*add_inc(2)
        end do
        end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(iy,iz)
        do iz=tmp%iez_x_sta(3),tmp%iez_x_end(3)
        do iy=tmp%iez_x_sta(2),tmp%iez_x_end(2)
          tmp%ez_x(ix,iy,iz)=tmp%ez_x(ix,iy,iz)+tmp%c2_inc_xyz(1)*add_inc(3)
        end do
        end do
!$omp end do
!$omp end parallel
      end if
    elseif(typ=='xz-plane') then !y propagation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(tmp%inc_pl_pe(iord,3)==1) then
        iy=tmp%inc_po_id(iord,2)
!$omp parallel
!$omp do private(ix,iz)
        do iz=tmp%iex_y_sta(3),tmp%iex_y_end(3)
        do ix=tmp%iex_y_sta(1),tmp%iex_y_end(1)
          tmp%ex_y(ix,iy,iz)=tmp%ex_y(ix,iy,iz)+tmp%c2_inc_xyz(2)*add_inc(1)
        end do
        end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iz)
        do iz=tmp%iez_y_sta(3),tmp%iez_y_end(3)
        do ix=tmp%iez_y_sta(1),tmp%iez_y_end(1)
          tmp%ez_y(ix,iy,iz)=tmp%ez_y(ix,iy,iz)+tmp%c2_inc_xyz(2)*add_inc(3)
        end do
        end do
!$omp end do
!$omp end parallel
      end if
    end if
    
  end subroutine eh_add_inc
  
  !=========================================================================================
  != add current ===========================================================================
  subroutine eh_add_curr(rjx,rjy,rjz)
    implicit none
    real(8),intent(in) :: rjx(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                              grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                              grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
                          rjy(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                              grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                              grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
                          rjz(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                              grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                              grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd)
    
    !ex
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=tmp%iex_y_sta(3),tmp%iex_y_end(3)
    do iy=tmp%iex_y_sta(2),tmp%iex_y_end(2)
    do ix=tmp%iex_y_sta(1),tmp%iex_y_end(1)
      tmp%ex_y(ix,iy,iz)=tmp%ex_y(ix,iy,iz)+tmp%c2_jx(ix,iy,iz)*rjx(ix,iy,iz)/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=tmp%iex_z_sta(3),tmp%iex_z_end(3)
    do iy=tmp%iex_z_sta(2),tmp%iex_z_end(2)
    do ix=tmp%iex_z_sta(1),tmp%iex_z_end(1)
      tmp%ex_z(ix,iy,iz)=tmp%ex_z(ix,iy,iz)+tmp%c2_jx(ix,iy,iz)*rjx(ix,iy,iz)/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    
    !ey
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=tmp%iey_z_sta(3),tmp%iey_z_end(3)
    do iy=tmp%iey_z_sta(2),tmp%iey_z_end(2)
    do ix=tmp%iey_z_sta(1),tmp%iey_z_end(1)
      tmp%ey_z(ix,iy,iz)=tmp%ey_z(ix,iy,iz)+tmp%c2_jy(ix,iy,iz)*rjy(ix,iy,iz)/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=tmp%iey_x_sta(3),tmp%iey_x_end(3)
    do iy=tmp%iey_x_sta(2),tmp%iey_x_end(2)
    do ix=tmp%iey_x_sta(1),tmp%iey_x_end(1)
      tmp%ey_x(ix,iy,iz)=tmp%ey_x(ix,iy,iz)+tmp%c2_jy(ix,iy,iz)*rjy(ix,iy,iz)/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    
    !ez
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=tmp%iez_x_sta(3),tmp%iez_x_end(3)
    do iy=tmp%iez_x_sta(2),tmp%iez_x_end(2)
    do ix=tmp%iez_x_sta(1),tmp%iez_x_end(1)
      tmp%ez_x(ix,iy,iz)=tmp%ez_x(ix,iy,iz)+tmp%c2_jz(ix,iy,iz)*rjz(ix,iy,iz)/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=tmp%iez_y_sta(3),tmp%iez_y_end(3)
    do iy=tmp%iez_y_sta(2),tmp%iez_y_end(2)
    do ix=tmp%iez_y_sta(1),tmp%iez_y_end(1)
      tmp%ez_y(ix,iy,iz)=tmp%ez_y(ix,iy,iz)+tmp%c2_jz(ix,iy,iz)*rjz(ix,iy,iz)/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    
  end subroutine eh_add_curr
  
  !=========================================================================================
  != check and update maximum of e and h ===================================================
  subroutine eh_update_max
    implicit none
    real(8) :: fe(grid%ng_sta(1):grid%ng_end(1),grid%ng_sta(2):grid%ng_end(2),grid%ng_sta(3):grid%ng_end(3)),&
               fh(grid%ng_sta(1):grid%ng_end(1),grid%ng_sta(2):grid%ng_end(2),grid%ng_sta(3):grid%ng_end(3))
    real(8) :: e_max_tmp(0:nproc_size_global-1), h_max_tmp(0:nproc_size_global-1),&
               e_max_tmp2(0:nproc_size_global-1),h_max_tmp2(0:nproc_size_global-1)
    
    e_max_tmp(:)=0.0d0; h_max_tmp(:)=0.0d0;
    do iz=grid%ng_sta(3),grid%ng_end(3)
    do iy=grid%ng_sta(2),grid%ng_end(2)
    do ix=grid%ng_sta(1),grid%ng_end(1)
      fe(ix,iy,iz)=sqrt( tmp%ex_s(ix,iy,iz)**2.0d0 + tmp%ey_s(ix,iy,iz)**2.0d0 + tmp%ez_s(ix,iy,iz)**2.0d0 )
      fh(ix,iy,iz)=sqrt( tmp%hx_s(ix,iy,iz)**2.0d0 + tmp%hy_s(ix,iy,iz)**2.0d0 + tmp%hz_s(ix,iy,iz)**2.0d0 )
      if(e_max_tmp(nproc_id_global)<fe(ix,iy,iz)) e_max_tmp(nproc_id_global)=fe(ix,iy,iz)
      if(h_max_tmp(nproc_id_global)<fh(ix,iy,iz)) h_max_tmp(nproc_id_global)=fh(ix,iy,iz)
    end do
    end do
    end do
    call comm_summation(e_max_tmp,e_max_tmp2,nproc_size_global,nproc_group_global)
    call comm_summation(h_max_tmp,h_max_tmp2,nproc_size_global,nproc_group_global)
    e_max_tmp2(:)=e_max_tmp2(:)*tmp%uVperm_from_au
    h_max_tmp2(:)=h_max_tmp2(:)*tmp%uAperm_from_au
    if(tmp%e_max<maxval(e_max_tmp2(:))) tmp%e_max=maxval(e_max_tmp2(:))
    if(tmp%h_max<maxval(h_max_tmp2(:))) tmp%h_max=maxval(h_max_tmp2(:))
    
  end subroutine eh_update_max
  
end subroutine eh_calc

!=========================================================================================
!= calculate finite difference in eh =====================================================
subroutine eh_fd(ista,iend,ng_sta,ng_end,Nd,c1,c2,f1,f2,f3,var,dir)
  implicit none
  integer,intent(in)      :: ista(3),iend(3),ng_sta(3),ng_end(3)
  integer,intent(in)      :: Nd
  real(8),intent(in)      :: c1(ng_sta(1)-Nd:ng_end(1)+Nd,&
                                ng_sta(2)-Nd:ng_end(2)+Nd,&
                                ng_sta(3)-Nd:ng_end(3)+Nd),&
                             c2(ng_sta(1)-Nd:ng_end(1)+Nd,&
                                ng_sta(2)-Nd:ng_end(2)+Nd,&
                                ng_sta(3)-Nd:ng_end(3)+Nd)
  real(8),intent(inout)   :: f1(ng_sta(1)-Nd:ng_end(1)+Nd,&
                                ng_sta(2)-Nd:ng_end(2)+Nd,&
                                ng_sta(3)-Nd:ng_end(3)+Nd)
  real(8),intent(in)      :: f2(ng_sta(1)-Nd:ng_end(1)+Nd,&
                                ng_sta(2)-Nd:ng_end(2)+Nd,&
                                ng_sta(3)-Nd:ng_end(3)+Nd),&
                             f3(ng_sta(1)-Nd:ng_end(1)+Nd,&
                                ng_sta(2)-Nd:ng_end(2)+Nd,&
                                ng_sta(3)-Nd:ng_end(3)+Nd)
  character(1),intent(in) :: var,dir
  integer :: ix,iy,iz
  
  if(var=='e') then
    if(dir=='x') then
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=ista(3),iend(3)
      do iy=ista(2),iend(2)
      do ix=ista(1),iend(1)
        f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                     +c2(ix,iy,iz)*( &
                     (f2(ix,iy,iz)+f3(ix,iy,iz))-(f2(ix-1,iy,iz)+f3(ix-1,iy,iz)) )
      end do
      end do
      end do
!$omp end do
!$omp end parallel
    elseif(dir=='y') then
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=ista(3),iend(3)
      do iy=ista(2),iend(2)
      do ix=ista(1),iend(1)
        f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                     +c2(ix,iy,iz)*( &
                     (f2(ix,iy,iz)+f3(ix,iy,iz))-(f2(ix,iy-1,iz)+f3(ix,iy-1,iz)) )
      end do
      end do
      end do
!$omp end do
!$omp end parallel
    elseif(dir=='z') then
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=ista(3),iend(3)
      do iy=ista(2),iend(2)
      do ix=ista(1),iend(1)
        f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                     +c2(ix,iy,iz)*( &
                     (f2(ix,iy,iz)+f3(ix,iy,iz))-(f2(ix,iy,iz-1)+f3(ix,iy,iz-1)) )
      end do
      end do
      end do
!$omp end do
!$omp end parallel
    end if
  elseif(var=='h') then
    if(dir=='x') then
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=ista(3),iend(3)
      do iy=ista(2),iend(2)
      do ix=ista(1),iend(1)
        f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                     +c2(ix,iy,iz)*( &
                     (f2(ix+1,iy,iz)+f3(ix+1,iy,iz))-(f2(ix,iy,iz)+f3(ix,iy,iz)) )
      end do
      end do
      end do
!$omp end do
!$omp end parallel
    elseif(dir=='y') then
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=ista(3),iend(3)
      do iy=ista(2),iend(2)
      do ix=ista(1),iend(1)
        f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                     +c2(ix,iy,iz)*( &
                     (f2(ix,iy+1,iz)+f3(ix,iy+1,iz))-(f2(ix,iy,iz)+f3(ix,iy,iz)) )
      end do
      end do
      end do
!$omp end do
!$omp end parallel
    elseif(dir=='z') then
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=ista(3),iend(3)
      do iy=ista(2),iend(2)
      do ix=ista(1),iend(1)
        f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                     +c2(ix,iy,iz)*( &
                     (f2(ix,iy,iz+1)+f3(ix,iy,iz+1))-(f2(ix,iy,iz)+f3(ix,iy,iz)) )
      end do
      end do
      end do
!$omp end do
!$omp end parallel
    end if
  end if
  
end subroutine eh_fd

!=========================================================================================
!= save plane data =======================================================================
subroutine eh_save_plane(id,ipl,conv,ng_sta,ng_end,lg_sta,lg_end,Nd,ifn,iobs,iter,f,var)
  use inputoutput,          only: directory
  use salmon_parallel,      only: nproc_id_global,nproc_group_global
  use salmon_communication, only: comm_is_root,comm_summation
  implicit none
  integer,intent(in)      :: id(3),ipl(3)
  real(8),intent(in)      :: conv
  integer,intent(in)      :: ng_sta(3),ng_end(3),lg_sta(3),lg_end(3)
  integer,intent(in)      :: Nd,ifn,iobs,iter
  real(8),intent(in)      :: f(ng_sta(1)-Nd:ng_end(1)+Nd,&
                               ng_sta(2)-Nd:ng_end(2)+Nd,&
                               ng_sta(3)-Nd:ng_end(3)+Nd)
  character(2),intent(in) :: var
  real(8),allocatable     :: save_pl(:,:),save_pl2(:,:)
  integer          :: ii,inum,i1,i1s,i2,i2s
  character(2)     :: plane_name
  character(128)   :: iobs_name,iter_name,save_name
  
  do ii=1,3
    !allocate
    if(ii==1) then     !xy
      i1s=1; i2s=2; plane_name='xy';
    elseif(ii==2) then !yz
      i1s=2; i2s=3; plane_name='yz';
    elseif(ii==3) then !xz
      i1s=1; i2s=3; plane_name='xz';
    end if
    allocate(save_pl(lg_sta(i1s):lg_end(i1s),lg_sta(i2s):lg_end(i2s)),&
             save_pl2(lg_sta(i1s):lg_end(i1s),lg_sta(i2s):lg_end(i2s)))
    save_pl(:,:)=0.0d0; save_pl2(:,:)=0.0d0
    inum=(lg_end(i1s)-lg_sta(i1s)+1)*(lg_end(i2s)-lg_sta(i2s)+1)
    
    !prepare save data
    if(ipl(ii)==1) then
      do i2=ng_sta(i2s),ng_end(i2s)
      do i1=ng_sta(i1s),ng_end(i1s)
        if(ii==1) then     !xy
          save_pl(i1,i2)=f(i1,i2,id(3))
        elseif(ii==2) then !yz
          save_pl(i1,i2)=f(id(1),i1,i2)
        elseif(ii==3) then !xz
          save_pl(i1,i2)=f(i1,id(2),i2)
        end if
      end do
      end do
    end if
    call comm_summation(save_pl,save_pl2,inum,nproc_group_global)
    
    !make save data
    if(comm_is_root(nproc_id_global)) then
      write(iobs_name,*) iobs
      write(iter_name,*) iter
      save_name=trim(adjustl(directory))//'/obs'//trim(adjustl(iobs_name))//'_'//var//&
                '_'//plane_name//'_'//trim(adjustl(iter_name))//'.data'
      open(ifn,file=save_name)
      do i2=lg_sta(i2s),lg_end(i2s)
      do i1=lg_sta(i1s),lg_end(i1s)
        write(ifn,'(I8,I8,E16.6e3)') i1,i2,save_pl2(i1,i2)*conv
      end do
      end do
      close(ifn)
    end if
    
    !deallocate
    deallocate(save_pl,save_pl2)
  end do
  
end subroutine eh_save_plane

!=========================================================================================
!= send and receive eh ===================================================================
!= (This routine is temporary) ===========================================================
!= (With unifying ARTED and GCEED, this routine will be removed) =========================
subroutine eh_sendrecv(grid,tmp,var)
  use scf_data,       only: iwk_size
  use sendrecvh_sub,  only: sendrecvh
  use salmon_maxwell, only: fdtd_grid,fdtd_tmp
  implicit none
  type(fdtd_grid)         :: grid
  type(fdtd_tmp)          :: tmp
  character(1),intent(in) :: var
  integer                 :: ix,iy,iz
  real(8),allocatable     :: f1(:,:,:),f2(:,:,:),f3(:,:,:)
  
  iwk_size=tmp%iwk_size_eh
  if(var=='e') then
    call sendrecvh(tmp%ex_y)
    call sendrecvh(tmp%ex_z)
    call sendrecvh(tmp%ey_z)
    call sendrecvh(tmp%ey_x)
    call sendrecvh(tmp%ez_x)
    call sendrecvh(tmp%ez_y)
  elseif(var=='h') then
    call sendrecvh(tmp%hx_y)
    call sendrecvh(tmp%hx_z)
    call sendrecvh(tmp%hy_z)
    call sendrecvh(tmp%hy_x)
    call sendrecvh(tmp%hz_x)
    call sendrecvh(tmp%hz_y)
  elseif(var=='r') then
    call sendrecvh(tmp%rmedia)
  elseif(var=='s') then
    !allocate temporary variable
    allocate(f1(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
             f2(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
             f3(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    f1(:,:,:)=0.0d0; f2(:,:,:)=0.0d0; f3(:,:,:)=0.0d0;
    
    !spatially adjust e for save
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=grid%ng_sta(3),grid%ng_end(3)
    do iy=grid%ng_sta(2),grid%ng_end(2)
    do ix=grid%ng_sta(1),grid%ng_end(1)
      f1(ix,iy,iz)=( tmp%ex_s(ix,iy,iz)+tmp%ex_s(ix-1,iy,iz) )/2.0d0
      f2(ix,iy,iz)=( tmp%ey_s(ix,iy,iz)+tmp%ey_s(ix,iy-1,iz) )/2.0d0
      f3(ix,iy,iz)=( tmp%ez_s(ix,iy,iz)+tmp%ez_s(ix,iy,iz-1) )/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    tmp%ex_s(:,:,:)=f1(:,:,:); tmp%ey_s(:,:,:)=f2(:,:,:); tmp%ez_s(:,:,:)=f3(:,:,:);
    f1(:,:,:)=0.0d0; f2(:,:,:)=0.0d0; f3(:,:,:)=0.0d0;
    
    !spatially adjust h for save
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=grid%ng_sta(3),grid%ng_end(3)
    do iy=grid%ng_sta(2),grid%ng_end(2)
    do ix=grid%ng_sta(1),grid%ng_end(1)
      f1(ix,iy,iz)=( tmp%hx_s(ix,iy,iz)+tmp%hx_s(ix,iy-1,iz) )/2.0d0
      f2(ix,iy,iz)=( tmp%hy_s(ix,iy,iz)+tmp%hy_s(ix,iy,iz-1) )/2.0d0
      f3(ix,iy,iz)=( tmp%hz_s(ix,iy,iz)+tmp%hz_s(ix-1,iy,iz) )/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    tmp%hx_s(:,:,:)=f1(:,:,:); tmp%hy_s(:,:,:)=f2(:,:,:); tmp%hz_s(:,:,:)=f3(:,:,:);
    f1(:,:,:)=0.0d0; f2(:,:,:)=0.0d0; f3(:,:,:)=0.0d0;
    call sendrecvh(tmp%hx_s)
    call sendrecvh(tmp%hy_s)
    call sendrecvh(tmp%hz_s)
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=grid%ng_sta(3),grid%ng_end(3)
    do iy=grid%ng_sta(2),grid%ng_end(2)
    do ix=grid%ng_sta(1),grid%ng_end(1)
      f1(ix,iy,iz)=( tmp%hx_s(ix,iy,iz)+tmp%hx_s(ix,iy,iz-1) )/2.0d0
      f2(ix,iy,iz)=( tmp%hy_s(ix,iy,iz)+tmp%hy_s(ix-1,iy,iz) )/2.0d0
      f3(ix,iy,iz)=( tmp%hz_s(ix,iy,iz)+tmp%hz_s(ix,iy-1,iz) )/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    tmp%hx_s(:,:,:)=f1(:,:,:); tmp%hy_s(:,:,:)=f2(:,:,:); tmp%hz_s(:,:,:)=f3(:,:,:);
    
    !deallocate temporary variable
    deallocate(f1,f2,f3)
  end if
  iwk_size=2
  
end subroutine eh_sendrecv
