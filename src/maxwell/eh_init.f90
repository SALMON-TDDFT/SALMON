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
subroutine eh_init(grid,tmp)
  use inputoutput,          only: nt_em,al_em,dl_em,dt_em,utime_from_au,ulength_from_au,uenergy_from_au,unit_system,&
                                  uenergy_to_au,ulength_to_au,ucharge_to_au,iperiodic,directory,&
                                  imedia_num,shape_file,epsilon,rmu,sigma,type_media,omega_p_d,gamma_d,&
                                  iobs_num_em,obs_loc_em,wave_input,&
                                  source_loc1,ek_dir1,epdir_re1,epdir_im1,ae_shape1,&
                                  phi_cep1,rlaser_int_wcm2_1,amplitude1,&
                                  source_loc2,ek_dir2,epdir_re2,epdir_im2,ae_shape2,&
                                  phi_cep2,rlaser_int_wcm2_2,amplitude2
  use salmon_parallel,      only: nproc_id_global, nproc_group_global
  use salmon_communication, only: comm_is_root, comm_bcast
  use salmon_maxwell,       only:fdtd_grid,fdtd_tmp
  implicit none
  type(fdtd_grid)     :: grid
  type(fdtd_tmp)      :: tmp
  integer             :: ii,ix,iy,iz,icount,icount_d
  real(8),parameter   :: pi=3.141592653589793d0
  real(8)             :: dt_cfl,diff_cep
  character(1)        :: dir
  character(128)      :: save_name
  
  !set initial parameter and value
  tmp%c_0        = 1.370359991378353d2
  tmp%Nd         = 4
  tmp%iter_sta   = 1
  tmp%iter_end   = nt_em
  grid%rlsize(:) = al_em(:)
  grid%hgs(:)    = dl_em(:)
  tmp%ifn         = 600
  tmp%ipml_l     = 8
  tmp%pml_m      = 4.0d0
  tmp%pml_r      = 1.0d-7
  if(iperiodic==0) then
    grid%i_bc(:,:)=1 !pml
  elseif(iperiodic==3) then
    grid%i_bc(:,:)=0
    stop 'invalid iperiodic. For theory = Maxwell, iperiodic is only allowed by 0.'
  end if
  select case(unit_system)
  case('au','a.u.')
    tmp%uVperm_from_au=1.0d0
    tmp%uAperm_from_au=1.0d0
  case('A_eV_fs')
    !see amplitude1 or amplitude2 in src/io/iunputoutput.f90
    tmp%uVperm_from_au=1/(uenergy_to_au/ulength_to_au/ucharge_to_au)
    tmp%uAperm_from_au=tmp%uVperm_from_au
  end select
  
  !prepare GCEED(set mpi condition, gird, coordinate, and sendrecv environment)
  call eh_prep_GCEED(grid,tmp)
  
  !set and check dt
  dt_cfl=1.0d0/( &
         tmp%c_0*sqrt( (1.0d0/grid%hgs(1))**2.0d0+(1.0d0/grid%hgs(2))**2.0d0+(1.0d0/grid%hgs(3))**2.0d0 ) &
         )
  if(dt_em==0.0d0) then
    grid%dt=dt_cfl*0.99d0
    if(comm_is_root(nproc_id_global)) then
      write(*,*) "**************************"
      write(*,*) "From CFL condition, dt_em is determined by", grid%dt*utime_from_au
      write(*,*) "in the unit system, ",trim(unit_system),"."
      write(*,*) "**************************"
    end if
  elseif(dt_em>=dt_cfl) then
    if(comm_is_root(nproc_id_global)) then
      write(*,*) "**************************"
      write(*,*) "To sufficient CFL condition, dt_em must be set smaller than", dt_cfl*utime_from_au
      write(*,*) "in the unit system, ",trim(unit_system),"."
      write(*,*) "**************************"
    end if
    stop
  else
    grid%dt=dt_em
    write(*,*) "**************************"
    write(*,*) "dt_em =", grid%dt*utime_from_au
    write(*,*) "in the unit system, ",trim(unit_system),"."
    write(*,*) "**************************"
  end if
  call comm_bcast(grid%dt,nproc_group_global)
  
  !allocate
  call eh_allocate
  allocate(tmp%c2_jx(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                     grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                     grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
           tmp%c2_jy(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                     grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                     grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
           tmp%c2_jz(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                     grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                     grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd) )
  tmp%c2_jx(:,:,:)=0.0d0; tmp%c2_jy(:,:,:)=0.0d0; tmp%c2_jz(:,:,:)=0.0d0;
  
  !input fdtd shape
  allocate(grid%imedia(grid%ng_sta(1):grid%ng_end(1),&
                          grid%ng_sta(2):grid%ng_end(2),&
                          grid%ng_sta(3):grid%ng_end(3)))
  grid%imedia(:,:,:)=0
  if(imedia_num>0) then
    !check file format and input shape file
    if(comm_is_root(nproc_id_global)) write(*,*)
    if(comm_is_root(nproc_id_global)) write(*,*) "**************************"
    if(index(shape_file,".cube", back=.true.)/=0) then
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "shape file is inputed by .cube format."
      end if
      call eh_input_shape(tmp%ifn,grid%ng_sta,grid%ng_end,grid%lg_sta,grid%lg_end,grid%imedia,'cu')
    elseif(index(shape_file,".mp", back=.true.)/=0) then
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "shape file is inputed by .mp format."
        write(*,*) "This version works for only .cube format.."
      end if
      stop
    else
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "shape file must be .cube or .mp formats."
      end if
      stop
    end if
    if(comm_is_root(nproc_id_global)) write(*,*) "**************************"
  end if
  
  !prepare Drude
  tmp%inum_d=0
  do ii=0,imedia_num
    select case(type_media(ii))
    case('DRUDE','Drude','drude','D','d')
      tmp%inum_d=tmp%inum_d+1
    end select
  end do
  if(tmp%inum_d>0) then
    !set counter
    icount_d=1;
    
    !allocate drude variable
    allocate(tmp%idx_d(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                       grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                       grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd,tmp%inum_d),&
             tmp%idy_d(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                       grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                       grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd,tmp%inum_d),&
             tmp%idz_d(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                       grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                       grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd,tmp%inum_d) )
    tmp%idx_d(:,:,:,:)=0; tmp%idy_d(:,:,:,:)=0; tmp%idz_d(:,:,:,:)=0;
    allocate( tmp%rjx_d(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                       grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                       grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd,tmp%inum_d),&
              tmp%rjy_d(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                       grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                       grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd,tmp%inum_d),&
              tmp%rjz_d(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                       grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                       grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd,tmp%inum_d) )
    tmp%rjx_d(:,:,:,:)=0.0d0; tmp%rjy_d(:,:,:,:)=0.0d0; tmp%rjz_d(:,:,:,:)=0.0d0;
    allocate(tmp%c1_j_d(tmp%inum_d),tmp%c2_j_d(tmp%inum_d))
    tmp%c1_j_d(:)=0.0d0; tmp%c2_j_d(:)=0.0d0;
  end if
  
  !set fdtd coeffient and write media information
  allocate(tmp%rep(0:imedia_num),tmp%rmu(0:imedia_num),tmp%sig(0:imedia_num))
  tmp%rep(:)=1.0d0; tmp%rmu(:)=1.0d0; tmp%sig(:)=0.0d0;
  do ii=0,imedia_num
    call eh_coeff
  end do
  deallocate(grid%imedia)
  if(comm_is_root(nproc_id_global)) then
    write(*,*)
    write(*,*) "**************************"
    write(*,'(A,I3)')          ' imedia_num = ',imedia_num
    do ii=0,imedia_num
      write(*,*) "=========================="
      write(*,'(A,I3,A)')      ' id ',ii, ':'
      select case(type_media(ii))
      case ('VACUUM','Vacuum','vacuum')
        if(epsilon(ii)/=1d0 .or. rmu(ii)/=1d0 .or. sigma(ii)/=0d0) then
          write(*,'(A)'  )     ' type_media =  constant media'
        else
          write(*,'(A,A)')       ' type_media =  ', trim(type_media(ii))
        end if
      case('DRUDE','Drude','drude','D','d')
        write(*,'(A,A)')       ' type_media =  ', trim(type_media(ii))
        write(*,'(A,ES12.5)')  ' omega_p_d  = ', omega_p_d(ii)*uenergy_from_au
        write(*,'(A,ES12.5)')  ' gamma_d    = ', gamma_d(ii)*uenergy_from_au
      case default
        write(*,'(A,A)')       ' type_media =  ', trim(type_media(ii))
      end select
      write(*,'(A,ES12.5)')    ' epsilon    = ', epsilon(ii)
      write(*,'(A,ES12.5)')    ' rmu        = ', rmu(ii)
      write(*,'(A,ES12.5)'   ) ' sigma      = ', sigma(ii)
    end do
    write(*,*) "**************************"
  end if
  
  !set calculation area
  tmp%iex_y_sta(:)=grid%ng_sta(:); tmp%iex_y_end(:)=grid%ng_end(:);
  tmp%iex_z_sta(:)=grid%ng_sta(:); tmp%iex_z_end(:)=grid%ng_end(:);
  tmp%iey_z_sta(:)=grid%ng_sta(:); tmp%iey_z_end(:)=grid%ng_end(:);
  tmp%iey_x_sta(:)=grid%ng_sta(:); tmp%iey_x_end(:)=grid%ng_end(:);
  tmp%iez_x_sta(:)=grid%ng_sta(:); tmp%iez_x_end(:)=grid%ng_end(:);
  tmp%iez_y_sta(:)=grid%ng_sta(:); tmp%iez_y_end(:)=grid%ng_end(:);
  tmp%ihx_y_sta(:)=grid%ng_sta(:); tmp%ihx_y_end(:)=grid%ng_end(:);
  tmp%ihx_z_sta(:)=grid%ng_sta(:); tmp%ihx_z_end(:)=grid%ng_end(:);
  tmp%ihy_z_sta(:)=grid%ng_sta(:); tmp%ihy_z_end(:)=grid%ng_end(:);
  tmp%ihy_x_sta(:)=grid%ng_sta(:); tmp%ihy_x_end(:)=grid%ng_end(:);
  tmp%ihz_x_sta(:)=grid%ng_sta(:); tmp%ihz_x_end(:)=grid%ng_end(:);
  tmp%ihz_y_sta(:)=grid%ng_sta(:); tmp%ihz_y_end(:)=grid%ng_end(:);
  if((grid%i_bc(1,1)==1).and.(grid%ng_sta(1)==grid%lg_sta(1))) then !x, bottom
    tmp%iey_x_sta(1)=grid%ng_sta(1)+1; tmp%iez_x_sta(1)=grid%ng_sta(1)+1;
  end if
  if((grid%i_bc(1,2)==1).and.(grid%ng_end(1)==grid%lg_end(1))) then !x, top
    tmp%iex_y_end(1)=grid%ng_end(1)-1; tmp%iex_z_end(1)=grid%ng_end(1)-1;
    tmp%iey_x_end(1)=grid%ng_end(1)-1; tmp%iez_x_end(1)=grid%ng_end(1)-1;
    tmp%ihy_z_end(1)=grid%ng_end(1)-1; tmp%ihy_x_end(1)=grid%ng_end(1)-1;
    tmp%ihz_x_end(1)=grid%ng_end(1)-1; tmp%ihz_y_end(1)=grid%ng_end(1)-1;
  end if
  if((grid%i_bc(2,1)==1).and.(grid%ng_sta(2)==grid%lg_sta(2))) then !y, bottom
    tmp%iex_y_sta(2)=grid%ng_sta(2)+1; tmp%iez_y_sta(2)=grid%ng_sta(2)+1;
  end if
  if((grid%i_bc(2,2)==1).and.(grid%ng_end(2)==grid%lg_end(2))) then !y, top
    tmp%iex_y_end(2)=grid%ng_end(2)-1; tmp%iey_z_end(2)=grid%ng_end(2)-1;
    tmp%iey_x_end(2)=grid%ng_end(2)-1; tmp%iez_y_end(2)=grid%ng_end(2)-1;
    tmp%ihx_y_end(2)=grid%ng_end(2)-1; tmp%ihx_z_end(2)=grid%ng_end(2)-1;
    tmp%ihz_x_end(2)=grid%ng_end(2)-1; tmp%ihz_y_end(2)=grid%ng_end(2)-1;
  end if
  if((grid%i_bc(3,1)==1).and.(grid%ng_sta(3)==grid%lg_sta(3))) then !z, bottom
    tmp%iex_z_sta(3)=grid%ng_sta(3)+1; tmp%iey_z_sta(3)=grid%ng_sta(3)+1;
  end if
  if((grid%i_bc(3,2)==1).and.(grid%ng_end(3)==grid%lg_end(3))) then !z, top
    tmp%iex_z_end(3)=grid%ng_end(3)-1; tmp%iey_z_end(3)=grid%ng_end(3)-1;
    tmp%iez_x_end(3)=grid%ng_end(3)-1; tmp%iez_y_end(3)=grid%ng_end(3)-1;
    tmp%ihx_y_end(3)=grid%ng_end(3)-1; tmp%ihx_z_end(3)=grid%ng_end(3)-1;
    tmp%ihy_z_end(3)=grid%ng_end(3)-1; tmp%ihy_x_end(3)=grid%ng_end(3)-1;
  end if
  
  !set pml
  call eh_set_pml(1,tmp%c1_ey_x,tmp%c2_ey_x,tmp%c1_ez_x,tmp%c2_ez_x,&
                    tmp%c1_hy_x,tmp%c2_hy_x,tmp%c1_hz_x,tmp%c2_hz_x) !x direction
  call eh_set_pml(2,tmp%c1_ez_y,tmp%c2_ez_y,tmp%c1_ex_y,tmp%c2_ex_y,&
                    tmp%c1_hz_y,tmp%c2_hz_y,tmp%c1_hx_y,tmp%c2_hx_y) !y direction
  call eh_set_pml(3,tmp%c1_ex_z,tmp%c2_ex_z,tmp%c1_ey_z,tmp%c2_ey_z,&
                    tmp%c1_hx_z,tmp%c2_hx_z,tmp%c1_hy_z,tmp%c2_hy_z) !z direction
  if(maxval(grid%i_bc(:,:))>0) then
    if(comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "**************************"
      do ii=1,3
        if(ii==1) then
          dir='x'
        elseif(ii==2) then
          dir='y'
        elseif(ii==3) then
          dir='z'
        end if
        if(grid%i_bc(ii,1)==1) write(*,'(A,A,A,ES12.5,A,ES12.5,A)') &
                               ' PML has been set for ',dir,'-direction: ',&
                               tmp%coo(grid%lg_sta(ii),ii)*ulength_from_au,' to ',&
                               tmp%coo(grid%lg_sta(ii)+tmp%ipml_l,ii)*ulength_from_au,'.'
        if(grid%i_bc(ii,2)==1) write(*,'(A,A,A,ES12.5,A,ES12.5,A)') &
                               ' PML has been set for ',dir,'-direction: ',&
                               tmp%coo(grid%lg_end(ii)-tmp%ipml_l,ii)*ulength_from_au,' to ',&
                               tmp%coo(grid%lg_end(ii),ii)*ulength_from_au,'.'
      end do
      write(*,*) "**************************"
    end if
  end if
  
  !prepare observation
  if(iobs_num_em>0) then
    !set initial
    allocate(tmp%ex_s(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
             tmp%ey_s(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
             tmp%ez_s(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
             tmp%hx_s(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
             tmp%hy_s(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
             tmp%hz_s(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    tmp%ex_s(:,:,:)=0; tmp%ey_s(:,:,:)=0; tmp%ez_s(:,:,:)=0; 
    tmp%hx_s(:,:,:)=0; tmp%hy_s(:,:,:)=0; tmp%hz_s(:,:,:)=0; 
    allocate(tmp%iobs_po_id(iobs_num_em,3)) !1:x,        2:y,        3:z
    allocate(tmp%iobs_po_pe(iobs_num_em))
    allocate(tmp%iobs_li_pe(iobs_num_em,3)) !1:x-line,   2:y-line,   3:z-line
    allocate(tmp%iobs_pl_pe(iobs_num_em,3)) !1:xy-plane, 2:yz-plane, 3:xz-plane
    tmp%iobs_po_id(:,:)=0; tmp%iobs_po_pe(:)=0; tmp%iobs_li_pe(:,:)=0; tmp%iobs_pl_pe(:,:)=0; 
    tmp%e_max=0.0d0; tmp%h_max=0.0d0;
    
    !search observation point
    do ii=1,iobs_num_em
      call eh_find_point(obs_loc_em(ii,:),&
                         tmp%iobs_po_id(ii,:),tmp%iobs_po_pe(ii),tmp%iobs_li_pe(ii,:),tmp%iobs_pl_pe(ii,:),&
                         grid%ng_sta,grid%ng_end,&
                         minval(grid%lg_sta)-tmp%Nd,maxval(grid%lg_end)+tmp%Nd,tmp%coo)
    end do
    
    !write information
    if(comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "**************************"
      if(iobs_num_em==1) then
        write(*,*) "Observation point is placed at"
      else
        write(*,*) "Observation points are placed at"
      end if
      do ii=1,iobs_num_em
        write(*,'(I3,A,3ES14.5)') ii,":",(tmp%coo(tmp%iobs_po_id(ii,ix),ix)*ulength_from_au,ix=1,3)
      end do
      write(*,*) "**************************"
      do ii=1,iobs_num_em
        write(save_name,*) ii
        save_name=trim(adjustl(directory))//'obs'//trim(adjustl(save_name))//'_at_point.data'
        open(tmp%ifn,file=save_name)
        select case(unit_system)
        case('au','a.u.')
          write(tmp%ifn,'(A)') "# time[a.u.], Ex[a.u.], Ey[a.u.], Ez[a.u.], Hx[a.u.], Hy[a.u.], Hz[a.u.]" 
        case('A_eV_fs')
          write(tmp%ifn,'(A)') "# time[fs], Ex[V/Ang.], Ey[V/Ang.], Ez[V/Ang.], Hx[A/Ang.], Hy[A/Ang.], Hz[A/Ang.]" 
        end select
        close(tmp%ifn)
      end do
    end if
  end if
  
  !check incident current source condition
  select case(wave_input)
  case('source')
    !source1
    if    (ek_dir1(1)==0.0d0.and.ek_dir1(2)==0.0d0.and.ek_dir1(3)==0.0d0) then 
      tmp%inc_dist1='none'
    elseif(ek_dir1(1)==1.0d0.and.ek_dir1(2)==0.0d0.and.ek_dir1(3)==0.0d0) then
      if(epdir_re1(1)/=0.0d0.or.epdir_im1(1)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re1(1) and epdir_im1(1):"
          write(*,*) "For theory = Maxwell and ek_dir1(1) = 1.0d0, epdir_re1(1) and epdir_im1(1) must be 0.0d0."
        end if
        stop
      else
        tmp%inc_dist1='yz-plane'
      end if
    elseif(ek_dir1(1)==0.0d0.and.ek_dir1(2)==1.0d0.and.ek_dir1(3)==0.0d0) then
      if(epdir_re1(2)/=0.0d0.or.epdir_im1(2)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re1(2) and epdir_im1(2):"
          write(*,*) "For theory = Maxwell and ek_dir1(2) = 1.0d0, epdir_re1(2) and epdir_im1(2) must be 0.0d0."
        end if
        stop
      else
        tmp%inc_dist1='xz-plane'
      end if
    elseif(ek_dir1(1)==0.0d0.and.ek_dir1(2)==0.0d0.and.ek_dir1(3)==1.0d0) then
      if(epdir_re1(3)/=0.0d0.or.epdir_im1(3)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re1(3) and epdir_im1(3):"
          write(*,*) "For theory = Maxwell and ek_dir1(3) = 1.0d0, epdir_re1(3) and epdir_im1(3) must be 0.0d0."
        end if
        stop
      else
        tmp%inc_dist1='xy-plane'
      end if
    else
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "invalid ek_dir1:"
        write(*,*) "For theory = Maxwell, ek_dir1 is only allowed by"
        write(*,*) "(0d0,0d0,0d0),(1d0,0d0,0d0),(0d0,1d0,0d0),or (0d0,0d0,1d0)."
      end if
      stop
    end if
    
    !source2
    if    (ek_dir2(1)==0.0d0.and.ek_dir2(2)==0.0d0.and.ek_dir2(3)==0.0d0) then 
      tmp%inc_dist2='none'
    elseif(ek_dir2(1)==1.0d0.and.ek_dir2(2)==0.0d0.and.ek_dir2(3)==0.0d0) then
      if(epdir_re2(1)/=0.0d0.or.epdir_im2(1)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re2(1) and epdir_im2(1):"
          write(*,*) "For theory = Maxwell and ek_dir2(1) = 1.0d0, epdir_re2(1) and epdir_im2(1) must be 0.0d0."
        end if
        stop
      else
        tmp%inc_dist2='yz-plane'
      end if
    elseif(ek_dir2(1)==0.0d0.and.ek_dir2(2)==1.0d0.and.ek_dir2(3)==0.0d0) then
      if(epdir_re2(2)/=0.0d0.or.epdir_im2(2)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re2(2) and epdir_im2(2):"
          write(*,*) "For theory = Maxwell and ek_dir2(2) = 1.0d0, epdir_re2(2) and epdir_im2(2) must be 0.0d0."
        end if
        stop
      else
        tmp%inc_dist2='xz-plane'
      end if
    elseif(ek_dir2(1)==0.0d0.and.ek_dir2(2)==0.0d0.and.ek_dir2(3)==1.0d0) then
      if(epdir_re2(3)/=0.0d0.or.epdir_im2(3)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re2(3) and epdir_im2(3):"
          write(*,*) "For theory = Maxwell and ek_dir2(3) = 1.0d0, epdir_re2(3) and epdir_im2(3) must be 0.0d0."
        end if
        stop
      else
        tmp%inc_dist2='xy-plane'
      end if
    else
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "invalid ek_dir2:"
        write(*,*) "For theory = Maxwell, ek_dir2 is only allowed by"
        write(*,*) "(0d0,0d0,0d0),(1d0,0d0,0d0),(0d0,1d0,0d0),or (0d0,0d0,1d0)."
      end if
      stop
    end if
  case('point','x-line','y-line','z-line')
    !these selection are for debug
    tmp%inc_dist1=wave_input;
    if(comm_is_root(nproc_id_global)) write(*,*) trim(wave_input), "source is used."
  case default
    stop 'invalid wave_input. For theory = Maxwell, wave_input is only allowed by source.'
  end select
  
  !prepare incident current source
  if((tmp%inc_dist1=='none').and.(tmp%inc_dist2=='none')) then
    tmp%inc_num=0
  else
    tmp%inc_num=2
  end if
  if(tmp%inc_num>0) then
    !set initial
    allocate(tmp%inc_po_id(iobs_num_em,3)) !1:x,        2:y,        3:z
    allocate(tmp%inc_po_pe(iobs_num_em))
    allocate(tmp%inc_li_pe(iobs_num_em,3)) !1:x-line,   2:y-line,   3:z-line
    allocate(tmp%inc_pl_pe(iobs_num_em,3)) !1:xy-plane, 2:yz-plane, 3:xz-plane
    tmp%inc_po_id(:,:)=0; tmp%inc_po_pe(:)=0; tmp%inc_li_pe(:,:)=0; tmp%inc_pl_pe(:,:)=0; 
    do ii=1,3
      tmp%c2_inc_xyz(ii)=(tmp%c_0/tmp%rep(0)*grid%dt) &
                         /(1.0d0+2.0d0*pi*tmp%sig(0)/tmp%rep(0)*grid%dt) &
                         *2.0d0/( grid%hgs(ii)*sqrt(tmp%rmu(0)/tmp%rep(0)) )
    end do
    
    !search incident current source point and check others
    if(tmp%inc_dist1/='none') then
      ii=1
      call eh_find_point(source_loc1(:),&
                         tmp%inc_po_id(ii,:),tmp%inc_po_pe(ii),tmp%inc_li_pe(ii,:),tmp%inc_pl_pe(ii,:),&
                         grid%ng_sta,grid%ng_end,&
                         minval(grid%lg_sta(:))-tmp%Nd,maxval(grid%lg_end(:))+tmp%Nd,tmp%coo(:,:))
      select case(ae_shape1)
      case("Ecos2","Acos2")
        continue
      case default
        stop 'set ae_shape1 to "Ecos2", or "Acos2"'
      end select
      diff_cep=(phi_cep1-0.25d0)*2.d0-int((phi_cep1-0.25d0)*2.d0)
      if(ae_shape1=="Ecos2".and.abs(diff_cep)>=1.d-12)then
        stop "phi_cep1 must be equal to 0.25+0.5*i when Ecos2 is specified for ae_shape1."
      end if
      if(rlaser_int_wcm2_1/=-1d0) &
        amplitude1=sqrt(rlaser_int_wcm2_1)*1.0d2*2.74492d1/(5.14223d11)!I[W/cm^2]->E[a.u.]
    end if
    if(tmp%inc_dist2/='none') then
      ii=2
      call eh_find_point(source_loc2(:),&
                         tmp%inc_po_id(ii,:),tmp%inc_po_pe(ii),tmp%inc_li_pe(ii,:),tmp%inc_pl_pe(ii,:),&
                         grid%ng_sta,grid%ng_end,&
                         minval(grid%lg_sta(:))-tmp%Nd,maxval(grid%lg_end(:))+tmp%Nd,tmp%coo(:,:))
      select case(ae_shape2)
      case("Ecos2","Acos2")
        continue
      case default
        stop 'set ae_shape2 to "Ecos2", or "Acos2"'
      end select
      diff_cep=(phi_cep2-0.25d0)*2.d0-int((phi_cep2-0.25d0)*2.d0)
      if(ae_shape2=="Ecos2".and.abs(diff_cep)>=1.d-12)then
        stop "phi_cep2 must be equal to 0.25+0.5*i when Ecos2 is specified for ae_shape2."
      end if
      if(rlaser_int_wcm2_2/=-1d0) &
        amplitude2=sqrt(rlaser_int_wcm2_2)*1.0d2*2.74492d1/(5.14223d11)!I[W/cm^2]->E[a.u.]
    end if
    
    !write information
    if(comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "**************************"
      if((tmp%inc_dist1=='none').or.(tmp%inc_dist2=='none')) then
        write(*,*) "Incident current source is placed at"
      else
        write(*,*) "Incident current sources are placed at"
      end if
      if(tmp%inc_dist1/='none') then
        ii=1
        write(*,'(I8,A,3ES14.5,A)') ii,":",(tmp%coo(tmp%inc_po_id(ii,ix),ix)*ulength_from_au,ix=1,3)
        write(*,'(A,3ES14.5)') " ek_dir1:",ek_dir1
      end if
      if(tmp%inc_dist2/='none') then
        ii=2
        write(*,'(I8,A,3ES14.5,A)') ii,":",(tmp%coo(tmp%inc_po_id(ii,ix),ix)*ulength_from_au,ix=1,3)
        write(*,'(A,3ES14.5)') " ek_dir2:",ek_dir2
      end if
      write(*,*) "**************************"
    end if
  end if
  
  !write strat
  if(comm_is_root(nproc_id_global)) then
    write(*,*)
    write(*,*) "**************************"
    write(*,*) "FDTD start"
    write(*,*) "**************************"
    write(*,*) "timestep"
    write(*,*) "-------------------------------------------------------"
  end if
  
contains
  
  !=========================================================================================
  != e and h allocation ====================================================================
  subroutine eh_allocate
    implicit none
    
    !e
    allocate(tmp%ex_y(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c1_ex_y(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c2_ex_y(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    tmp%ex_y(:,:,:)=0.0d0; tmp%c1_ex_y(:,:,:)=0.0d0; tmp%c2_ex_y(:,:,:)=0.0d0;
    allocate(tmp%ex_z(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c1_ex_z(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c2_ex_z(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    tmp%ex_z(:,:,:)=0.0d0; tmp%c1_ex_z(:,:,:)=0.0d0; tmp%c2_ex_z(:,:,:)=0.0d0;
    allocate(tmp%ey_z(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c1_ey_z(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c2_ey_z(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    tmp%ey_z(:,:,:)=0.0d0; tmp%c1_ey_z(:,:,:)=0.0d0; tmp%c2_ey_z(:,:,:)=0.0d0;
    allocate(tmp%ey_x(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c1_ey_x(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c2_ey_x(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    tmp%ey_x(:,:,:)=0.0d0; tmp%c1_ey_x(:,:,:)=0.0d0; tmp%c2_ey_x(:,:,:)=0.0d0;
    allocate(tmp%ez_x(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c1_ez_x(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c2_ez_x(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    tmp%ez_x(:,:,:)=0.0d0; tmp%c1_ez_x(:,:,:)=0.0d0; tmp%c2_ez_x(:,:,:)=0.0d0;
    allocate(tmp%ez_y(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c1_ez_y(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c2_ez_y(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    tmp%ez_y(:,:,:)=0.0d0; tmp%c1_ez_y(:,:,:)=0.0d0; tmp%c2_ez_y(:,:,:)=0.0d0;
    
    !h
    allocate(tmp%hx_y(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c1_hx_y(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c2_hx_y(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    tmp%hx_y(:,:,:)=0.0d0; tmp%c1_hx_y(:,:,:)=0.0d0; tmp%c2_hx_y(:,:,:)=0.0d0;
    allocate(tmp%hx_z(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c1_hx_z(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c2_hx_z(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    tmp%hx_z(:,:,:)=0.0d0; tmp%c1_hx_z(:,:,:)=0.0d0; tmp%c2_hx_z(:,:,:)=0.0d0;
    allocate(tmp%hy_z(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c1_hy_z(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c2_hy_z(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    tmp%hy_z(:,:,:)=0.0d0; tmp%c1_hy_z(:,:,:)=0.0d0; tmp%c2_hy_z(:,:,:)=0.0d0;
    allocate(tmp%hy_x(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c1_hy_x(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c2_hy_x(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    tmp%hy_x(:,:,:)=0.0d0; tmp%c1_hy_x(:,:,:)=0.0d0; tmp%c2_hy_x(:,:,:)=0.0d0;
    allocate(tmp%hz_x(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c1_hz_x(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c2_hz_x(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    tmp%hz_x(:,:,:)=0.0d0; tmp%c1_hz_x(:,:,:)=0.0d0; tmp%c2_hz_x(:,:,:)=0.0d0;
    allocate(tmp%hz_y(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                      grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                      grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c1_hz_y(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    allocate(tmp%c2_hz_y(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                         grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                         grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd))
    tmp%hz_y(:,:,:)=0.0d0; tmp%c1_hz_y(:,:,:)=0.0d0; tmp%c2_hz_y(:,:,:)=0.0d0;
    
  end subroutine eh_allocate
  
  !=========================================================================================
  != set fdtd coefficient ==================================================================
  subroutine eh_coeff
    implicit none
    real(8)  :: c1_e,c2_e_x,c2_e_y,c2_e_z,c1_h,c2_h_x,c2_h_y,c2_h_z,c2_j
    
    !set constatn parameter
    tmp%rep(ii)=epsilon(ii); tmp%rmu(ii)=rmu(ii); tmp%sig(ii)=sigma(ii);
    
    !prepare coefficient
    c1_e  =(1.0d0-2.0d0*pi*tmp%sig(ii)/tmp%rep(ii)*grid%dt) &
           /(1.0d0+2.0d0*pi*tmp%sig(ii)/tmp%rep(ii)*grid%dt)
    c2_e_x=(tmp%c_0/tmp%rep(ii)*grid%dt) &
           /(1.0d0+2.0d0*pi*tmp%sig(ii)/tmp%rep(ii)*grid%dt)/grid%hgs(1)
    c2_e_y=(tmp%c_0/tmp%rep(ii)*grid%dt) &
           /(1.0d0+2.0d0*pi*tmp%sig(ii)/tmp%rep(ii)*grid%dt)/grid%hgs(2)
    c2_e_z=(tmp%c_0/tmp%rep(ii)*grid%dt) &
           /(1.0d0+2.0d0*pi*tmp%sig(ii)/tmp%rep(ii)*grid%dt)/grid%hgs(3)
    call comm_bcast(c1_e,  nproc_group_global)
    call comm_bcast(c2_e_x,nproc_group_global)
    call comm_bcast(c2_e_y,nproc_group_global)
    call comm_bcast(c2_e_z,nproc_group_global)
    c1_h=1.0d0
    c2_h_x=tmp%c_0/tmp%rmu(ii)*grid%dt/grid%hgs(1)
    c2_h_y=tmp%c_0/tmp%rmu(ii)*grid%dt/grid%hgs(2)
    c2_h_z=tmp%c_0/tmp%rmu(ii)*grid%dt/grid%hgs(3)
    call comm_bcast(c1_h,  nproc_group_global)
    call comm_bcast(c2_h_x,nproc_group_global)
    call comm_bcast(c2_h_y,nproc_group_global)
    call comm_bcast(c2_h_z,nproc_group_global)
    c2_j=(4.0d0*pi/tmp%rep(ii)*grid%dt) &
         /(1.0d0+2.0d0*pi*tmp%sig(ii)/tmp%rep(ii)*grid%dt)
    call comm_bcast(c2_j,nproc_group_global)
    
    !check type_media
    select case(type_media(ii))
    case('PEC','Pec','pec')
      c1_e=0.0d0; c2_e_x=0.0d0; c2_e_y=0.0d0; c2_e_z=0.0d0;
    case('DRUDE','Drude','drude','D','d')
      do iz=grid%ng_sta(3),grid%ng_end(3)
      do iy=grid%ng_sta(2),grid%ng_end(2)
      do ix=grid%ng_sta(1),grid%ng_end(1)
        if(grid%imedia(ix,iy,iz)==ii) then
          tmp%idx_d(ix-1:ix,iy,iz,icount_d)=1;
          tmp%idy_d(ix,iy-1:iy,iz,icount_d)=1;
          tmp%idz_d(ix,iy,iz-1:iz,icount_d)=1;
        end if
      end do
      end do
      end do
      tmp%c1_j_d(icount_d)=(1.0d0-gamma_d(ii)*grid%dt/2.0d0)           / (1.0d0+gamma_d(ii)*grid%dt/2.0d0);
      tmp%c2_j_d(icount_d)=((omega_p_d(ii)**2.0d0)*grid%dt/(4.0d0*pi)) / (1.0d0+gamma_d(ii)*grid%dt/2.0d0);
      icount_d=icount_d+1
    end select
    
    !set coefficient
    if(ii==0) then
      tmp%c1_ex_y(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c1_e
      tmp%c2_ex_y(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c2_e_y
      tmp%c1_ex_z(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c1_e
      tmp%c2_ex_z(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=-c2_e_z
          
      tmp%c1_ey_z(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c1_e
      tmp%c2_ey_z(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c2_e_z
      tmp%c1_ey_x(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c1_e
      tmp%c2_ey_x(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=-c2_e_x
        
      tmp%c1_ez_x(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c1_e
      tmp%c2_ez_x(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c2_e_x
      tmp%c1_ez_y(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c1_e
      tmp%c2_ez_y(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=-c2_e_y
        
      tmp%c1_hx_y(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c1_h
      tmp%c2_hx_y(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=-c2_h_y
      tmp%c1_hx_z(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c1_h
      tmp%c2_hx_z(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c2_h_z
        
      tmp%c1_hy_z(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c1_h
      tmp%c2_hy_z(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=-c2_h_z
      tmp%c1_hy_x(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c1_h
      tmp%c2_hy_x(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c2_h_x
      
      tmp%c1_hz_x(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c1_h
      tmp%c2_hz_x(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=-c2_h_x
      tmp%c1_hz_y(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c1_h
      tmp%c2_hz_y(grid%ng_sta(1):grid%ng_end(1),&
                  grid%ng_sta(2):grid%ng_end(2),&
                  grid%ng_sta(3):grid%ng_end(3))=c2_h_y
                  
      tmp%c2_jx(grid%ng_sta(1):grid%ng_end(1),&
                grid%ng_sta(2):grid%ng_end(2),&
                grid%ng_sta(3):grid%ng_end(3))=-c2_j
      tmp%c2_jy(grid%ng_sta(1):grid%ng_end(1),&
                grid%ng_sta(2):grid%ng_end(2),&
                grid%ng_sta(3):grid%ng_end(3))=-c2_j
      tmp%c2_jz(grid%ng_sta(1):grid%ng_end(1),&
                grid%ng_sta(2):grid%ng_end(2),&
                grid%ng_sta(3):grid%ng_end(3))=-c2_j
    else
      do iz=grid%ng_sta(3),grid%ng_end(3)
      do iy=grid%ng_sta(2),grid%ng_end(2)
      do ix=grid%ng_sta(1),grid%ng_end(1)
        if(grid%imedia(ix,iy,iz)==ii) then
          !ex
          tmp%c1_ex_y(ix-1:ix,iy,iz)=c1_e; tmp%c2_ex_y(ix-1:ix,iy,iz)= c2_e_y;
          tmp%c1_ex_z(ix-1:ix,iy,iz)=c1_e; tmp%c2_ex_z(ix-1:ix,iy,iz)=-c2_e_z;
          
          !ey
          tmp%c1_ey_z(ix,iy-1:iy,iz)=c1_e; tmp%c2_ey_z(ix,iy-1:iy,iz)= c2_e_z;
          tmp%c1_ey_x(ix,iy-1:iy,iz)=c1_e; tmp%c2_ey_x(ix,iy-1:iy,iz)=-c2_e_x;
          
          !ez
          tmp%c1_ez_x(ix,iy,iz-1:iz)=c1_e; tmp%c2_ez_x(ix,iy,iz-1:iz)= c2_e_x;
          tmp%c1_ez_y(ix,iy,iz-1:iz)=c1_e; tmp%c2_ez_y(ix,iy,iz-1:iz)=-c2_e_y;
          
          !hx
          tmp%c1_hx_y(ix,iy-1:iy,iz-1:iz)=c1_h; tmp%c2_hx_y(ix,iy-1:iy,iz-1:iz)=-c2_h_y;
          tmp%c1_hx_z(ix,iy-1:iy,iz-1:iz)=c1_h; tmp%c2_hx_z(ix,iy-1:iy,iz-1:iz)= c2_h_z;
          
          !hy
          tmp%c1_hy_z(ix-1:ix,iy,iz-1:iz)=c1_h; tmp%c2_hy_z(ix-1:ix,iy,iz-1:iz)=-c2_h_z;
          tmp%c1_hy_x(ix-1:ix,iy,iz-1:iz)=c1_h; tmp%c2_hy_x(ix-1:ix,iy,iz-1:iz)= c2_h_x;
          
          !hz
          tmp%c1_hz_x(ix-1:ix,iy-1:iy,iz)=c1_h; tmp%c2_hz_x(ix-1:ix,iy-1:iy,iz)=-c2_h_x;
          tmp%c1_hz_y(ix-1:ix,iy-1:iy,iz)=c1_h; tmp%c2_hz_y(ix-1:ix,iy-1:iy,iz)= c2_h_y;
          
          !j
          tmp%c2_jx(ix-1:ix,iy,iz)=-c2_j; tmp%c2_jy(ix,iy-1:iy,iz)=-c2_j; tmp%c2_jz(ix,iy,iz-1:iz)=-c2_j;
        end if
      end do
      end do
      end do
    end if
    
  end subroutine eh_coeff
  
  !=========================================================================================
  != set pml ===============================================================================
  subroutine eh_set_pml(idir,c1_e1,c2_e1,c1_e2,c2_e2,c1_h1,c2_h1,c1_h2,c2_h2)
    implicit none
    integer,intent(in)  :: idir
    real(8),intent(out) :: c1_e1(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                                 grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                                 grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
                           c2_e1(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                                 grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                                 grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
                           c1_e2(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                                 grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                                 grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
                           c2_e2(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                                 grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                                 grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
                           c1_h1(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                                 grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                                 grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
                           c2_h1(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                                 grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                                 grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
                           c1_h2(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                                 grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                                 grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd),&
                           c2_h2(grid%ng_sta(1)-tmp%Nd:grid%ng_end(1)+tmp%Nd,&
                                 grid%ng_sta(2)-tmp%Nd:grid%ng_end(2)+tmp%Nd,&
                                 grid%ng_sta(3)-tmp%Nd:grid%ng_end(3)+tmp%Nd)
    integer :: ista,iend
    real(8) :: pml_del,s_max
    real(8) :: s_l(tmp%ipml_l+1),sh_l(tmp%ipml_l), &
               c1_pml(tmp%ipml_l+1),c2_pml(tmp%ipml_l+1),c1_pml_h(tmp%ipml_l),c2_pml_h(tmp%ipml_l)

    !set pml conductivity
    pml_del=grid%hgs(idir)
    s_max=-(tmp%pml_m+1.0d0)*log(tmp%pml_r)/(2.0d0*dble(tmp%ipml_l)*pml_del) &
          *tmp%c_0/(4.0d0*pi)*sqrt(tmp%rep(0)/tmp%rmu(0));
    do ii=1,(tmp%ipml_l+1)
      s_l(ii)=s_max*(&
                    (dble(tmp%ipml_l)*pml_del-(dble(ii)-1.0d0)*pml_del)/(dble(tmp%ipml_l)*pml_del)&
                    )**tmp%pml_m;
    end do
    do ii=1,tmp%ipml_l
      sh_l(ii)=(tmp%rmu(0)/tmp%rep(0)) &
               *s_max*(&
                      (dble(tmp%ipml_l)*pml_del-(dble(ii)-0.5d0)*pml_del)/(dble(tmp%ipml_l)*pml_del)&
                      )**tmp%pml_m;
    end do
    
    !set pml coefficient
    do ii=1,(tmp%ipml_l+1)
      c1_pml(ii)=(1.0d0-2.0d0*pi*s_l(ii)/tmp%rep(0)*grid%dt) &
                 /(1.0d0+2.0d0*pi*s_l(ii)/tmp%rep(0)*grid%dt);
      c2_pml(ii)=(tmp%c_0/tmp%rep(0)*grid%dt) &
                 /(1.0d0+2.0d0*pi*s_l(ii)/tmp%rep(0)*grid%dt)/pml_del
    end do
    call comm_bcast(c1_pml,nproc_group_global)
    call comm_bcast(c2_pml,nproc_group_global)
    do ii=1,tmp%ipml_l
      c1_pml_h(ii)=(1.0d0-2.0d0*pi*sh_l(ii)/tmp%rmu(0)*grid%dt) &
                   /(1.0d0+2.0d0*pi*sh_l(ii)/tmp%rmu(0)*grid%dt);
      c2_pml_h(ii)=(tmp%c_0/tmp%rmu(0)*grid%dt) &
                   /(1.0d0+2.0d0*pi*sh_l(ii)/tmp%rmu(0)*grid%dt)/pml_del
    end do
    call comm_bcast(c1_pml_h,nproc_group_global)
    call comm_bcast(c2_pml_h,nproc_group_global)
    
    !set pml(bottom)
    if((grid%i_bc(idir,1)==1).and.(grid%ng_sta(idir)<=(grid%lg_sta(idir)+tmp%ipml_l))) then
      !e
      iend=grid%lg_sta(idir)+tmp%ipml_l
      if(grid%ng_end(idir)<iend) then
        iend=grid%ng_end(idir)
      end if
      icount=1
      do ii=grid%ng_sta(idir),iend
        if(idir==1) then
          c1_e1(ii,:,:)= c1_pml(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c2_e1(ii,:,:)=-c2_pml(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c1_e2(ii,:,:)= c1_pml(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c2_e2(ii,:,:)= c2_pml(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
        elseif(idir==2) then
          c1_e1(:,ii,:)= c1_pml(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c2_e1(:,ii,:)=-c2_pml(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c1_e2(:,ii,:)= c1_pml(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c2_e2(:,ii,:)= c2_pml(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
        elseif(idir==3) then
          c1_e1(:,:,ii)= c1_pml(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c2_e1(:,:,ii)=-c2_pml(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c1_e2(:,:,ii)= c1_pml(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c2_e2(:,:,ii)= c2_pml(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
        end if
        icount=icount+1
      end do
      
      !h
      if(iend==(grid%lg_sta(idir)+tmp%ipml_l)) then
        iend=iend-1
      end if
      icount=1
      do ii=grid%ng_sta(idir),iend
        if(idir==1) then
          c1_h1(ii,:,:)= c1_pml_h(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c2_h1(ii,:,:)= c2_pml_h(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c1_h2(ii,:,:)= c1_pml_h(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c2_h2(ii,:,:)=-c2_pml_h(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
        elseif(idir==2) then
          c1_h1(:,ii,:)= c1_pml_h(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c2_h1(:,ii,:)= c2_pml_h(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c1_h2(:,ii,:)= c1_pml_h(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c2_h2(:,ii,:)=-c2_pml_h(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
        elseif(idir==3) then
          c1_h1(:,:,ii)= c1_pml_h(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c2_h1(:,:,ii)= c2_pml_h(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c1_h2(:,:,ii)= c1_pml_h(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
          c2_h2(:,:,ii)=-c2_pml_h(grid%ng_sta(idir)-grid%lg_sta(idir)+icount)
        end if
        icount=icount+1
      end do
    end if
    
    !set pml(top)
    if((grid%i_bc(idir,2)==1).and.(grid%ng_end(idir)>=(grid%lg_end(idir)-tmp%ipml_l))) then
      !e
      ista=grid%lg_end(idir)-tmp%ipml_l
      if(grid%ng_sta(idir)>ista) then
        ista=grid%ng_sta(idir)
      end if
      icount=1
      do ii=ista,grid%ng_end(idir)
        if(idir==1) then
          c1_e1(ii,:,:)= c1_pml((tmp%ipml_l+1)-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c2_e1(ii,:,:)=-c2_pml((tmp%ipml_l+1)-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c1_e2(ii,:,:)= c1_pml((tmp%ipml_l+1)-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c2_e2(ii,:,:)= c2_pml((tmp%ipml_l+1)-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
        elseif(idir==2) then
          c1_e1(:,ii,:)= c1_pml((tmp%ipml_l+1)-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c2_e1(:,ii,:)=-c2_pml((tmp%ipml_l+1)-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c1_e2(:,ii,:)= c1_pml((tmp%ipml_l+1)-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c2_e2(:,ii,:)= c2_pml((tmp%ipml_l+1)-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
        elseif(idir==3) then
          c1_e1(:,:,ii)= c1_pml((tmp%ipml_l+1)-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c2_e1(:,:,ii)=-c2_pml((tmp%ipml_l+1)-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c1_e2(:,:,ii)= c1_pml((tmp%ipml_l+1)-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c2_e2(:,:,ii)= c2_pml((tmp%ipml_l+1)-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
        end if
        icount=icount+1
      end do
      
      !h
      if(grid%ng_end(idir)==grid%lg_end(idir)) then
        iend=grid%ng_end(idir)-1
      else
        iend=grid%ng_end(idir)
      end if
      icount=1
      do ii=ista,iend
        if(idir==1) then
          c1_h1(ii,:,:)= c1_pml_h(tmp%ipml_l-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c2_h1(ii,:,:)= c2_pml_h(tmp%ipml_l-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c1_h2(ii,:,:)= c1_pml_h(tmp%ipml_l-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c2_h2(ii,:,:)=-c2_pml_h(tmp%ipml_l-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
        elseif(idir==2) then
          c1_h1(:,ii,:)= c1_pml_h(tmp%ipml_l-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c2_h1(:,ii,:)= c2_pml_h(tmp%ipml_l-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c1_h2(:,ii,:)= c1_pml_h(tmp%ipml_l-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c2_h2(:,ii,:)=-c2_pml_h(tmp%ipml_l-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
        elseif(idir==3) then
          c1_h1(:,:,ii)= c1_pml_h(tmp%ipml_l-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c2_h1(:,:,ii)= c2_pml_h(tmp%ipml_l-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c1_h2(:,:,ii)= c1_pml_h(tmp%ipml_l-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
          c2_h2(:,:,ii)=-c2_pml_h(tmp%ipml_l-(ista-(grid%lg_end(idir)-tmp%ipml_l)+(icount-1)))
        end if
        icount=icount+1
      end do
    end if
    
  end subroutine eh_set_pml
  
end subroutine eh_init

!=========================================================================================
!= input fdtd shape data =================================================================
subroutine eh_input_shape(ifn,ng_sta,ng_end,lg_sta,lg_end,imat,format)
  use salmon_parallel,      only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use inputoutput,          only: shape_file
  implicit none
  integer,intent(in)      :: ifn
  integer,intent(in)      :: ng_sta(3),ng_end(3),lg_sta(3),lg_end(3)
  integer,intent(out)     :: imat(ng_sta(1):ng_end(1),&
                                  ng_sta(2):ng_end(2),&
                                  ng_sta(3):ng_end(3))
  character(2),intent(in) :: format
  integer,allocatable     :: itmp1d(:)
  integer                 :: inum(3),inum_check(3)
  integer                 :: ii,ix,iy,iz,iflag_x,iflag_y
  
  !open file
  open(ifn,file=trim(shape_file), status='old')
  
  if(trim(format)=='cu') then
    !check grid information
    inum(:)=lg_end(:)-lg_sta(:)+1
    read(ifn,*); read (ifn,*); read (ifn,*); !skip
    allocate(itmp1d(4))
    read (ifn,*) itmp1d; inum_check(1)=itmp1d(1);
    read (ifn,*) itmp1d; inum_check(2)=itmp1d(1);
    read (ifn,*) itmp1d; inum_check(3)=itmp1d(1);
    deallocate(itmp1d)
    do ii=1,3
      if(inum(ii)/=inum_check(ii)) then
        if(comm_is_root(nproc_id_global)) write(*,*) "al_em or dl_em does not mutch shape file."
        stop
      end if
    end do
    read (ifn,*); !skip
    
    !input shape
    allocate(itmp1d(lg_sta(3):lg_end(3)))
    do ix=lg_sta(1),lg_end(1)
    do iy=lg_sta(2),lg_end(2)
      read (ifn,*) itmp1d
      iflag_x=0; iflag_y=0;
      if(ix>=ng_sta(1) .and. ix<=ng_end(1)) iflag_x=1
      if(iy>=ng_sta(2) .and. iy<=ng_end(2)) iflag_y=1
      if(iflag_x==1 .and. iflag_y==1) then
        do iz=ng_sta(3),ng_end(3)
          imat(ix,iy,iz)=itmp1d(iz)
        end do
      end if
    end do
    end do
  elseif(trim(format)=='mp') then
  end if
  
  !close file
  close(ifn)

end subroutine eh_input_shape

!=========================================================================================
!= find point and set corresponding processor element ====================================
subroutine eh_find_point(rloc,id,ipo,ili,ipl,ista,iend,icoo_sta,icoo_end,coo)
  use salmon_parallel,      only: nproc_id_global,nproc_size_global,nproc_group_global
  use salmon_communication, only: comm_summation
  implicit none
  real(8),intent(in)    :: rloc(3)
  integer,intent(inout) :: id(3)
  integer,intent(out)   :: ipo
  integer,intent(out)   :: ili(3),ipl(3)
  integer,intent(in)    :: ista(3),iend(3)
  integer,intent(in)    :: icoo_sta,icoo_end
  real(8),intent(in)    :: coo(icoo_sta:icoo_end,3)
  integer               :: ii,ix,iy,iz,ipe_tmp,i1,i1s,i2,i2s
  integer               :: id_tmp(3),id_tmp2(3)
  real(8)               :: err(0:nproc_size_global-1),err2(0:nproc_size_global-1)
  real(8)               :: err_tmp
  
  !set initial value
  err(:)              =0.0d0
  err(nproc_id_global)=1.0d9
  err2(:)             =0.0d0
  id_tmp(:)           =0
  id_tmp2(:)          =0
  
  !find observation point in each PE
  do iz=ista(3),iend(3)
  do iy=ista(2),iend(2)
  do ix=ista(1),iend(1)
    err_tmp=sqrt( (coo(ix,1)-rloc(1))**2.0d0 &
                 +(coo(iy,2)-rloc(2))**2.0d0 &
                 +(coo(iz,3)-rloc(3))**2.0d0 )
    if(err(nproc_id_global)>err_tmp) then
      err(nproc_id_global)=err_tmp
      id_tmp(1)=ix; id_tmp(2)=iy; id_tmp(3)=iz;
    end if
  end do
  end do
  end do
  call comm_summation(err,err2,nproc_size_global,nproc_group_global)
  
  !determine and share id + determine pe including the point
  ipe_tmp=-1; err_tmp=1.0d9;
  do ii=0,nproc_size_global-1
    if(err_tmp>err2(ii)) then
      err_tmp=err2(ii)
      ipe_tmp=ii
    end if
  end do
  if(nproc_id_global==ipe_tmp) then
    ipo=1;
  else
    ipo=0; id_tmp(:)=0;
  end if
  call comm_summation(id_tmp,id_tmp2,3,nproc_group_global)
  id(:)=id_tmp2(:)
  
  !determine pe including the line
  do ii=1,3
    if(ii==1) then     !x-line(searching at yz-plane)
      i1s=3; i2s=2;
    elseif(ii==2) then !x-line(searching at xz-plane)
      i1s=3; i2s=1;
    elseif(ii==3) then !z-line(searching at xy-plane)
      i1s=2; i2s=1;
    end if
    do i2=ista(i2s),iend(i2s)
    do i1=ista(i1s),iend(i1s)
      if( (i1==id(i1s)).and.(i2==id(i2s)) ) ili(ii)=1
    end do
    end do
  end do
  
  !determine pe including the plane
  do ii=1,3
    if(ii==1) then     !xy-plane(searching at z-line)
      i1s=3;
    elseif(ii==2) then !yz-plane(searching at x-line)
      i1s=1;
    elseif(ii==3) then !xz-plane(searching at y-line)
      i1s=2;
    end if
    do i1=ista(i1s),iend(i1s)
      if(i1==id(i1s)) ipl(ii)=1
    end do
  end do
  
end subroutine eh_find_point

!=========================================================================================
!= prepare GCEED =========================================================================
!= (This routine is temporary) ===========================================================
!= (With unifying ARTED and GCEED, this routine will be removed) =========================
subroutine eh_prep_GCEED(grid,tmp)
  use inputoutput,       only: nproc_domain,nproc_domain_s,num_kgrid,nproc_k,nproc_ob,isequential,iperiodic
  use salmon_parallel,   only: nproc_id_orbitalgrid,nproc_id_global,nproc_size_global
  use scf_data,          only: nproc_Mxin,nproc_Mxin_s,nproc_Mxin_mul,nproc_Mxin_mul_s_dm,nproc_Mxin_s_dm,&
                               k_sta,k_end,k_num,num_kpoints_3d,num_kpoints_rd,&
                               rLsize,Harray,Hgs,Hvol,imesh_oddeven,&
                               lg_sta,lg_end,lg_num, &
                               mg_sta,mg_end,mg_num, &
                               ng_sta,ng_end,ng_num,&
                               ista_Mx_ori,iend_Mx_ori,inum_Mx_ori,Nd, &
                               ista_Mxin,iend_Mxin,inum_Mxin,&
                               ista_Mxin_s,iend_Mxin_s,inum_Mxin_s, &
                               gridcoo,iwk_size,make_iwksta_iwkend
  use new_world_sub,     only: make_new_world
  use init_sendrecv_sub, only: init_updown,init_itype,init_sendrecv_matrix
  use persistent_comm,   only: init_persistent_requests
  use salmon_maxwell,    only: fdtd_grid,fdtd_tmp
  implicit none
  type(fdtd_grid)  :: grid
  type(fdtd_tmp)   :: tmp
  
  !set mpi condition
  num_kpoints_3d(1:3)=num_kgrid(1:3)
  num_kpoints_rd=num_kpoints_3d(1)*num_kpoints_3d(2)*num_kpoints_3d(3)
  nproc_Mxin=nproc_domain
  nproc_Mxin_s=nproc_domain_s
  call set_numcpu_scf
  nproc_Mxin_mul=nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)
  nproc_Mxin_mul_s_dm=nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)
  call make_new_world
  call setk(k_sta,k_end,k_num,num_kpoints_rd,nproc_k,nproc_id_orbitalgrid)
  
  !set grid
  rLsize(:,1)=grid%rlsize(:); Harray(:,1)=grid%hgs(:);
  Hgs(:)=Harray(:,1); Hvol=Hgs(1)*Hgs(2)*Hgs(3);
  call set_imesh_oddeven(1)
  call setlg(lg_sta,lg_end,lg_num,ista_Mx_ori,iend_Mx_ori,inum_Mx_ori, &
             Hgs,Nd,rLsize(:,1),imesh_oddeven,iperiodic)
  allocate(ista_Mxin(3,0:nproc_size_global-1),iend_Mxin(3,0:nproc_size_global-1), &
           inum_Mxin(3,0:nproc_size_global-1))
  call setmg(mg_sta,mg_end,mg_num,ista_Mxin,iend_Mxin,inum_Mxin, &
             lg_sta,lg_num,nproc_size_global,nproc_id_global,nproc_Mxin,nproc_k,nproc_ob,isequential)
  allocate(ista_Mxin_s(3,0:nproc_size_global-1),iend_Mxin_s(3,0:nproc_size_global-1))
  allocate(inum_Mxin_s(3,0:nproc_size_global-1))
  call setng(ng_sta,ng_end,ng_num,ista_Mxin_s,iend_Mxin_s,inum_Mxin_s, &
             nproc_size_global,nproc_id_global,nproc_Mxin,nproc_Mxin_s_dm,ista_Mxin,iend_Mxin,isequential)
  grid%lg_sta(:)=lg_sta(:); grid%lg_end(:)=lg_end(:);
  grid%ng_sta(:)=ng_sta(:); grid%ng_end(:)=ng_end(:);
  
  !set coordinate
  allocate(tmp%coo(minval(grid%lg_sta(:))-tmp%Nd:maxval(grid%lg_end(:))+tmp%Nd,3))
  call set_gridcoo
  tmp%coo(:,:)=gridcoo(:,:)
  
  !set sendrecv environment
  call init_updown
  call init_itype
  call init_sendrecv_matrix
  call init_persistent_requests
  iwk_size=12
  call make_iwksta_iwkend
  iwk_size=2
  
end subroutine eh_prep_GCEED
