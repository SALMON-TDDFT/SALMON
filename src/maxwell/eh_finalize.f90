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
subroutine eh_finalize(grid,tmp)
  use inputoutput,          only: utime_from_au,ulength_from_au,uenergy_from_au,unit_system,iperiodic,&
                                  ae_shape1,ae_shape2,e_impulse,sysname,nt_em,nenergy,de, &
                                  directory,iobs_num_em,iobs_samp_em
  use salmon_parallel,      only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use salmon_maxwell, only:fdtd_grid,fdtd_tmp
  implicit none
  type(fdtd_grid)     :: grid
  type(fdtd_tmp)      :: tmp
  integer             :: ii
  real(8),parameter   :: pi=3.141592653589793d0
  character(128)      :: save_name
  
  !output linear response(matter dipole pm and current jm are outputted: pm = -dip and jm = -curr)
  if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
    if(iperiodic==0) then
      !output time-dependent dipole data
      if(comm_is_root(nproc_id_global)) then
        save_name=trim(adjustl(directory))//'/'//trim(adjustl(sysname))//'_p.data'
        open(tmp%ifn,file=save_name)
        select case(unit_system)
        case('au','a.u.')
          write(tmp%ifn,'(A)') "# time[a.u.], dipoleMoment(x,y,z)[a.u.]" 
        case('A_eV_fs')
          write(tmp%ifn,'(A)') "# time[fs], dipoleMoment(x,y,z)[Ang.]" 
        end select
        do ii=1,nt_em
          write(tmp%ifn, '(E13.5)',advance="no")     tmp%time_lr(ii)*utime_from_au
          write(tmp%ifn, '(3E16.6e3)',advance="yes") -tmp%dip_lr(ii,:)*ulength_from_au
        end do
        close(tmp%ifn)
      end if
      
      !output lr data
      call eh_fourier(nt_em,nenergy,grid%dt,de,tmp%time_lr,tmp%dip_lr(:,1),tmp%fr_lr(:,1),tmp%fi_lr(:,1))
      call eh_fourier(nt_em,nenergy,grid%dt,de,tmp%time_lr,tmp%dip_lr(:,2),tmp%fr_lr(:,2),tmp%fi_lr(:,2))
      call eh_fourier(nt_em,nenergy,grid%dt,de,tmp%time_lr,tmp%dip_lr(:,3),tmp%fr_lr(:,3),tmp%fi_lr(:,3))
      if(comm_is_root(nproc_id_global)) then
        save_name=trim(adjustl(directory))//'/'//trim(adjustl(sysname))//'_lr.data'
        open(tmp%ifn,file=save_name)
        select case(unit_system)
        case('au','a.u.')
          write(tmp%ifn,'(A)') "# energy[a.u.], Re[alpha](x,y,z)[a.u.], Im[alpha](x,y,z)[a.u.], df/dE(x,y,z)[a.u.]"
        case('A_eV_fs')
          write(tmp%ifn,'(A)') "# energy[eV], Re[alpha](x,y,z)[Ang.**3], Im[alpha](x,y,z)[Ang.**3], df/dE(x,y,z)[1/eV]"
        end select
        do ii=0,nenergy
          write(tmp%ifn, '(E13.5)',advance="no")     dble(ii)*de*uenergy_from_au
          write(tmp%ifn, '(3E16.6e3)',advance="no")  tmp%fr_lr(ii,:)/(-e_impulse)*(ulength_from_au**3.0d0)
          write(tmp%ifn, '(3E16.6e3)',advance="no")  tmp%fi_lr(ii,:)/(-e_impulse)*(ulength_from_au**3.0d0)
          write(tmp%ifn, '(3E16.6e3)',advance="yes") 2.0d0*dble(ii)*de/pi*tmp%fi_lr(ii,:)/(-e_impulse)/uenergy_from_au
        end do
        close(tmp%ifn)
      end if
    elseif(iperiodic==3) then
      !output time-dependent dipole data
      if(comm_is_root(nproc_id_global)) then
        save_name=trim(adjustl(directory))//'/'//trim(adjustl(sysname))//'_current.data'
        open(tmp%ifn,file=save_name)
        select case(unit_system)
        case('au','a.u.')
          write(tmp%ifn,'(A)') "# time[a.u.],  current(x,y,z)[a.u.]" 
        case('A_eV_fs')
          write(tmp%ifn,'(A)') "# time[fs],    current(x,y,z)[A/Ang.^2]" 
        end select
        do ii=1,nt_em
          write(tmp%ifn, '(E13.5)',advance="no")     tmp%time_lr(ii)*utime_from_au
          write(tmp%ifn, '(3E16.6e3)',advance="yes") -tmp%curr_lr(ii,:)*tmp%uAperm_from_au/ulength_from_au
        end do
        close(tmp%ifn)
      end if
      
      !output lr data
      call eh_fourier(nt_em,nenergy,grid%dt,de,tmp%time_lr,tmp%curr_lr(:,1),tmp%fr_lr(:,1),tmp%fi_lr(:,1))
      call eh_fourier(nt_em,nenergy,grid%dt,de,tmp%time_lr,tmp%curr_lr(:,2),tmp%fr_lr(:,2),tmp%fi_lr(:,2))
      call eh_fourier(nt_em,nenergy,grid%dt,de,tmp%time_lr,tmp%curr_lr(:,3),tmp%fr_lr(:,3),tmp%fi_lr(:,3))
      if(comm_is_root(nproc_id_global)) then
        save_name=trim(adjustl(directory))//'/'//trim(adjustl(sysname))//'_lr.data'
        open(tmp%ifn,file=save_name)
        select case(unit_system)
        case('au','a.u.')
          write(tmp%ifn,'(A)') "# energy[a.u.], Re[epsilon](x,y,z), Im[epsilon](x,y,z)"
        case('A_eV_fs')
          write(tmp%ifn,'(A)') "# energy[eV], Re[epsilon](x,y,z), Im[epsilon](x,y,z)"
        end select
        do ii=1,nenergy
          write(tmp%ifn, '(E13.5)',advance="no")     dble(ii)*de*uenergy_from_au
          write(tmp%ifn, '(3E16.6e3)',advance="no")  1.0d0-4.0d0*pi*tmp%fi_lr(ii,:)/(-e_impulse)/(dble(ii)*de)
          write(tmp%ifn, '(3E16.6e3)',advance="yes") 4.0d0*pi*tmp%fr_lr(ii,:)/(-e_impulse)/(dble(ii)*de)
        end do
      end if
    end if
  end if
  
  !observation
  if(iobs_num_em>0) then
    if(comm_is_root(nproc_id_global)) then
      !make information file
      open(tmp%ifn,file=trim(directory)//"/obs0_info.data")
      write(tmp%ifn,'(A,A14)')                      'unit_system   =',trim(unit_system)
      write(tmp%ifn,'(A,I14)')                      'iperiodic     =',iperiodic
      write(tmp%ifn,'(A,ES14.5)')                   'dt_em         =',grid%dt*utime_from_au
      write(tmp%ifn,'(A,I14)')                      'nt_em         =',(tmp%iter_end-tmp%iter_sta+1)
      write(tmp%ifn,'(A,ES14.5,A,ES14.5,A,ES14.5)') 'al_em         =',&
            grid%rlsize(1)*ulength_from_au,', ',&
            grid%rlsize(2)*ulength_from_au,', ',&
            grid%rlsize(3)*ulength_from_au
      write(tmp%ifn,'(A,ES14.5,A,ES14.5,A,ES14.5)') 'dl_em         =',&
            grid%hgs(1)*ulength_from_au,', ',&
            grid%hgs(2)*ulength_from_au,', ',&
            grid%hgs(3)*ulength_from_au
      write(tmp%ifn,'(A,I14,A,I14,A,I14)')          'lg_sta        =',&
            grid%lg_sta(1),', ',grid%lg_sta(2),', ',grid%lg_sta(3)
      write(tmp%ifn,'(A,I14,A,I14,A,I14)')          'lg_end        =',&
            grid%lg_end(1),', ',grid%lg_end(2),', ',grid%lg_end(3)
      write(tmp%ifn,'(A,I14)')                      'iobs_num_em   =',iobs_num_em
      write(tmp%ifn,'(A,I14)')                      'iobs_samp_em  =',iobs_samp_em
      write(tmp%ifn,'(A,ES14.5)')                   'e_max         =',tmp%e_max
      write(tmp%ifn,'(A,ES14.5)')                   'h_max         =',tmp%h_max
      close(tmp%ifn)
    end if
  end if
  
  !deallocate
  deallocate(tmp%ex_y,tmp%c1_ex_y,tmp%c2_ex_y,tmp%ex_z,tmp%c1_ex_z,tmp%c2_ex_z,&
             tmp%ey_z,tmp%c1_ey_z,tmp%c2_ey_z,tmp%ey_x,tmp%c1_ey_x,tmp%c2_ey_x,&
             tmp%ez_x,tmp%c1_ez_x,tmp%c2_ez_x,tmp%ez_y,tmp%c1_ez_y,tmp%c2_ez_y,&
             tmp%hx_y,tmp%c1_hx_y,tmp%c2_hx_y,tmp%hx_z,tmp%c1_hx_z,tmp%c2_hx_z,&
             tmp%hy_z,tmp%c1_hy_z,tmp%c2_hy_z,tmp%hy_x,tmp%c1_hy_x,tmp%c2_hy_x,&
             tmp%hz_x,tmp%c1_hz_x,tmp%c2_hz_x,tmp%hz_y,tmp%c1_hz_y,tmp%c2_hz_y)
  
  !write end
  if(comm_is_root(nproc_id_global)) then
    write(*,*) "-------------------------------------------------------"
    write(*,*) "**************************"
    write(*,*) "FDTD end"
    write(*,*) "**************************"
  end if
  
end subroutine eh_finalize

!=========================================================================================
!= Fourier transformation in eh ==========================================================
subroutine eh_fourier(nt,ne,dt,de,ti,ft,fr,fi)
  use inputoutput, only: wf_em
  implicit none
  integer,intent(in)   :: nt,ne
  real(8),intent(in)   :: dt,de
  real(8),intent(in)   :: ti(nt),ft(nt)
  real(8),intent(out)  :: fr(0:ne),fi(0:ne)
  integer              :: ie,it
  real(8)              :: ft_wf(nt)
  real(8)              :: hw
  complex(8),parameter :: zi=(0.d0,1.d0)
  complex(8)           :: zf
  
  !apply window function
  if(wf_em=='y') then
    do it=1,nt
      ft_wf(it)=ft(it)*( 1.0d0 -3.0d0*(ti(it)/maxval(ti(:)))**2.0d0 +2.0d0*(ti(it)/maxval(ti(:)))**3.0d0 )
    end do
  else
    ft_wf(:)=ft(:)
  end if
  
  !Fourier transformation
  do ie=0,ne
    hw=dble(ie)*de; zf=(0.0d0,0.0d0);
!$omp parallel
!$omp do private(it) reduction( + : zf )
    do it=1,nt
      zf=zf+exp(zi*hw*ti(it))*ft_wf(it)
    end do
!$omp end do
!$omp end parallel
    zf=zf*dt; fr(ie)=real(zf,8); fi(ie)=aimag(zf)
  end do
  
end subroutine eh_fourier
