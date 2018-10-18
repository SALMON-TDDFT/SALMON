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
  use inputoutput,          only: utime_from_au,ulength_from_au,unit_system,&
                                  directory,iobs_num_em,iobs_samp_em
  use salmon_parallel,      only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use salmon_maxwell, only:fdtd_grid,fdtd_tmp
  implicit none
  type(fdtd_grid)  :: grid
  type(fdtd_tmp)   :: tmp
  
  !observation
  if(iobs_num_em>0) then
    if(comm_is_root(nproc_id_global)) then
      !make information file
      open(tmp%ifn,file=trim(directory)//"obs0_info.data")
      write(tmp%ifn,'(A,A14)')                      'unit_system  =',trim(unit_system)
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
      write(tmp%ifn,'(A,I14,A,I14,A,I14)')          'lg_sta       =',&
            grid%lg_sta(1),', ',grid%lg_sta(2),', ',grid%lg_sta(3)
      write(tmp%ifn,'(A,I14,A,I14,A,I14)')          'lg_end       =',&
            grid%lg_end(1),', ',grid%lg_end(2),', ',grid%lg_end(3)
      write(tmp%ifn,'(A,I14)')                      'iobs_num_em  =',iobs_num_em
      write(tmp%ifn,'(A,I14)')                      'iobs_samp_em =',iobs_samp_em
      write(tmp%ifn,'(A,ES14.5)')                   'e_max        =',tmp%e_max
      write(tmp%ifn,'(A,ES14.5)')                   'h_max        =',tmp%h_max
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
