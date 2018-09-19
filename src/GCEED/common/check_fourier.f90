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
!=======================================================================
subroutine check_fourier
  use salmon_global, only: fourier
  use salmon_parallel, only: nproc_id_global, end_parallel
  use salmon_communication, only: comm_is_root
  use scf_data
  implicit none
  integer :: j,ii
  integer :: lg_num_tmp

  if(iflag_hartree==4)then
    ! this code treats the situation that lg_num(1:3) is less than or equal to 48,828,125
    do j=1,3
      lg_num_tmp=lg_num(j)
      do ii=1,26
        if(mod(lg_num_tmp,2)==0)then
          lg_num_tmp=lg_num_tmp/2
        end if
      end do
    
      do ii=1,17
        if(mod(lg_num_tmp,3)==0)then
          lg_num_tmp=lg_num_tmp/3
        end if
      end do
    
      do ii=1,11
        if(mod(lg_num_tmp,5)==0)then
          lg_num_tmp=lg_num_tmp/5
        end if
      end do
  
      if(lg_num_tmp/=1)then
        if (comm_is_root(nproc_id_global)) then
          write(*,*) "When using FFTE, prime factors for number of grids must be combination of 2, 3 or 5."
        end if
        call end_parallel
        stop
      end if
    end do 
  end if
  
end subroutine check_fourier
