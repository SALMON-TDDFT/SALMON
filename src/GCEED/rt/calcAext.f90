!
!  Copyright 2017-2019 SALMON developers
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

subroutine calcAext
!$ use omp_lib
use salmon_parallel, only: nproc_id_global
use salmon_communication, only: comm_is_root, comm_summation
use scf_data
implicit none
integer :: j2
integer :: ipulse
real(8) :: env_trigon_1
real(8) :: env_trigon_2

select case(ikind_eext)
  case(0) 
    A_ext(1:3,:)=0.d0
    A_ext(1,:)=Fst*epdir_re1(1)
    A_ext(2,:)=Fst*epdir_re1(2)
    A_ext(3,:)=Fst*epdir_re1(3)
  case(1,15)
    A_ext(1:3,:)=0.d0
    do itt=Miter_rt+1,itotNtime+1
      if(dt*dble(itt) <= pulse_tw1)then
        ipulse=1
        call calc_env_trigon(ipulse,env_trigon_1)
        do j2=1,3
          A_ext(j2,itt)=A_ext(j2,itt)+amplitude1*epdir_re1(j2)*env_trigon_1
        end do
      end if
      if(dt*dble(itt)-t1_t2 >= 1.d-12 .and. dt*dble(itt)-t1_t2 <= pulse_tw2)then
        ipulse=2
        call calc_env_trigon(ipulse,env_trigon_2)
        do j2=1,3
          A_ext(j2,itt)=A_ext(j2,itt)+amplitude1*epdir_re1(j2)*env_trigon_2
        enddo
      end if
    end do
end select 

return

end subroutine calcAext
