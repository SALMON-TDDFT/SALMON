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
module sendrecv_sub

  use scf_data
  use new_world_sub
  use init_sendrecv_sub
  
  interface sendrecv
  
    module procedure R_sendrecv, C_sendrecv
  
  end interface 

contains

!==================================================================================================

subroutine R_sendrecv(tpsi)
  use salmon_communication, only: comm_proc_null, comm_start_all, comm_wait_all
  use persistent_comm,      only: ireq => nreqs_rorbital, &
                                  nrange => nrange_groupob, nshape => nshape_orbital
  use pack_unpack
  implicit none
  real(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd,mg_sta(2)-Nd:mg_end(2)+Nd, &
                  mg_sta(3)-Nd:mg_end(3)+Nd)
  integer :: iup,idw,jup,jdw,kup,kdw

  iup=iup_array(1)
  idw=idw_array(1)
  jup=jup_array(1)
  jdw=jdw_array(1)
  kup=kup_array(1)
  kdw=kdw_array(1)

  !send from idw to iup
  if(iup/=comm_proc_null)then
    call pack_data(nshape(1:3), nrange(1:3,1), tpsi, srmatbox1_x_3d)
  end if
  call comm_start_all(ireq(1:2))

  !send from iup to idw
  if(idw/=comm_proc_null)then
    call pack_data(nshape(1:3), nrange(1:3,2), tpsi, srmatbox3_x_3d)
  end if
  call comm_start_all(ireq(3:4))

  !send from jdw to jup
  if(jup/=comm_proc_null)then
    call pack_data(nshape(1:3), nrange(1:3,3), tpsi, srmatbox1_y_3d)
  end if
  call comm_start_all(ireq(5:6))

  !send from jup to jdw
  if(jdw/=comm_proc_null)then
    call pack_data(nshape(1:3), nrange(1:3,4), tpsi, srmatbox3_y_3d)
  end if
  call comm_start_all(ireq(7:8))

  !send from kdw to kup
  if(kup/=comm_proc_null)then
    call pack_data(nshape(1:3), nrange(1:3,5), tpsi, srmatbox1_z_3d)
  end if
  call comm_start_all(ireq(9:10))

  !send from kup to kdw
  if(kdw/=comm_proc_null)then
    call pack_data(nshape(1:3), nrange(1:3,6), tpsi, srmatbox3_z_3d)
  end if
  call comm_start_all(ireq(11:12))


  call comm_wait_all(ireq(1:2))
  if(idw/=comm_proc_null)then
    call unpack_data(nshape(1:3), nrange(1:3,7), srmatbox2_x_3d, tpsi)
  end if

  call comm_wait_all(ireq(3:4))
  if(iup/=comm_proc_null)then
    call unpack_data(nshape(1:3), nrange(1:3,8), srmatbox4_x_3d, tpsi)
  end if

  call comm_wait_all(ireq(5:6))
  if(jdw/=comm_proc_null)then
    call unpack_data(nshape(1:3), nrange(1:3,9), srmatbox2_y_3d, tpsi)
  end if

  call comm_wait_all(ireq(7:8))
  if(jup/=comm_proc_null)then
    call unpack_data(nshape(1:3), nrange(1:3,10), srmatbox4_y_3d, tpsi)
  end if

  call comm_wait_all(ireq(9:10))
  if(kdw/=comm_proc_null)then
    call unpack_data(nshape(1:3), nrange(1:3,11), srmatbox2_z_3d, tpsi)
  end if

  call comm_wait_all(ireq(11:12))
  if(kup/=comm_proc_null)then
    call unpack_data(nshape(1:3), nrange(1:3,12), srmatbox4_z_3d, tpsi)
  end if

end subroutine R_sendrecv

!==================================================================================================

subroutine C_sendrecv(tpsi)
  use salmon_communication, only: comm_proc_null, comm_start_all, comm_wait_all
  use persistent_comm,      only: ireq => nreqs_corbital, &
                                  nrange => nrange_groupob, nshape => nshape_groupob
  use pack_unpack
  implicit none
  complex(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd,mg_sta(2)-Nd:mg_end(2)+Nd, &
                     mg_sta(3)-Nd:mg_end(3)+Nd)
  integer :: iup,idw,jup,jdw,kup,kdw
  
  iup=iup_array(1)
  idw=idw_array(1)
  jup=jup_array(1)
  jdw=jdw_array(1)
  kup=kup_array(1)
  kdw=kdw_array(1)
  
  !send from idw to iup
  if(iup/=comm_proc_null)then
    call pack_data(nshape(1:3), nrange(1:3,1), tpsi, scmatbox1_x_3d)
  end if
  call comm_start_all(ireq(1:2))

  !send from iup to idw
  if(idw/=comm_proc_null)then
    call pack_data(nshape(1:3), nrange(1:3,2), tpsi, scmatbox3_x_3d)
  end if
  call comm_start_all(ireq(3:4))

  !send from jdw to jup
  if(jup/=comm_proc_null)then
    call pack_data(nshape(1:3), nrange(1:3,3), tpsi, scmatbox1_y_3d)
  end if
  call comm_start_all(ireq(5:6))

  !send from jup to jdw
  if(jdw/=comm_proc_null)then
    call pack_data(nshape(1:3), nrange(1:3,4), tpsi, scmatbox3_y_3d)
  end if
  call comm_start_all(ireq(7:8))

  !send from kdw to kup
  if(kup/=comm_proc_null)then
    call pack_data(nshape(1:3), nrange(1:3,5), tpsi, scmatbox1_z_3d)
  end if
  call comm_start_all(ireq(9:10))

  !send from kup to kdw
  if(kdw/=comm_proc_null)then
    call pack_data(nshape(1:3), nrange(1:3,6), tpsi, scmatbox3_z_3d)
  end if
  call comm_start_all(ireq(11:12))


  call comm_wait_all(ireq(1:2))
  if(idw/=comm_proc_null)then
    call unpack_data(nshape(1:3), nrange(1:3,7), scmatbox2_x_3d, tpsi)
  end if

  call comm_wait_all(ireq(3:4))
  if(iup/=comm_proc_null)then
    call unpack_data(nshape(1:3), nrange(1:3,8), scmatbox4_x_3d, tpsi)
  end if

  call comm_wait_all(ireq(5:6))
  if(jdw/=comm_proc_null)then
    call unpack_data(nshape(1:3), nrange(1:3,9), scmatbox2_y_3d, tpsi)
  end if

  call comm_wait_all(ireq(7:8))
  if(jup/=comm_proc_null)then
    call unpack_data(nshape(1:3), nrange(1:3,10), scmatbox4_y_3d, tpsi)
  end if

  call comm_wait_all(ireq(9:10))
  if(kdw/=comm_proc_null)then
    call unpack_data(nshape(1:3), nrange(1:3,11), scmatbox2_z_3d, tpsi)
  end if

  call comm_wait_all(ireq(11:12))
  if(kup/=comm_proc_null)then
    call unpack_data(nshape(1:3), nrange(1:3,12), scmatbox4_z_3d, tpsi)
  end if
  
end subroutine C_sendrecv

!==================================================================================================

end module sendrecv_sub
