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
subroutine deallocate_sendrecv
use scf_data
implicit none

if(iSCFRT==1)then
  deallocate(srmatbox1_x_3d,srmatbox1_y_3d,srmatbox1_z_3d)
  deallocate(srmatbox2_x_3d,srmatbox2_y_3d,srmatbox2_z_3d)
  deallocate(srmatbox3_x_3d,srmatbox3_y_3d,srmatbox3_z_3d)
  deallocate(srmatbox4_x_3d,srmatbox4_y_3d,srmatbox4_z_3d)

  deallocate(scmatbox1_x_3d,scmatbox1_y_3d,scmatbox1_z_3d)
  deallocate(scmatbox2_x_3d,scmatbox2_y_3d,scmatbox2_z_3d)
  deallocate(scmatbox3_x_3d,scmatbox3_y_3d,scmatbox3_z_3d)
  deallocate(scmatbox4_x_3d,scmatbox4_y_3d,scmatbox4_z_3d)
end if

if(iSCFRT==1.and.icalcforce==1)then
  deallocate(srmatbox1_x_5d,srmatbox1_y_5d,srmatbox1_z_5d)
  deallocate(srmatbox2_x_5d,srmatbox2_y_5d,srmatbox2_z_5d)
  deallocate(srmatbox3_x_5d,srmatbox3_y_5d,srmatbox3_z_5d)
  deallocate(srmatbox4_x_5d,srmatbox4_y_5d,srmatbox4_z_5d)
end if

end subroutine deallocate_sendrecv

