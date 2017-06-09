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
subroutine set_filename
  use salmon_global
  use inputoutput
  use scf_data
  implicit none
  
  file_IN=trim(sysname)//".data"
  file_OUT=trim(sysname)//".data"
  LDA_Info=trim(sysname)//".info"

  file_OUT_rt=trim(sysname)//"_rtreent.data"
  file_IN_rt=trim(sysname)//"_rtreent.data"

  file_RT=trim(sysname)//"-RT.data"
  file_alpha=trim(sysname)//"-ALP.data"

  file_RT_q=trim(sysname)//"-RT_q.data"
  file_alpha_q=trim(sysname)//"-ALP_q.data"
  file_RT_e=trim(sysname)//"-RT_e.data"
  file_RT_dip2=trim(sysname)//"-RT-dip2.data"
  file_alpha_dip2=trim(sysname)//"-ALP-dip2.data"
  file_RT_dip2_q=trim(sysname)//"-RT-dip2-q.data"
  file_alpha_dip2_q=trim(sysname)//"-ALP-dip2-q.data"
  file_RT_dip2_e=trim(sysname)//"-RT-dip2-e.data"
  file_external=trim(sysname)//"-ext.data"
  file_Projection=trim(sysname)//"-proj.data"
  
  file_ini=trim(sysname)//"_ini.data"

  return

end subroutine set_filename

