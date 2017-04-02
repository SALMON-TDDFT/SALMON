!
!  Copyright 2016 ARTED developers
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
module environment
  implicit none

  real(8) :: CPU_TASK_RATIO
  real(8) :: MIC_TASK_RATIO

  integer :: CPU_PROCESS_PER_NODE
  integer :: MIC_PROCESS_PER_NODE

  integer :: ENABLE_LOAD_BALANCER

contains
  subroutine load_environments
    implicit none

    call get_cpu_task_ratio_internal(CPU_TASK_RATIO)
    MIC_TASK_RATIO = 2.0 - CPU_TASK_RATIO

    call get_cpu_ppn_internal(CPU_PROCESS_PER_NODE)
    call get_mic_ppn_internal(MIC_PROCESS_PER_NODE)

    call get_load_balancer_flag_internal(ENABLE_LOAD_BALANCER)
  end subroutine
end module environment
