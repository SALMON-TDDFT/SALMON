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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
module salmon_global
  implicit none

!Parameters for pseudo-potential
  integer, parameter :: maxMKI=10
  integer :: MI,MKI
  integer :: Lmax_ps(maxMKI),Lloc_ps(maxMKI)
  integer :: iZatom(maxMKI)   ! Charge of ion
  character(10) :: ps_format(maxMKI) !shinohara
  integer :: ipsfileform(maxMKI)   ! file format for pseudo potential
! List of pseudopotential file formats
  integer,parameter :: n_Yabana_Bertsch_psformat = 1 !.rps
  integer,parameter :: n_ABINIT_psformat = 2 ! .pspnc
  integer,parameter :: n_FHI_psformat = 3 ! .cpi
  integer,parameter :: n_ABINITFHI_psformat = 4 ! .fhi


end module salmon_global
