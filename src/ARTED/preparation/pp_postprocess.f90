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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine pp_postprocess
  use global_variables, only: pp,Mass,NRps,NE,Rps,Lref,Mlps,Zps,Lmax,anorm,inorm,Nrmax,rad,radnl, &
                              vloctbl,dvloctbl,udVtbl,dudVtbl,rho_nlcc_tbl,tau_nlcc_tbl,NRloc,Rloc,flag_nlcc
  implicit none

  Mass(1:NE)=pp%rmass(1:NE)

  Lref(1:NE)=pp%lref(1:NE)
  NRps(1:NE)=pp%nrps(1:NE)
  Rps(1:NE)=pp%rps(1:NE)
  Mlps(1:NE)=pp%mlps(1:NE)
  Zps(1:NE)=pp%zps(1:NE)
  NRloc(1:NE)=pp%nrloc(1:NE)
  Rloc(1:NE)=pp%rloc(1:NE)

  anorm(0:Lmax,1:NE)=pp%anorm(0:Lmax,1:NE)
  inorm(0:Lmax,1:NE)=pp%inorm(0:Lmax,1:NE)
  
  rad(1:NRmax,1:NE)=pp%rad(1:Nrmax,1:NE)
  radnl(1:NRmax,1:NE)=pp%radnl(1:Nrmax,1:NE)
  
  vloctbl(1:Nrmax,1:NE)=pp%vloctbl(1:Nrmax,1:NE)
  dvloctbl(1:Nrmax,1:NE)=pp%dvloctbl(1:Nrmax,1:NE)
  udvtbl(1:Nrmax,0:Lmax,1:NE)=pp%udvtbl(1:Nrmax,0:Lmax,1:NE)
  dudvtbl(1:Nrmax,0:Lmax,1:NE)=pp%dudvtbl(1:Nrmax,0:Lmax,1:NE)
  
  rho_nlcc_tbl(1:Nrmax,1:NE)=pp%rho_nlcc_tbl(1:Nrmax,1:NE)
  tau_nlcc_tbl(1:Nrmax,1:NE)=pp%tau_nlcc_tbl(1:Nrmax,1:NE)

  flag_nlcc=pp%flag_nlcc

end subroutine pp_postprocess
