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
Subroutine Fermi_Dirac_distribution
  use Global_Variables
  use communication
  use misc_routines, only: get_wtime
  implicit none
  real(8) :: chemical_potential
  real(8) :: chem_max,chem_min,elec_num
  real(8) :: beta_FD
  real(8),allocatable :: occ_l(:,:),esp_l(:,:),esp_temp(:,:)
  real(8) s,st
  integer :: ik,ib
  real(8) :: timer1,timer2
  
  timer1=get_wtime()
  beta_FD=1d0/(KbTev/(2d0*Ry))
  allocate(occ_l(NB,NK),esp_l(NB,NK),esp_temp(NB,NK))
  occ_l=0d0 ; esp_l=0d0
  esp_l(:,NK_s:NK_e)=esp(:,NK_s:NK_e)
  call comm_summation(esp_l,esp_temp,NB*NK,proc_group(2))
  chem_max=maxval(esp_temp)
  chem_min=minval(esp_temp)
  if(comm_is_root())then
    write(*,*)'max esp =',chem_max
    write(*,*)'min esp =',chem_min
  end if
  chemical_potential=0.5d0*(chem_max+chem_min)
  do 

    do ik=NK_s,NK_e
      do ib=1,NB
        occ_l(ib,ik)=2.d0/(NKxyz)*wk(ik) &
          &/(exp(beta_FD*(esp(ib,ik)-chemical_potential))+1d0)
      end do
    end do

    st=sum(occ_l(:,NK_s:NK_e))
    call comm_summation(st,s,proc_group(2))
    elec_num=s

    if(abs(elec_num-dble(Nelec)) <= 1d-6)exit
    if(elec_num-dble(Nelec) > 0d0)then
      chem_max=chemical_potential
      chemical_potential=0.5d0*(chem_max+chem_min)
    else
      chem_min=chemical_potential
      chemical_potential=0.5d0*(chem_max+chem_min)
    end if
    if(chem_max-chem_min < 1d-10)exit

  end do

  call comm_summation(occ_l,occ,NB*NK,proc_group(2))
  st=sum(occ_l(Nelec/2+1:NB,NK_s:NK_e))
  call comm_summation(st,s,proc_group(2))

  timer2=get_wtime()
  if(comm_is_root())then
    write(*,*)'Fermi-Dirac dist. time=',timer2-timer1,'sec'
    write(*,*)'chemical potential =',chemical_potential
    write(*,*)'elec_num =',elec_num
    write(*,*)'excited electron =',s
!    open(126,file='occ.dat')
!    do ik=1,NK
!      do ib=1,NB
!        write(126,*)ik,ib,occ(ib,ik)
!      end do
!    end do
!    close(126)
  end if

  return
end Subroutine Fermi_Dirac_distribution
