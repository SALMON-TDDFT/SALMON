! Copyright 2017 Katsuyuki Nobusada, Masashi Noda, Kazuya Ishimura, Kenji Iida, Maiku Yamaguchi
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!=======================================================================
!=================================================== LDA (Perdew-Zunger)

SUBROUTINE Exc_Cor_fast(trho,Ex,Ec)
!$ use omp_lib
use scf_data

implicit none
integer :: ix,iy,iz
real(8) :: trho(mg_sta(1):mg_end(1),  &
                mg_sta(2):mg_end(2),  &
                mg_sta(3):mg_end(3))

real(8) :: Ex(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
             mg_sta(3):mg_end(3))
real(8) :: Ec(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
             mg_sta(3):mg_end(3))

!$OMP parallel do
do iz=ng_sta(3),ng_end(3)
do iy=ng_sta(2),ng_end(2)
do ix=ng_sta(1),ng_end(1)
  if(rho(ix,iy,iz)>=0.d0)then
    trho(ix,iy,iz)=rho(ix,iy,iz)+1.d-20
  else if(rho(ix,iy,iz)<-1.d-5)then
    write(*,*) "rho has large negative value", ix,iy,iz,rho(ix,iy,iz)
    stop
  else
    trho(ix,iy,iz)=1.d-20
  end if
end do
end do
end do

call xc_LDA_fast(trho,Ex,Ec)

END SUBROUTINE Exc_Cor_fast

!======================================================================

SUBROUTINE xc_LDA_fast(trho,Ex,Ec)
!$ use omp_lib
use scf_data
use new_world_sub
implicit none

real(8),parameter :: gm=-0.1423d0
! letter s means small, letter l means large
real(8),parameter :: sbu1=1.0529d0, sbu2=0.3334d0
real(8),parameter :: sbp1=1.3981d0, sbp2=0.2611d0
real(8),parameter :: rlau=0.0311d0 , rlbu=-0.048d0
real(8),parameter :: rlcu=0.002d0 , rldu=-0.0116d0

real(8) :: Cx
real(8) :: rs
real(8) :: zeta
real(8) :: sf
real(8) :: dsf

integer :: ix,iy,iz
real(8) :: sum1

real(8) :: trho(mg_sta(1):mg_end(1),  &
                mg_sta(2):mg_end(2),  &
                mg_sta(3):mg_end(3))

real(8) :: Ex(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
             mg_sta(3):mg_end(3))
real(8) :: Ec(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
             mg_sta(3):mg_end(3))

!$OMP parallel do &
!$OMP private(Cx,rs,zeta,sf,dsf)
do iz=ng_sta(3),ng_end(3)
do iy=ng_sta(2),ng_end(2)
do ix=ng_sta(1),ng_end(1)

  rs=(3.d0/(4*Pi*trho(ix,iy,iz)))**(1.d0/3.d0)

! Exchange

  Cx=3.d0/4.d0*(3.d0/(2*Pi))**(2.d0/3.d0)

  if(ilsda == 0)then
    Ex(ix,iy,iz)=-Cx/rs
  end if

! Correlation
! J. P. Perdew and A. Zunger, Phys. Rev. B, vol. 23, 5048 (1981).

  if(ilsda == 0)then
    if ( rs > 1.d0 ) then
      Ec(ix,iy,iz)=gm/( 1.d0 + sbu1*sqrt(rs) + sbu2*rs )
      Vxc(ix,iy,iz)=-4.d0/3.d0*Cx/rs+gm*( 1.d0 + 7.d0/6.d0*sbu1*sqrt(rs) + 4.d0/3.d0*sbu2*rs ) &
                     /( 1.d0 + sbu1*sqrt(rs) + sbu2*rs )**2
    else 
      Ec(ix,iy,iz)=rlau*log(rs) + rlbu   &
                     + rlcu*rs*log(rs) + rldu*rs
      Vxc(ix,iy,iz)=-4.d0/3.d0*Cx/rs + rlau*log(rs) + (rlbu-rlau/3.d0)   &
                     + 2.d0/3.d0*rlcu*rs*log(rs)      &
                     + (2.d0*rldu-rlcu)/3.d0*rs
    end if
  end if

end do
end do
end do

sum1=0.d0
!$omp parallel do reduction(+ : sum1)
do iz=ng_sta(3),ng_end(3)
do iy=ng_sta(2),ng_end(2)
do ix=ng_sta(1),ng_end(1)
  sum1=sum1+trho(ix,iy,iz)*(Ex(ix,iy,iz)+Ec(ix,iy,iz))
end do
end do
end do
sum1=sum1*Hvol
call MPI_ALLREDUCE(sum1,Exc,1,MPI_DOUBLE_PRECISION,      &
                  MPI_SUM,newworld_comm_h,IERR)

END SUBROUTINE xc_LDA_fast

!======================================================================
