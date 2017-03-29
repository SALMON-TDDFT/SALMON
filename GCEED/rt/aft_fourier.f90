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

program aft_fourier
implicit none
integer :: t,j,iene,Nenergy
integer,parameter :: Ntime=4000
integer :: ibox
real(8) :: rbox
real(8), allocatable :: curr(:,:)
real(8), allocatable :: curr2(:,:)
complex(8),parameter :: zi=(0.d0,1.d0)
real(8),parameter :: Pi=3.141592653589793d0
real(8),parameter :: Ry=13.6058d0
real(8) :: F
real(8) :: hw,t2,TT,dt,dE
complex(8) :: zalpha(1:3)

open(10,file="current.data")
open(11,file="aft-fourier.data")

allocate(curr(1:3,0:Ntime))
curr(1:3,0)=0.d0
allocate(curr2(1:3,0:Ntime))
curr2(1:3,0)=0.d0

read(10,*)
read(10,*)

do t=1,Ntime
  read(10,*) rbox,(curr(j,t),j=1,3)
end do

dt=1.935104d-3/2.418884d-2
dE=0.01d0/2.d0/Ry
F=0.0001741555d0/0.348310d0

TT = dt*Ntime
Nenergy=2000

!curr2(1:3,1:Ntime)=-curr(1:3,1:Ntime)/F
curr2(1:3,1:Ntime)=curr(1:3,1:Ntime)/F

do iene=1,Nenergy
  hw=iene*dE ; zalpha=(0.d0,0.d0)  ! [a.u.]
  do t=1,Ntime
    t2=t*dt ; zalpha(:)=zalpha(:)+exp(zi*hw*t2)*curr2(:,t) & !hw*t is unitless      
                      *(1-3*(t2/TT)**2+2*(t2/TT)**3)
  end do
!  zalpha(1:3)=1.d0/(1.d0+4.d0*Pi*zi*zalpha(1:3)/hw)
!  zalpha(1:3)=1.d0+4.d0*Pi*zi*zalpha(1:3)/hw
!  zalpha(1:3)=1.d0+4.d0*Pi*zi*zalpha(1:3)/hw/F
  write(11,'(f12.4,6e14.6)') hw*2.d0*Ry,(real(zalpha(j),8),j=1,3),(imag(zalpha(j)),j=1,3)
end do


end program aft_fourier
