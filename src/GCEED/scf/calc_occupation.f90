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
subroutine calc_occupation
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use scf_data
  implicit none
  real(8),parameter :: bltz=8.6173214d-5
  real(8) :: factor
  integer :: ii,iob,p1,p2,p5,ik
  integer :: is
  integer :: is_sta,is_end
  integer :: iobsta(2),iobend(2)
  real(8) :: temperature_au

  if(ilsda==0)then
    is_sta=1
    is_end=1
  else
    is_sta=1
    is_end=2
    iobsta(1)=1
    iobend(1)=MST(1)
    iobsta(2)=MST(1)+1
    iobend(2)=itotMST
  end if

  rocc(1:itotMST,:num_kpoints_rd)=0.d0
  if(ilsda==0)then
    rocc(1:ifMST(1),:num_kpoints_rd)=2.d0
  else
    rocc(1:ifMST(1),:num_kpoints_rd)=1.d0
    rocc(MST(1)+1:MST(1)+ifMST(2),:num_kpoints_rd)=1.d0
  end if

end subroutine calc_occupation
SUBROUTINE ne2mu !(nein,muout)
use salmon_parallel, only: nproc_id_global
use salmon_communication, only: comm_is_root
use scf_data
implicit none

!      real(8), intent(in)  :: nein
!      real(8), intent(out) :: muout
real(8) :: nein,muout
real(8) :: mu1,mu2,mu3,ne1,ne3,ne3o,diff_ne
integer :: ii,p5,p1,p2,iob,ik
integer :: nc
real(8) :: diff_mu,muo
real(8) :: diff_ne2

!      diff_ne = 100.d0

      nein=rNetot

      mu1 = esp(1,1)       !-100.d0
      mu2 = esp(itotMST,1)
      nc=0

200   continue
      diff_ne = 100.d0
      diff_ne2 = 100.d0
!      mu2 = 0.4d0   !esp(itotMST)       !100.d0

      ii = 1

      ne1  = 0d0
      ne3  = 0d0
      ne3o = 0d0
!
      diff_mu =100.d0
      muo=0.0d0
!
      call mu2ne(mu1,ne1)

!      do while ( diff_ne > 1d-10 .and. ii < 1000 )
      do ii=1,1000
         if ( ii .eq. 1000 ) then
            if ( nc .le. 50)  then
!               print *,'=================================='
!               print *,'  Start next const Ne iterration  '
!               print *,'=================================='
               nc= nc + 1
               mu1 = esp(1,1) - 0.2d0*dble(nc)
               mu2 = esp(itotMST,1) + 0.2d0*dble(nc)
               go to 200
            else
               print *,'=================================='
               print *,'Const Ne does not converged!!!!!!!'
               print *,'=================================='
               goto 100
            endif
         endif
         if ( diff_ne < 1d-10  .and. diff_mu < 1d-9  &
             .and. diff_ne2 < 1.d-9 ) goto 100
         mu3 = mu1 + (mu2-mu1)/2.d0
!         if (ii.eq.1) mu3=mu2
         ne3=0
         call mu2ne(mu3,ne3)
!         if (ii.eq.1) call mu2ne(mu2,ne3)
         diff_ne = abs((ne3-ne3o)/ne3)
         diff_ne2 = abs(nein-ne3)
!         print *, ii,'diff_ne2=',diff_ne2,diff_ne,mu3,ne1,ne3
         if ( (ne1-nein)*(ne3-nein) > 0 ) then
            mu1=mu3
            ne1=ne3
         else
            mu2=mu3
         end if

!         ii = ii+1
         ne3o = ne3
         diff_mu = mu3 - muo
         muo  = mu3

      end do
100   continue
!
      muout = mu3

      call mu2ne(muout,ne3)

      if(comm_is_root(nproc_id_global))then       !.and.iscf==1) then
         write(*,*)
         write(*,*) 'Fractional Occupation Numbers :'
         write(*,*)
         do ik=1,num_kpoints_rd
           do p5=1,(itotMST+4)/5
             p1=5*(p5-1)+1
             p2=5*p5 ; if ( p2 > itotMST ) p2=itotMST
             write(*,'(1x,5(i3,f8.4,1x))')  (iob,rocc(iob,ik),iob=p1,p2)
           end do
         end do
         write(*,*)
         write(*,'(a,f15.8)') ' Fermi level         = ',muout*2d0*Ry
         write(*,'(a,f15.8)') ' Number of Electrons = ',ne3
         write(*,*)
      end if
!
END SUBROUTINE ne2mu
!
!
SUBROUTINE mu2ne(muin,neout)
use scf_data
!use global_variables
implicit none

real(8),parameter :: bltz=8.6173214d-5
real(8), intent(in)  :: muin
real(8), intent(out) :: neout
real(8) :: fact
integer :: iob
integer :: ik
real(8) :: temperature_au

temperature_au=temperature_k*bltz/2.d0/Ry

neout=0d0
do ik=1,num_kpoints_rd
  do iob=1,itotMST
    fact=(dble(esp(iob,ik))-muin)/temperature_au
    if(fact.ge.40.d0) then
       rocc(iob,ik)=0.d0
    else
       rocc(iob,ik)=2d0/(1.d0+dexp((dble(esp(iob,ik))-muin)/temperature_au))
    endif
    neout=neout+rocc(iob,ik)
  end do
end do


END SUBROUTINE mu2ne

!
SUBROUTINE ne2mu_p
use salmon_parallel, only: nproc_id_global
use salmon_communication, only: comm_is_root
use scf_data
implicit none

real(8) :: nein,muout
real(8) :: mu1,mu2,mu3,ne1,ne3,ne3o,diff_ne
integer :: ii,p5,p1,p2,iob,iik
integer :: nc
real(8) :: diff_mu,muo
real(8) :: diff_ne2

!      diff_ne = 100.d0
      nein=rNetot

      mu1 = esp(1,1)-0.3d0
      mu2 = esp(itotMST,1)+0.3d0
      nc=0

200   continue
      diff_ne = 100.d0
      diff_ne2 = 100.d0
!      mu2 = 0.4d0   !esp(itotMST)       !100.d0

      ii = 1

      ne1  = 0d0
      ne3  = 0d0
      ne3o = 0d0
!
      diff_mu =100.d0
      muo=0.0d0
!
      call mu2ne_p(mu1,ne1)

!      do while ( diff_ne > 1d-10 .and. ii < 1000 )
      do ii=1,1000
         if ( ii .eq. 1000 ) then
            if ( nc .le. 50)  then
!               print *,'=================================='
!               print *,'  Start next const Ne iterration  '
!               print *,'=================================='
               nc= nc + 1
               mu1 = esp(1,1) - 0.2d0*dble(nc)
               mu2 = esp(itotMST,1) + 0.2d0*dble(nc)
               go to 200
            else
               print *,'=================================='
               print *,'Const Ne does not converged!!!!!!!'
               print *,'=================================='
               goto 100
            endif
         endif
         if ( diff_ne < 1d-10  .and. diff_mu < 1d-9  &
             .and. diff_ne2 < 1.d-9 ) goto 100
         mu3 = mu1 + (mu2-mu1)/2.d0
!         if (ii.eq.1) mu3=mu2
         ne3=0
         call mu2ne_p(mu3,ne3)
!         if (ii.eq.1) call mu2ne(mu2,ne3)
         diff_ne = abs((ne3-ne3o)/ne3)
         diff_ne2 = abs(nein-ne3)
!         print *, ii,'diff_ne2=',diff_ne2,diff_ne,mu3,ne1,ne3
         if ( (ne1-nein)*(ne3-nein) > 0 ) then
            mu1=mu3
            ne1=ne3
         else
            mu2=mu3
         end if

!         ii = ii+1
         ne3o = ne3
         diff_mu = mu3 - muo
         muo  = mu3

      end do
100   continue
!
      muout = mu3

      call mu2ne_p(muout,ne3)

      if(comm_is_root(nproc_id_global))then
         write(*,*)
         write(*,*) 'Fractional Occupation Numbers :'
         write(*,*)
         do iik=1,num_kpoints_rd
            if(iik<=10)then
              print *, ' iik = ',iik
              do p5=1,(itotMST+4)/5
                 p1=5*(p5-1)+1
                 p2=5*p5 ; if ( p2 > itotMST ) p2=itotMST
                 write(*,'(1x,5(i3,f8.4,1x))')  (iob,rocc(iob,iik),iob=p1,p2)
              end do
            endif
         enddo
         write(*,*)
         write(*,'(a,f15.8)') ' Fermi level         = ',muout*2d0*Ry
         write(*,'(a,f15.8)') ' Number of Electrons = ',ne3
         write(*,*)
      end if
!
END SUBROUTINE ne2mu_p
!
SUBROUTINE mu2ne_p(muin,neout)
use scf_data
implicit none

real(8),parameter :: bltz=8.6173214d-5
real(8), intent(in)  :: muin
real(8), intent(out) :: neout
real(8) :: fact
integer :: iob,iik
real(8) :: temperature_au

temperature_au=temperature_k*bltz/2.d0/Ry

neout=0d0

do iik=1,num_kpoints_rd
do iob=1,itotMST
   fact=(dble(esp(iob,iik))-muin)/temperature_au
   if(fact.ge.40.d0) then
      rocc(iob,iik)=0.d0
   else
      rocc(iob,iik)=2d0/(1.d0+dexp((dble(esp(iob,iik))-muin)/temperature_au))
   endif
   neout=neout+rocc(iob,iik)*wtk(iik)
end do
end do

END SUBROUTINE mu2ne_p



