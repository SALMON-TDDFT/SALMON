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

! ifunc : +1:  calclation of esp / +0: no calculation of esp
!         +2:  2nd or larger loop (zpsi+htpsi) / +0: 1st loop (tpsi+htpsi)
!         +4:  calculation of rhobox / +0: no calculation of rhobox
!         -1:  special version for N_hamil=4 and nn=3
subroutine add_polynomial(tpsi,htpsi,tpsi_out,iobmax,nn,ifunc)
use scf_data
implicit none
complex(8) :: tpsi(iwk2sta(1):iwk2end(1)+1,iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3),   &
                   1:iobnum,1)
complex(8) :: htpsi(iwk2sta(1):iwk2end(1)+1,  &
                    iwk2sta(2):iwk2end(2),      &
                    iwk2sta(3):iwk2end(3),     &
                   1:iobnum,1)
complex(8) :: tpsi_out(iwk2sta(1):iwk2end(1)+1,iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3),   &
                   1:iobnum,1)
integer :: iobmax
integer :: ifunc
integer :: nn
integer :: iob,ix,iy,iz
complex(8) :: cbox
complex(8), parameter :: zi=(0.d0,1.d0)
integer :: iob_allob

if(ifunc==0)then
  do iob=1,iobmax
!$OMP parallel do
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
      tpsi_out(ix,iy,iz,iob,1)=tpsi(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
else if(ifunc==1)then
  do iob=1,iobmax
    cbox=0.d0
!$OMP parallel do reduction(+:cbox)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      cbox=cbox+conjg(tpsi(ix,iy,iz,iob,1))*htpsi(ix,iy,iz,iob,1)
      htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
      tpsi_out(ix,iy,iz,iob,1)=tpsi(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
    end do
    end do
    end do
    esp2(iob,1)=dble(cbox)*Hvol
  end do
else if(ifunc==2)then
  do iob=1,iobmax
!$OMP parallel do
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
      tpsi_out(ix,iy,iz,iob,1)=tpsi_out(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
else if(ifunc==3)then
  do iob=1,iobmax
    cbox=0.d0
!$OMP parallel do reduction(+:cbox)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      cbox=cbox+conjg(tpsi(ix,iy,iz,iob,1))*htpsi(ix,iy,iz,iob,1)
      htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
      tpsi_out(ix,iy,iz,iob,1)=tpsi_out(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
    end do
    end do
    end do
    esp2(iob,1)=dble(cbox)*Hvol
  end do
else if(ifunc==4)then
  if(ilsda==0)then
    do iob=1,iobmax
      call calc_allob(iob,iob_allob)
!$OMP parallel do
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
        tpsi_out(ix,iy,iz,iob,1)=tpsi(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
        rhobox(ix,iy,iz)=rhobox(ix,iy,iz)+tpsi_out(ix,iy,iz,iob,1)*conjg(tpsi_out(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
      end do
      end do
      end do
    end do
  else
    do iob=1,iobmax
      call calc_allob(iob,iob_allob)
      if(iob_allob<=MST(1))then
!$OMP parallel do
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
          tpsi_out(ix,iy,iz,iob,1)=tpsi(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
          rhobox_s(ix,iy,iz,1)=rhobox_s(ix,iy,iz,1)+tpsi_out(ix,iy,iz,iob,1)*conjg(tpsi_out(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
        end do
        end do
        end do
      else
!$OMP parallel do
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
          tpsi_out(ix,iy,iz,iob,1)=tpsi(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
          rhobox_s(ix,iy,iz,2)=rhobox_s(ix,iy,iz,2)+tpsi_out(ix,iy,iz,iob,1)*conjg(tpsi_out(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
        end do
        end do
        end do
      end if
    end do
  end if 
else if(ifunc==5)then
  if(ilsda==0)then
    do iob=1,iobmax
      call calc_allob(iob,iob_allob)
      cbox=0.d0
!$OMP parallel do reduction(+:cbox)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        cbox=cbox+conjg(tpsi(ix,iy,iz,iob,1))*htpsi(ix,iy,iz,iob,1)
        htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
        tpsi_out(ix,iy,iz,iob,1)=tpsi(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
        rhobox(ix,iy,iz)=rhobox(ix,iy,iz)+tpsi_out(ix,iy,iz,iob,1)*conjg(tpsi_out(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
      end do
      end do
      end do
      esp2(iob,1)=dble(cbox)*Hvol
    end do
  else if(ilsda==1)then
    do iob=1,iobmax
      call calc_allob(iob,iob_allob)
      cbox=0.d0
      if(iob_allob<=MST(1))then
!$OMP parallel do reduction(+:cbox)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cbox=cbox+conjg(tpsi(ix,iy,iz,iob,1))*htpsi(ix,iy,iz,iob,1)
          htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
          tpsi_out(ix,iy,iz,iob,1)=tpsi(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
          rhobox_s(ix,iy,iz,1)=rhobox_s(ix,iy,iz,1)+tpsi_out(ix,iy,iz,iob,1)*conjg(tpsi_out(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
        end do
        end do
        end do
      else
!$OMP parallel do reduction(+:cbox)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cbox=cbox+conjg(tpsi(ix,iy,iz,iob,1))*htpsi(ix,iy,iz,iob,1)
          htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
          tpsi_out(ix,iy,iz,iob,1)=tpsi(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
          rhobox_s(ix,iy,iz,2)=rhobox_s(ix,iy,iz,2)+tpsi_out(ix,iy,iz,iob,1)*conjg(tpsi_out(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
        end do
        end do
        end do
      end if
      esp2(iob,1)=dble(cbox)*Hvol
    end do
  end if
else if(ifunc==6)then
  if(ilsda==0)then
    do iob=1,iobmax
      call calc_allob(iob,iob_allob)
!$OMP parallel do
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
        tpsi_out(ix,iy,iz,iob,1)=tpsi_out(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
        rhobox(ix,iy,iz)=rhobox(ix,iy,iz)+tpsi_out(ix,iy,iz,iob,1)*conjg(tpsi_out(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
      end do
      end do
      end do
    end do
  else if(ilsda==1)then
    do iob=1,iobmax
      call calc_allob(iob,iob_allob)
      if(iob_allob<=MST(1))then
!$OMP parallel do
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
          tpsi_out(ix,iy,iz,iob,1)=tpsi_out(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
          rhobox_s(ix,iy,iz,1)=rhobox_s(ix,iy,iz,1)+tpsi_out(ix,iy,iz,iob,1)*conjg(tpsi_out(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
        end do
        end do
        end do
      else
!$OMP parallel do
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
          tpsi_out(ix,iy,iz,iob,1)=tpsi_out(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
          rhobox_s(ix,iy,iz,2)=rhobox_s(ix,iy,iz,2)+tpsi_out(ix,iy,iz,iob,1)*conjg(tpsi_out(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
        end do
        end do
        end do
      end if
    end do
  end if
else if(ifunc==7)then
  if(ilsda==0)then
    do iob=1,iobmax
      cbox=0.d0
!$OMP parallel do reduction(+:cbox)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        cbox=cbox+conjg(tpsi(ix,iy,iz,iob,1))*htpsi(ix,iy,iz,iob,1)
        htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
        tpsi_out(ix,iy,iz,iob,1)=tpsi_out(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
        rhobox(ix,iy,iz)=rhobox(ix,iy,iz)+tpsi_out(ix,iy,iz,iob,1)*conjg(tpsi_out(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
      end do
      end do
      end do
      esp2(iob,1)=dble(cbox)*Hvol
    end do
  else if(ilsda==1)then
    do iob=1,iobmax
      call calc_allob(iob,iob_allob)
      cbox=0.d0
      if(iob_allob<=MST(1))then
!$OMP parallel do reduction(+:cbox)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cbox=cbox+conjg(tpsi(ix,iy,iz,iob,1))*htpsi(ix,iy,iz,iob,1)
          htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
          tpsi_out(ix,iy,iz,iob,1)=tpsi_out(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
          rhobox_s(ix,iy,iz,1)=rhobox_s(ix,iy,iz,1)+tpsi_out(ix,iy,iz,iob,1)*conjg(tpsi_out(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
        end do
        end do
        end do
      else
!$OMP parallel do reduction(+:cbox)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cbox=cbox+conjg(tpsi(ix,iy,iz,iob,1))*htpsi(ix,iy,iz,iob,1)
          htpsi(ix,iy,iz,iob,1)=-zi*dt*htpsi(ix,iy,iz,iob,1)/dble(nn)
          tpsi_out(ix,iy,iz,iob,1)=tpsi_out(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
          rhobox_s(ix,iy,iz,2)=rhobox_s(ix,iy,iz,2)+tpsi_out(ix,iy,iz,iob,1)*conjg(tpsi_out(ix,iy,iz,iob,1))*rocc(iob_allob,1)*wtk(1)
        end do
        end do
        end do
      end if
      esp2(iob,1)=dble(cbox)*Hvol
    end do
  end if
else if(ifunc==-1)then
  do iob=1,iobmax
!$OMP parallel do
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      tpsi(ix,iy,iz,iob,1)=-zi*dt*tpsi(ix,iy,iz,iob,1)/dble(nn-1)
      htpsi(ix,iy,iz,iob,1)=(-zi*dt)**2*htpsi(ix,iy,iz,iob,1)/dble(nn-1)/dble(nn)
      tpsi_out(ix,iy,iz,iob,1)=tpsi_out(ix,iy,iz,iob,1)+tpsi(ix,iy,iz,iob,1)+htpsi(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
end if

end subroutine add_polynomial
