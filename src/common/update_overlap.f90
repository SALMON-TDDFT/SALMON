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
module update_overlap_sub

contains

!===================================================================================================================================

subroutine update_overlap_R(tpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb,Nd &
                     ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,irank_overlap,icomm)
  use salmon_communication, only: comm_proc_null, comm_isend, comm_irecv, comm_wait_all
  implicit none
  integer,intent(in) :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb,Nd &
                       ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,irank_overlap(6),icomm
  real(8) :: tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb)
  !
  integer :: ix,iy,iz,iorb
  integer :: iup,idw,jup,jdw,kup,kdw
  integer :: ireq(12)

  real(8),allocatable :: commbuf_x(:,:,:,:,:),commbuf_y(:,:,:,:,:),commbuf_z(:,:,:,:,:)

  iup = irank_overlap(1)
  idw = irank_overlap(2)
  jup = irank_overlap(3)
  jdw = irank_overlap(4)
  kup = irank_overlap(5)
  kdw = irank_overlap(6)

  allocate(commbuf_x(Nd,iy_end-iy_sta+1,iz_end-iz_sta+1,Norb,4))
  allocate(commbuf_y(ix_end-ix_sta+1,Nd,iz_end-iz_sta+1,Norb,4))
  allocate(commbuf_z(ix_end-ix_sta+1,iy_end-iy_sta+1,Nd,Norb,4))

  !send from idw to iup

  if(iup/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iz,iy,ix) 
      do iz=1,iz_end-iz_sta+1
      do iy=1,iy_end-iy_sta+1
      do ix=1,Nd
        commbuf_x(ix,iy,iz,iorb,1)=tpsi(ix_end-Nd+ix,iy+iy_sta-1,iz+iz_sta-1,iorb)
      end do
      end do
      end do
    end do
  end if
  ireq(1) = comm_isend(commbuf_x(:,:,:,:,1:1),iup,3,icomm)
  ireq(2) = comm_irecv(commbuf_x(:,:,:,:,2:2),idw,3,icomm)

  !send from iup to idw

  if(idw/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iz,iy,ix)
      do iz=1,iz_end-iz_sta+1
      do iy=1,iy_end-iy_sta+1
      do ix=1,Nd
        commbuf_x(ix,iy,iz,iorb,3)=tpsi(ix_sta+ix-1,iy+iy_sta-1,iz+iz_sta-1,iorb)
      end do
      end do
      end do
    end do
  end if
  ireq(3) = comm_isend(commbuf_x(:,:,:,:,3:3),idw,4,icomm)
  ireq(4) = comm_irecv(commbuf_x(:,:,:,:,4:4),iup,4,icomm)

  !send from jdw to jup

  if(jup/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iz,iy,ix) 
      do iz=1,iz_end-iz_sta+1
      do iy=1,Nd
      do ix=1,ix_end-ix_sta+1
        commbuf_y(ix,iy,iz,iorb,1)=tpsi(ix+ix_sta-1,iy_end-Nd+iy,iz+iz_sta-1,iorb)
      end do
      end do
      end do
    end do
  end if
  ireq(5) = comm_isend(commbuf_y(:,:,:,:,1:1),jup,5,icomm)
  ireq(6) = comm_irecv(commbuf_y(:,:,:,:,2:2),jdw,5,icomm)

  !send from jup to jdw

  if(jdw/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iz,iy,ix)
      do iz=1,iz_end-iz_sta+1
      do iy=1,Nd
      do ix=1,ix_end-ix_sta+1
        commbuf_y(ix,iy,iz,iorb,3)=tpsi(ix+ix_sta-1,iy_sta+iy-1,iz+iz_sta-1,iorb)
      end do
      end do
      end do
    end do
  end if
  ireq(7) = comm_isend(commbuf_y(:,:,:,:,3:3),jdw,6,icomm)
  ireq(8) = comm_irecv(commbuf_y(:,:,:,:,4:4),jup,6,icomm)

  !send from kdw to kup

  if(kup/=comm_proc_null)then
    do iorb=1,Norb
      do iz=1,Nd
  !$OMP parallel do private(iy,ix)
      do iy=1,iy_end-iy_sta+1
      do ix=1,ix_end-ix_sta+1
        commbuf_z(ix,iy,iz,iorb,1)=tpsi(ix+ix_sta-1,iy+iy_sta-1,iz_end-Nd+iz,iorb)
      end do
      end do
      end do
    end do
  end if
  ireq( 9) = comm_isend(commbuf_z(:,:,:,:,1:1),kup,7,icomm)
  ireq(10) = comm_irecv(commbuf_z(:,:,:,:,2:2),kdw,7,icomm)

  !send from kup to kdw

  if(kdw/=comm_proc_null)then
    do iorb=1,Norb
      do iz=1,Nd
  !$OMP parallel do private(iy,ix) 
      do iy=1,iy_end-iy_sta+1
      do ix=1,ix_end-ix_sta+1
        commbuf_z(ix,iy,iz,iorb,3)=tpsi(ix+ix_sta-1,iy+iy_sta-1,iz_sta+iz-1,iorb)
      end do
      end do
      end do
    end do
  end if
  ireq(11) = comm_isend(commbuf_z(:,:,:,:,3:3),kdw,8,icomm)
  ireq(12) = comm_irecv(commbuf_z(:,:,:,:,4:4),kup,8,icomm)


  call comm_wait_all(ireq(1:2))
  if(idw/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iy,ix)
      do iz=1,iz_end-iz_sta+1
      do iy=1,iy_end-iy_sta+1
      do ix=1,Nd
        tpsi(ix_sta-1-Nd+ix,iy+iy_sta-1,iz+iz_sta-1,iorb)=commbuf_x(ix,iy,iz,iorb,2)
      end do
      end do
      end do
    end do
  end if

  call comm_wait_all(ireq(3:4))
  if(iup/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iz,iy,ix)
      do iz=1,iz_end-iz_sta+1
      do iy=1,iy_end-iy_sta+1
      do ix=1,Nd
        tpsi(ix_end+ix,iy+iy_sta-1,iz+iz_sta-1,iorb)=commbuf_x(ix,iy,iz,iorb,4)
      end do
      end do
      end do
    end do
  end if

  call comm_wait_all(ireq(5:6))
  if(jdw/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iz,iy,ix)
      do iz=1,iz_end-iz_sta+1
      do iy=1,Nd
      do ix=1,ix_end-ix_sta+1
        tpsi(ix+ix_sta-1,iy_sta-1-Nd+iy,iz+iz_sta-1,iorb)=commbuf_y(ix,iy,iz,iorb,2)
      end do
      end do
      end do
    end do
  end if

  call comm_wait_all(ireq(7:8))
  if(jup/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iz,iy,ix) 
      do iz=1,iz_end-iz_sta+1
      do iy=1,Nd
      do ix=1,ix_end-ix_sta+1
        tpsi(ix+ix_sta-1,iy_end+iy,iz+iz_sta-1,iorb)=commbuf_y(ix,iy,iz,iorb,4)
      end do
      end do
      end do
    end do
  end if

  call comm_wait_all(ireq(9:10))
  if(kdw/=comm_proc_null)then
    do iorb=1,Norb
      do iz=1,Nd
  !$OMP parallel do private(iy,ix) 
      do iy=1,iy_end-iy_sta+1
      do ix=1,ix_end-ix_sta+1
        tpsi(ix+ix_sta-1,iy+iy_sta-1,iz_sta-1-Nd+iz,iorb)=commbuf_z(ix,iy,iz,iorb,2)
      end do
      end do
      end do
    end do
  end if

  call comm_wait_all(ireq(11:12))
  if(kup/=comm_proc_null)then
    do iorb=1,Norb
      do iz=1,Nd
  !$OMP parallel do private(iy,ix)
      do iy=1,iy_end-iy_sta+1
      do ix=1,ix_end-ix_sta+1
        tpsi(ix+ix_sta-1,iy+iy_sta-1,iz_end+iz,iorb)=commbuf_z(ix,iy,iz,iorb,4)
      end do
      end do
      end do
    end do
  end if

  deallocate(commbuf_x,commbuf_y,commbuf_z)

  return
end subroutine update_overlap_R

!===================================================================================================================================

subroutine update_overlap_C(tpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb,Nd &
                     ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,irank_overlap,icomm)
  use salmon_communication, only: comm_proc_null, comm_isend, comm_irecv, comm_wait_all
  implicit none
  integer,intent(in) :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,Norb,Nd &
                       ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,irank_overlap(6),icomm
  complex(8) :: tpsi(ipx_sta:ipx_end,ipy_sta:ipy_end,ipz_sta:ipz_end,1:Norb)
  !
  integer :: ix,iy,iz,iorb
  integer :: iup,idw,jup,jdw,kup,kdw
  integer :: ireq(12)

  complex(8),allocatable :: commbuf_x(:,:,:,:,:),commbuf_y(:,:,:,:,:),commbuf_z(:,:,:,:,:)

  iup = irank_overlap(1)
  idw = irank_overlap(2)
  jup = irank_overlap(3)
  jdw = irank_overlap(4)
  kup = irank_overlap(5)
  kdw = irank_overlap(6)

  allocate(commbuf_x(Nd,iy_end-iy_sta+1,iz_end-iz_sta+1,Norb,4))
  allocate(commbuf_y(ix_end-ix_sta+1,Nd,iz_end-iz_sta+1,Norb,4))
  allocate(commbuf_z(ix_end-ix_sta+1,iy_end-iy_sta+1,Nd,Norb,4))

  !send from idw to iup

  if(iup/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iz,iy,ix) 
      do iz=1,iz_end-iz_sta+1
      do iy=1,iy_end-iy_sta+1
      do ix=1,Nd
        commbuf_x(ix,iy,iz,iorb,1)=tpsi(ix_end-Nd+ix,iy+iy_sta-1,iz+iz_sta-1,iorb)
      end do
      end do
      end do
    end do
  end if
  ireq(1) = comm_isend(commbuf_x(:,:,:,:,1:1),iup,3,icomm)
  ireq(2) = comm_irecv(commbuf_x(:,:,:,:,2:2),idw,3,icomm)

  !send from iup to idw

  if(idw/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iz,iy,ix)
      do iz=1,iz_end-iz_sta+1
      do iy=1,iy_end-iy_sta+1
      do ix=1,Nd
        commbuf_x(ix,iy,iz,iorb,3)=tpsi(ix_sta+ix-1,iy+iy_sta-1,iz+iz_sta-1,iorb)
      end do
      end do
      end do
    end do
  end if
  ireq(3) = comm_isend(commbuf_x(:,:,:,:,3:3),idw,4,icomm)
  ireq(4) = comm_irecv(commbuf_x(:,:,:,:,4:4),iup,4,icomm)

  !send from jdw to jup

  if(jup/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iz,iy,ix) 
      do iz=1,iz_end-iz_sta+1
      do iy=1,Nd
      do ix=1,ix_end-ix_sta+1
        commbuf_y(ix,iy,iz,iorb,1)=tpsi(ix+ix_sta-1,iy_end-Nd+iy,iz+iz_sta-1,iorb)
      end do
      end do
      end do
    end do
  end if
  ireq(5) = comm_isend(commbuf_y(:,:,:,:,1:1),jup,5,icomm)
  ireq(6) = comm_irecv(commbuf_y(:,:,:,:,2:2),jdw,5,icomm)

  !send from jup to jdw

  if(jdw/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iz,iy,ix)
      do iz=1,iz_end-iz_sta+1
      do iy=1,Nd
      do ix=1,ix_end-ix_sta+1
        commbuf_y(ix,iy,iz,iorb,3)=tpsi(ix+ix_sta-1,iy_sta+iy-1,iz+iz_sta-1,iorb)
      end do
      end do
      end do
    end do
  end if
  ireq(7) = comm_isend(commbuf_y(:,:,:,:,3:3),jdw,6,icomm)
  ireq(8) = comm_irecv(commbuf_y(:,:,:,:,4:4),jup,6,icomm)

  !send from kdw to kup

  if(kup/=comm_proc_null)then
    do iorb=1,Norb
      do iz=1,Nd
  !$OMP parallel do private(iy,ix)
      do iy=1,iy_end-iy_sta+1
      do ix=1,ix_end-ix_sta+1
        commbuf_z(ix,iy,iz,iorb,1)=tpsi(ix+ix_sta-1,iy+iy_sta-1,iz_end-Nd+iz,iorb)
      end do
      end do
      end do
    end do
  end if
  ireq( 9) = comm_isend(commbuf_z(:,:,:,:,1:1),kup,7,icomm)
  ireq(10) = comm_irecv(commbuf_z(:,:,:,:,2:2),kdw,7,icomm)

  !send from kup to kdw

  if(kdw/=comm_proc_null)then
    do iorb=1,Norb
      do iz=1,Nd
  !$OMP parallel do private(iy,ix) 
      do iy=1,iy_end-iy_sta+1
      do ix=1,ix_end-ix_sta+1
        commbuf_z(ix,iy,iz,iorb,3)=tpsi(ix+ix_sta-1,iy+iy_sta-1,iz_sta+iz-1,iorb)
      end do
      end do
      end do
    end do
  end if
  ireq(11) = comm_isend(commbuf_z(:,:,:,:,3:3),kdw,8,icomm)
  ireq(12) = comm_irecv(commbuf_z(:,:,:,:,4:4),kup,8,icomm)


  call comm_wait_all(ireq(1:2))
  if(idw/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iy,ix)
      do iz=1,iz_end-iz_sta+1
      do iy=1,iy_end-iy_sta+1
      do ix=1,Nd
        tpsi(ix_sta-1-Nd+ix,iy+iy_sta-1,iz+iz_sta-1,iorb)=commbuf_x(ix,iy,iz,iorb,2)
      end do
      end do
      end do
    end do
  end if

  call comm_wait_all(ireq(3:4))
  if(iup/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iz,iy,ix)
      do iz=1,iz_end-iz_sta+1
      do iy=1,iy_end-iy_sta+1
      do ix=1,Nd
        tpsi(ix_end+ix,iy+iy_sta-1,iz+iz_sta-1,iorb)=commbuf_x(ix,iy,iz,iorb,4)
      end do
      end do
      end do
    end do
  end if

  call comm_wait_all(ireq(5:6))
  if(jdw/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iz,iy,ix)
      do iz=1,iz_end-iz_sta+1
      do iy=1,Nd
      do ix=1,ix_end-ix_sta+1
        tpsi(ix+ix_sta-1,iy_sta-1-Nd+iy,iz+iz_sta-1,iorb)=commbuf_y(ix,iy,iz,iorb,2)
      end do
      end do
      end do
    end do
  end if

  call comm_wait_all(ireq(7:8))
  if(jup/=comm_proc_null)then
    do iorb=1,Norb
  !$OMP parallel do private(iz,iy,ix) 
      do iz=1,iz_end-iz_sta+1
      do iy=1,Nd
      do ix=1,ix_end-ix_sta+1
        tpsi(ix+ix_sta-1,iy_end+iy,iz+iz_sta-1,iorb)=commbuf_y(ix,iy,iz,iorb,4)
      end do
      end do
      end do
    end do
  end if

  call comm_wait_all(ireq(9:10))
  if(kdw/=comm_proc_null)then
    do iorb=1,Norb
      do iz=1,Nd
  !$OMP parallel do private(iy,ix) 
      do iy=1,iy_end-iy_sta+1
      do ix=1,ix_end-ix_sta+1
        tpsi(ix+ix_sta-1,iy+iy_sta-1,iz_sta-1-Nd+iz,iorb)=commbuf_z(ix,iy,iz,iorb,2)
      end do
      end do
      end do
    end do
  end if

  call comm_wait_all(ireq(11:12))
  if(kup/=comm_proc_null)then
    do iorb=1,Norb
      do iz=1,Nd
  !$OMP parallel do private(iy,ix)
      do iy=1,iy_end-iy_sta+1
      do ix=1,ix_end-ix_sta+1
        tpsi(ix+ix_sta-1,iy+iy_sta-1,iz_end+iz,iorb)=commbuf_z(ix,iy,iz,iorb,4)
      end do
      end do
      end do
    end do
  end if

  deallocate(commbuf_x,commbuf_y,commbuf_z)

  return
end subroutine update_overlap_C

!===================================================================================================================================

end module update_overlap_sub
