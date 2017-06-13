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
subroutine calc_pdos
use salmon_parallel, only: nproc_id_global, nproc_group_grid, nproc_group_orbital
use salmon_communication, only: comm_is_root
use mpi, only: mpi_double_precision, mpi_sum
use inputoutput
use scf_data
use allocate_psl_sub
use new_world_sub
implicit none
integer :: iob,iobmax,iob_allob,iatom,L,ix,iy,iz,iene
integer :: ikoa
integer :: intr
real(8) :: phi_r
real(8) :: rr
real(8) :: ratio1,ratio2
real(8) :: xx,yy,zz
real(8) :: Ylm
integer :: lm
real(8) :: rbox_pdos(25,MI)
real(8) :: rbox_pdos2(25,MI)
real(8) :: rbox_pdos3(-300:300,0:4,MI)
real(8) :: pdos(-300:300,0:4,MI)
real(8),parameter :: sigma_gd=0.01d0
character(100) :: Outfile
integer :: ierr

call calc_pmax(iobmax)

rbox_pdos3=0.d0

do iob=1,iobmax
  call calc_allob(iob,iob_allob)
  rbox_pdos=0.d0
  do iatom=1,MI
    ikoa=Kion(iatom)
    do L=0,Mlps(ikoa)
      do lm=L**2+1,(L+1)**2
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          xx=gridcoo(ix,1)-Rion(1,iatom)
          yy=gridcoo(iy,2)-Rion(2,iatom)
          zz=gridcoo(iz,3)-Rion(3,iatom)
          rr=sqrt(xx**2+yy**2+zz**2)+1.d-50
          call bisection(rr,intr,ikoa)
          ratio1=(rr-rad_psl(intr,ikoa))/(rad_psl(intr+1,ikoa)-rad_psl(intr,ikoa)) ; ratio2=1.d0-ratio1
          phi_r= ratio1*uppr(intr+1,Lref(ikoa),ikoa)+ratio2*uppr(intr,Lref(ikoa),ikoa)
          call Ylm_sub(xx,yy,zz,lm,Ylm)
          rbox_pdos(lm,iatom)=rbox_pdos(lm,iatom)+psi(ix,iy,iz,iob,1)*phi_r*Ylm*Hvol
        end do
        end do
        end do
      end do
    end do
  end do
  call MPI_Allreduce(rbox_pdos,rbox_pdos2,25*MI,MPI_DOUBLE_PRECISION,MPI_SUM,nproc_group_orbital,ierr) 
  do iatom=1,MI
    ikoa=Kion(iatom)
    do L=0,Mlps(ikoa)
      do lm=L**2+1,(L+1)**2
        do iene=-300,300
          rbox_pdos3(iene,L,iatom)=rbox_pdos3(iene,L,iatom)  &
            +abs(rbox_pdos2(lm,iatom))**2*  &
             exp(-(dble(iene)/10d0/au_energy_ev-esp(iob_allob,1))**2/(2.d0*sigma_gd**2))/sqrt(2.d0*Pi*sigma_gd**2)
        end do
      end do
    end do
  end do
end do
call MPI_Allreduce(rbox_pdos3,pdos,601*5*MI,MPI_DOUBLE_PRECISION,MPI_SUM,nproc_group_grid,ierr) 

if(comm_is_root(nproc_id_global))then
  do iatom=1,MI
    ikoa=Kion(iatom)
    write(fileNumber, '(i8)') iatom
    OutFile = "pdos"//trim(adjustl(fileNumber))//".data"
    open(101,file=OutFile)
    write(101,'("# Projected Density of States")') 
    select case(unit_energy)
    case('au','a.u.')
      if(Mlps(ikoa)==0)then
        write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.]")') 
      else if(Mlps(ikoa)==1)then
        write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.]")') 
      else if(Mlps(ikoa)==2)then
        write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.] PDOS(l=2)[a.u.]")') 
      else if(Mlps(ikoa)==3)then
        write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.] PDOS(l=2)[a.u.] PDOS(l=3)[a.u.]")')
      end if 
    case('ev','eV')
      if(Mlps(ikoa)==0)then
        write(101,'("# Energy[eV]  PDOS(l=0)[1/eV]")')
      else if(Mlps(ikoa)==1)then
        write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV]")')
      else if(Mlps(ikoa)==2)then
        write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV] PDOS(l=2)[1/eV]")') 
      else if(Mlps(ikoa)==3)then
        write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV] PDOS(l=2)[1/eV] PDOS(l=3)[1/eV]")')
      end if 
    end select
    write(101,'("#-----------------------")') 
    if(Mlps(ikoa)==0)then
      do iene=-300,300
        write(101,'(f10.5,f14.8)') dble(iene)/10.d0/au_energy_ev*uenergy_from_au,   &
                                   (pdos(iene,L,iatom)*au_energy_ev/uenergy_from_au,L=0,Mlps(ikoa))
      end do
    else if(Mlps(ikoa)==1)then
      do iene=-300,300
        write(101,'(f10.5,2f14.8)') dble(iene)/10.d0/au_energy_ev*uenergy_from_au,   &
                                    (pdos(iene,L,iatom)*au_energy_ev/uenergy_from_au,L=0,Mlps(ikoa))
      end do
    else if(Mlps(ikoa)==2)then
      do iene=-300,300
        write(101,'(f10.5,3f14.8)') dble(iene)/10.d0/au_energy_ev*uenergy_from_au,   &
                                    (pdos(iene,L,iatom)*au_energy_ev/uenergy_from_au,L=0,Mlps(ikoa))
      end do
    else if(Mlps(ikoa)==3)then
      do iene=-300,300
        write(101,'(f10.5,4f14.8)') dble(iene)/10.d0/au_energy_ev*uenergy_from_au,   &
                                    (pdos(iene,L,iatom)*au_energy_ev/uenergy_from_au,L=0,Mlps(ikoa))
      end do
    end if
    close(101)
  end do
end if

end subroutine calc_pdos
