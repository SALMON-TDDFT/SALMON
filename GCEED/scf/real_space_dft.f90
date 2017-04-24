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

MODULE global_variables_scf

use scf_data
use hpsi2_sub
use allocate_mat_sub
use deallocate_mat_sub
use new_world_sub
use init_sendrecv_sub
use copy_psi_mesh_sub
use Total_Energy_sub
use calc_density_sub
use change_order_sub
use allocate_sendrecv_groupob_sub
use read_pslfile_sub
use allocate_psl_sub
implicit none

integer :: iDiterYBCG

END MODULE global_variables_scf

!=======================================================================

subroutine Real_Space_DFT(nprocs,nprocid)
!$ use omp_lib
use global_variables_scf
implicit none

integer :: ix,iy,iz,ik,ikoa
integer :: is
integer :: iter,iatom,iob,p1,p2,p5,ii,jj,iflag,ibox
integer :: ibox1,ibox3
integer :: ibox_array(3)
real(8) :: box
real(8) :: sum0,sum1
real(8),allocatable :: rbox_array(:)
real(8),allocatable :: rbox_array2(:)
character(100) :: file_OUT,file_IN,LDA_Info,file_ini
character(100) :: file_atoms_coo
real(8),allocatable :: tpsi_old(:,:,:)
real(8),allocatable :: tpsi(:,:,:)
complex(8),allocatable :: ztpsi_old(:,:,:)
complex(8),allocatable :: ztpsi(:,:,:)
integer :: ibox2
character(100)::fileRho
real(8) :: rNebox1,rNebox2
complex(8),allocatable :: shtpsi(:,:,:,:,:)
complex(8),allocatable :: zpsi_tmp(:,:,:,:,:)

integer :: nprocs,nprocid

fileRho = "density.cube"

nproc=nprocs
myrank=nprocid

iSCFRT=1
ihpsieff=0
iflag_comm_rho=1

iblacsinit=0

elp3(:)=0.d0
elp3(101)=MPI_Wtime()

inumcpu_check=0

call setbN
call setcN

call read_input_scf(file_IN,file_OUT,LDA_Info,file_ini,iDiterYBCG,file_atoms_coo)

if(ilsda==0)then
  call calc_iobnum(itotMST,nproc_ob,newrank_comm_grid,iobnum,nproc_ob,iparaway_ob)
else if(ilsda==1)then
  if(nproc_ob==1)then
    iobnum=itotMST
  else
    if(newrank_comm_spin<nproc_ob_spin(1))then
      call calc_iobnum(MST(1),nproc_ob_spin(1),newrank_comm_grid,iobnum,nproc_ob_spin(1),iparaway_ob)
    else
      call calc_iobnum(MST(2),nproc_ob_spin(2),newrank_comm_grid,iobnum,nproc_ob_spin(2),iparaway_ob)
    end if
  end if
end if

Structure_Optimization_Iteration : do istopt=1,iter_stopt
Multigrid_Iteration : do img=1,ntmg

elp3(102)=MPI_Wtime()

if(istopt==1)then
  select case( IC )

!------------------------------ New calculation

  case default

    Hgs(1:3)=Harray(1:3,1)
    Hvol=Hgs(1)*Hgs(2)*Hgs(3)
    Miter = 0        ! Miter: Iteration counter set to zero
    call init_mesh(img)
    call set_gridcoo
    call init_mesh_s

    call init_updown
    call init_itype
    call init_sendrecv_matrix
    if(MEO==2.or.MEO==3) call make_corr_pole
    call make_icoobox_bound
        
    call allocate_mat
    call allocate_sendrecv_groupob
    allocate( Vpsl(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
    if(icalcforce==1)then
      allocate( Vpsl_atom(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),MI) )
    end if

    if(iflag_ps.eq.0)then
      Vpsl=0d0
    else
      call read_pslfile
      call allocate_psl
      call init_ps
    end if

    if(iobnum >= 1)then
      allocate( psi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum,1) ) 
    end if
    if(iswitch_orbital_mesh==1.or.iflag_subspace_diag==1)then
      allocate( psi_mesh(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:itotMST,1) ) 
    end if

    call init_wf_ns(1)

    call Gram_Schmidt_ns

    allocate( rho(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )  
    allocate( rho_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:num_rho_stock+1) )  
    allocate( rho_out(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:num_rho_stock) ) 
    rho_in=0.d0
    rho_out=0.d0
                                
    if(ilsda == 1)then
      allocate( rho_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2) )  
      allocate( rho_s_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2,1:num_rho_stock+1) )  
      allocate( rho_s_out(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2,1:num_rho_stock) )  
      rho_s_in=0.d0
      rho_s_out=0.d0
    end if
    rho=0.d0 

    call calc_density(psi,1)

    if(ilsda==0)then
      allocate (Vlocal(mg_sta(1):mg_end(1),  &
              mg_sta(2):mg_end(2),  &
              mg_sta(3):mg_end(3),1))
    else if(ilsda==1)then
      allocate (Vlocal(mg_sta(1):mg_end(1),  &
              mg_sta(2):mg_end(2),  &
              mg_sta(3):mg_end(3),2))
    end if

    allocate( Vh(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )  
    Vh=0.d0

    call Hartree_ns

    
    if(ilsda == 0) then
      allocate( Vxc(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )  
    else if(ilsda == 1) then
      allocate( Vxc_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2) )  
    end if
    allocate( esp(itotMST,1) )

    if(ilsda==0)then
      call conv_core_exc_cor
    else if(ilsda==1)then
      call Exc_cor_ns
    end if

    call mpi_allgatherv_vlocal

    call Total_Energy(psi)
      
!------------------------------ Continue the previous calculation

  case(1,3)
    call IN_data(file_IN)

    call allocate_mat
    call allocate_sendrecv_groupob

    if(iflag_ps/=0) then
      call read_pslfile
      call allocate_psl
      call init_ps
    end if



    call init_updown
    call init_itype
    call init_sendrecv_matrix
    if(MEO==2.or.MEO==3) call make_corr_pole
    call make_icoobox_bound
  end select

else if(istopt>=2)then
  Miter = 0        ! Miter: Iteration counter set to zero
  if(iflag_ps/=0) then
    call init_ps
  end if
end if

elp3(103)=MPI_Wtime()

if(myrank.eq.0) then
  write(*,*) '-----------------------------------------------'
  if(iflag_diisjump == 0) then
    write(*,'(1x,"iter =",i6,5x,"Total Energy =",f19.8,5x,"Vh iteration =",i4)') 0,Etot*2d0*Ry,iterVh
  else if(iflag_diisjump == 1) then
    write(*,'(1x,"iter =",i6,5x,"Total Energy =",f19.8,5x,"Vh iteration =",i4,"**")') 0,Etot*2d0*Ry,iterVh
  end if
  do p5=1,(itotMST+3)/4
    p1=4*(p5-1)+1
    p2=4*p5 ; if ( p2 > itotMST ) p2=itotMST
    write(*,'(1x,4(i5,f15.4,2x))') (iob,esp(iob,1)*2d0*Ry,iob=p1,p2)
  end do
end if
!---------------------------------------- Iteration

iflag=1
iterVh=1000
sum1=1.0d9

iflag_diisjump=0

allocate(idiis_sd(itotMST))
idiis_sd=0

if(img==1.and.istopt==1) allocate(norm_diff_psi_stock(itotMST,1))
norm_diff_psi_stock=1.0d9

if(img>=2.or.istopt>=2) deallocate(rho_stock,Vlocal_stock)
allocate(rho_stock(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1))
if(ilsda==0)then
  allocate(Vlocal_stock(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1))
else
  allocate(Vlocal_stock(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:2))
end if

if(ilsda==0)then
!$OMP parallel do
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    rho_stock(ix,iy,iz,1)=rho(ix,iy,iz)
    Vlocal_stock(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)
  end do
  end do
  end do
else
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    rho_stock(ix,iy,iz,1)=rho(ix,iy,iz)
    Vlocal_stock(ix,iy,iz,1:2)=Vlocal(ix,iy,iz,1:2)
  end do
  end do
  end do
end if

DFT_Iteration : do iter=1,iDiter(img)

  elp3(111)=MPI_Wtime()

  if(iflag_convergence==1)then
    if(iterVh <= ithresholdVh(img)) cycle DFT_Iteration
  else if(iflag_convergence==2)then
    if(sum1<threshold_norm_diff_rho(img)) cycle DFT_Iteration
  else if(iflag_convergence==3)then
    if(sum1<threshold_square_norm_diff_Vlocal(img)) cycle DFT_Iteration
  end if

  elp3(112)=MPI_Wtime()
  elp3(122)=elp3(122)+elp3(112)-elp3(111)

  Miter=Miter+1

  call calc_occupation(iter)

  call copy_density

  if(iscf_order==1)then
   
    if( minroutine == 1 .or.       &
   (minroutine == 4 .and. Miter <= iDiterYBCG) ) then
      elp3(181)=MPI_Wtime()
      call DTcg(psi,iflag)
      elp3(182)=MPI_Wtime()
      elp3(183)=elp3(183)+elp3(182)-elp3(181)
    else if( minroutine == 3 .or. minroutine == 4 ) then
      elp3(181)=MPI_Wtime()
      call rmmdiis(psi,iflag)
      elp3(182)=MPI_Wtime()
      elp3(184)=elp3(184)+elp3(182)-elp3(181)
    end if
  
  
    elp3(113)=MPI_Wtime()
    elp3(123)=elp3(123)+elp3(113)-elp3(112)
  
    call Gram_Schmidt_ns
  
    if(iflag_subspace_diag==1)then
      if(Miter>iDiter_nosubspace_diag)then
        call subspace_diag
      end if
    end if
  
    elp3(114)=MPI_Wtime()
    elp3(124)=elp3(124)+elp3(114)-elp3(113)
     
    call calc_density(psi,2)
  
    call simple_mixing(1.d0-mixrate,mixrate)
    
    call calc_rho_in
  
    elp3(115)=MPI_Wtime()
    elp3(125)=elp3(125)+elp3(115)-elp3(114)
  
    if(imesh_s_all==1.or.(imesh_s_all==0.and.myrank<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
      call Hartree_ns
    end if
  
    elp3(116)=MPI_Wtime()
    elp3(126)=elp3(126)+elp3(116)-elp3(115)
  
    if(imesh_s_all==1.or.(imesh_s_all==0.and.myrank<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
      if(ilsda==0)then
        call conv_core_exc_cor
      else if(ilsda==1)then
        call Exc_cor_ns
      end if
    end if
   
    call mpi_allgatherv_vlocal
    
    elp3(117)=MPI_Wtime()
    elp3(127)=elp3(127)+elp3(117)-elp3(116)
  
    call Total_Energy(psi)
  
    elp3(118)=MPI_Wtime()
    elp3(128)=elp3(128)+elp3(118)-elp3(117)
    elp3(131)=MPI_Wtime()
  
    elp3(132)=MPI_Wtime()
    elp3(142)=elp3(142)+elp3(132)-elp3(131)
    
    elp3(118)=MPI_Wtime()
  
    call change_order(psi)
  
  else if(iscf_order==2)then

    call Gram_Schmidt_ns

    if(Miter>iDiter_nosubspace_diag)then
      call subspace_diag
    end if

    call Gram_Schmidt_ns

    if( minroutine == 1 .or. (minroutine == 4 .and. Miter <= iDiterYBCG) ) then
      call DTcg(psi,iflag)
    else if( minroutine == 3 .or. minroutine == 4 ) then
      call rmmdiis(psi,iflag)
    end if

    call Gram_Schmidt_ns

    call calc_density(psi,2)
  
    call simple_mixing(1.d0-mixrate,mixrate)
    
    call calc_rho_in

    if(imesh_s_all==1.or.(imesh_s_all==0.and.myrank<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
      call Hartree_ns
    end if
  
    if(imesh_s_all==1.or.(imesh_s_all==0.and.myrank<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
      if(ilsda==0)then
        call conv_core_exc_cor
      else if(ilsda==1)then
        call Exc_cor_ns
      end if
    end if
   
    call mpi_allgatherv_vlocal
    
    call Gram_Schmidt_ns

    call Total_Energy(psi)
  end if

  if(iflag_convergence==2)then
    sum0=0.d0
!$OMP parallel do reduction(+:sum0)
    do iz=ng_sta(3),ng_end(3) 
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      sum0=sum0+(rho(ix,iy,iz)-rho_stock(ix,iy,iz,1))**2
    end do
    end do
    end do
    call MPI_Allreduce(sum0,sum1,1,MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_h,ierr)
    sum1=sum1*Hvol
  else if(iflag_convergence==3)then
    sum0=0.d0
!$OMP parallel do reduction(+:sum0)
    do iz=ng_sta(3),ng_end(3) 
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      sum0=sum0+(Vlocal(ix,iy,iz,1)-Vlocal_stock(ix,iy,iz,1))**2
    end do
    end do
    end do
    call MPI_Allreduce(sum0,sum1,1,MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_h,ierr)
    sum1=sum1*Hvol
  end if

  if(myrank.eq.0) then
    write(*,*) '-----------------------------------------------'
    if(iflag_diisjump == 1) then
      write(*,'("Diisjump occured. Steepest descent was used.")')
    end if
    write(*,'(1x,"iter =",i6,5x,"Total Energy =",f19.8,5x,"Vh iteration =",i4)') Miter,Etot*2d0*Ry,iterVh
    do p5=1,(itotMST+3)/4
      p1=4*(p5-1)+1
      p2=4*p5 ; if ( p2 > itotMST ) p2=itotMST
      write(*,'(1x,4(i5,f15.4,2x))') (iob,esp(iob,1)*2d0*Ry,iob=p1,p2)
    end do
    if(iflag_convergence==2)then
      write(*,'("iter and square of norm of difference of rho = ",i6,e15.8)') Miter,sum1/a_B**3
    else if(iflag_convergence==3)then
      write(*,'("iter and square of norm of difference of Vlocal = ",i6,e15.8)') Miter, sum1*(2.d0*Ry)**2*a_B**3
    end if
  end if 
  rNebox1=0.d0 
!$OMP parallel do reduction(+:rNebox1)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    rNebox1=rNebox1+rho(ix,iy,iz)
  end do
  end do
  end do
  call MPI_Allreduce(rNebox1,rNebox2,1,  &
&             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(myrank==0)then
    write(*,*) "Ne=",rNebox2*Hvol
  end if

  elp3(119)=MPI_Wtime()
  elp3(129)=elp3(129)+elp3(119)-elp3(118)
  elp3(130)=elp3(130)+elp3(119)-elp3(111)

if(ilsda==0)then
!$OMP parallel do
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    rho_stock(ix,iy,iz,1)=rho(ix,iy,iz)
    Vlocal_stock(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)
  end do
  end do
  end do
else if(ilsda==1)then
!$OMP parallel do
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    rho_stock(ix,iy,iz,1)=rho(ix,iy,iz)
    Vlocal_stock(ix,iy,iz,1:2)=Vlocal(ix,iy,iz,1:2)
  end do
  end do
  end do
end if

end do DFT_Iteration
elp3(104)=MPI_Wtime()

deallocate(idiis_sd)

if(icalcforce==1) call calc_force

if(iflag_stopt==1) then
  call structure_opt
  if(myrank==0)then
    write(*,*) "atomic coordinate"
    do iatom=1,MI
      write(*,'(a3,3f16.8,2i3)') AtomName(iatom), (Rion(jj,iatom)*a_B,jj=1,3),Kion(iatom),istopt_a(iatom)
    end do
  end if
end if

end do Multigrid_Iteration

end do Structure_Optimization_Iteration

!---------------------------------------- Output

if(iflag_writepsi==1) then
  call writepsi
end if

if(iflag_dos==1) then
  call calc_dos
end if

if(iflag_pdos==1) then
  call calc_pdos
end if

if(OC==2)then
  call prep_ini(file_ini)
end if

if(iflag_ELF==1)then
  allocate(elf(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),      &
               lg_sta(3):lg_end(3)))
  iSCFRT=1
  call calcELF
  deallocate(elf)
end if

! LDA data
! subroutines in scf_data.f90
if ( OC==1.or.OC==2.or.OC==3 ) then
  call OUT_data(file_OUT) !
  call outRho(fileRho)
end if
elp3(105)=MPI_Wtime()

! LDA information

if(myrank.eq.0) then
  open(1,file=LDA_info)

  write(1,*) "Total number of iteration = ", Miter
  write(1,*)
  if(ilsda == 0) then
     write(1,*) "Number of orbitals = ", MST(1)
     write(1,*) "Number of filled orbitals = ", ifMST(1)
  else if(ilsda == 1) then
     write(1,*) "Number of orbitals = ", (MST(is),is=1,2)
     write(1,*) "Number of filled orbitals = ", (ifMST(is),is=1,2)
  end if
  write(1,*)
  write(1,*) "Total energy = ", Etot*2d0*Ry
  write(1,*) "1-particle energies"
  do p5=1,(itotMST+3)/4
    p1=4*(p5-1)+1
    p2=4*p5 ; if ( p2 > itotMST ) p2=itotMST
    write(1,'(1x,4(i5,f15.4,2x))') (iob,esp(iob,1)*2d0*Ry,iob=p1,p2)
  end do
  write(1,*)

  do ii=1,ntmg
     write(1,*) "Size of the box     = ", rLsize(:,ii)*a_B
  end do

  write(1,*) "Grid spacing          = ", (Hgs(jj)*a_B,jj=1,3)
  write(1,*)
  write(1,*) "Number of atoms = ", MI
  do ik=1,MKI
    write(1,'(1x,"iZatom(",i2,")       = ",i12)') ik, iZatom(ik)
  end do
  write(1,*)
  write(1,*) "Ref. and max angular momentum",      &
        " and pseudo-core radius of PP"
  do ikoa=1,MKI
     write(1,'(1x,"(",i2,")  "," Ref, Max, Rps =",2i4,f8.3)')      &
                              ikoa,Lref(ikoa),Mlps(ikoa),Rps(ikoa)*a_B
  end do

  close(1)

end if
elp3(106)=MPI_Wtime()

if(myrank==0)then
   write(*,'(a)') "==================== elapsed time ===================="
   if(IC==0)then
     write(*,'(a,f16.8)') "elapsed time before initializing [s]  = ", elp3(102)-elp3(101)
     write(*,'(a,f16.8)') "elapsed time for initializing [s]     = ", elp3(103)-elp3(102)
   else if(IC==1)then
     write(*,'(a,f16.8)') "elapsed time before reading data [s]  = ", elp3(102)-elp3(101)
     write(*,'(a,f16.8)') "elapsed time for reading data [s]     = ", elp3(103)-elp3(102)
   end if
   write(*,'(a,f16.8)') "elapsed time for scf iterations [s]   = ", elp3(104)-elp3(103)
   write(*,'(a,f16.8)') "elapsed time for writing data [s]     = ", elp3(105)-elp3(104)
   write(*,'(a,f16.8)') "elapsed time after writing data [s]   = ", elp3(106)-elp3(105)
   write(*,'(a,f16.8)') "total time [s]                        = ", elp3(106)-elp3(101)
   write(*,'(a)') "======================================================"
   write(*,'(a)') "================== in scf iterations ================="
   write(*,'(a,f16.8)') "elapsed time for copying psi [s]      = ", elp3(122)
   write(*,'(a,f16.8)') "elapsed time for optimizing psi [s]   = ", elp3(123)
   write(*,'(a,f16.8)') "elapsed time for Gram Schmidt [s]     = ", elp3(124)
   write(*,'(a,f16.8)') "elapsed time for calculating rho  [s] = ", elp3(125)
   write(*,'(a,f16.8)') "elapsed time for Hartree routine  [s] = ", elp3(126)
   write(*,'(a,f16.8)') "elapsed time for Exc_Cor routine  [s] = ", elp3(127)
   write(*,'(a,f16.8)') "elapsed time for calculating Etot [s] = ", elp3(128)
   write(*,'(a,f16.8)') "elapsed time for subspace-diag. [s]   = ", elp3(142)
   write(*,'(a,f16.8)') "elapsed time for writing info. [s]    = ", elp3(129)
   write(*,'(a,f16.8)') "total time for scf iterations [s]     = ", elp3(130)
   write(*,'(a)') "======================================================"
   write(*,'(a)') "================== in subspace-diag. ================="
   write(*,'(a,f16.8)') "elapsed time for initialization [s]   = ", elp3(352)
   write(*,'(a,f16.8)') "elapsed time for Vlocal [s]           = ", elp3(353)
   write(*,'(a,f16.8)') "elapsed time for Amat [s]             = ", elp3(354)
   write(*,'(a,f16.8)') "elapsed time for eigen [s]            = ", elp3(355)
   write(*,'(a,f16.8)') "elapsed time for psi_mesh_box (1) [s] = ", elp3(356)
   write(*,'(a,f16.8)') "elapsed time for psi_mesh_box (2) [s] = ", elp3(357)
   write(*,'(a)') "======================================================"
   write(*,'(a)') "================== in eigen. ========================="
   write(*,'(a,f16.8)') "elapsed time for initialization [s]   = ", elp3(362)
   write(*,'(a,f16.8)') "elapsed time for BLACS, DESCINIT [s]  = ", elp3(363)
   write(*,'(a,f16.8)') "elapsed time for PDELSET (1) [s]      = ", elp3(364)
   write(*,'(a,f16.8)') "elapsed time for PDSYEVX [s]          = ", elp3(365)
   write(*,'(a,f16.8)') "elapsed time for PDELSET (2) [s]      = ", elp3(366)
   write(*,'(a)') "======================================================"
   write(*,'(a,f16.8)') "elapsed time for CG [s]               = ", elp3(183)
   write(*,'(a,f16.8)') "elapsed time for DIIS [s]             = ", elp3(184)
   write(*,'(a,f16.8)') "comm. for inner product in CG, DIIS[s]= ", elp3(190)
   write(*,'(a,f16.8)') "elapsed time for Gram Schmidt [s]     = ", elp3(124)
   write(*,'(a,f16.8)') "comm. for inner product in GS (1) [s] = ", elp3(188)
   write(*,'(a,f16.8)') "comm. for inner product in GS (2) [s] = ", elp3(189)
   write(*,'(a,f16.8)') "hpsi2 in CG [s]                       = ", elp3(193)
   write(*,'(a,f16.8)') "gk in CG [s]                          = ", elp3(194)
   write(*,'(a)') "======================================================"
end if

deallocate(Vlocal)

END subroutine Real_Space_DFT

!=======================================================================
!========================================= Grid generation and labelling

SUBROUTINE init_mesh
use global_variables_scf
implicit none

real(8) :: rLsize1(3)

if(myrank.eq.0)      &
    print *,"----------------------------------- init_mesh"

rLsize1(:)=rLsize(:,img)
call setlg(lg_sta,lg_end,lg_num,ista_Mx_ori,iend_Mx_ori,inum_Mx_ori,    &
           Hgs,Nd,rLsize1,imesh_oddeven)

allocate(ista_Mxin(3,0:nproc-1),iend_Mxin(3,0:nproc-1))
allocate(inum_Mxin(3,0:nproc-1))

call setmg(mg_sta,mg_end,mg_num,ista_Mxin,iend_Mxin,inum_Mxin,  &
           lg_sta,lg_end,lg_num,nproc,myrank,nproc_Mxin,nproc_ob,isequential)

if(myrank.eq.0) write(*,*) "Mx     =", iend_Mx_ori

return

END SUBROUTINE init_mesh

