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
!=======================================================================

MODULE global_variables_scf

use inputoutput
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
use read_pslfile_sub
use allocate_psl_sub
use persistent_comm
use structure_opt_sub
implicit none

END MODULE global_variables_scf

!=======================================================================

subroutine Real_Space_DFT
use salmon_parallel, only: nproc_id_global, nproc_size_global, nproc_group_global, &
                           nproc_group_h, nproc_id_kgrid, nproc_id_spin, nproc_id_orbitalgrid, &
                           nproc_group_spin
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_xc, only: init_xc, finalize_xc
use misc_routines, only: get_wtime
use global_variables_scf
implicit none

integer :: ix,iy,iz,ik,ikoa
integer :: is
integer :: iter,iatom,iob,p1,p2,p5,ii,jj,iflag
real(8) :: sum0,sum1
character(100) :: file_atoms_coo
complex(8),allocatable :: zpsi_tmp(:,:,:,:,:)
real(8) :: rNebox1,rNebox2
integer :: itmg

call init_xc(xc_func, ispin, cval, xcname=xc, xname=xname, cname=cname)

iSCFRT=1
ihpsieff=0
iflag_comm_rho=1

iblacsinit=0

elp3(:)=0.d0
elp3(101)=get_wtime()

inumcpu_check=0

call setbN
call setcN

call check_dos_pdos

call convert_input_scf(file_atoms_coo)

call set_filename

call setk(k_sta, k_end, k_num, num_kpoints_rd, nproc_k, nproc_id_orbitalgrid)

if(ilsda==0)then
  call calc_iobnum(itotMST,nproc_ob,nproc_id_kgrid,iobnum,nproc_ob,iparaway_ob)
else if(ilsda==1)then
  if(nproc_ob==1)then
    iobnum=itotMST
  else
    if(nproc_id_spin<nproc_ob_spin(1))then
      call calc_iobnum(MST(1),nproc_ob_spin(1),nproc_id_kgrid,iobnum,nproc_ob_spin(1),iparaway_ob)
    else
      call calc_iobnum(MST(2),nproc_ob_spin(2),nproc_id_kgrid,iobnum,nproc_ob_spin(2),iparaway_ob)
    end if
  end if
end if

if(iflag_stopt==1)then
  call structure_opt_ini(MI)
  istopt_tranc=0
end if

Structure_Optimization_Iteration : do istopt=1,iter_stopt
Multigrid_Iteration : do img=1,ntmg

elp3(102)=get_wtime()

if(istopt==1)then
  select case( IC )

!------------------------------ New calculation

  case default

    Hgs(1:3)=Harray(1:3,1)
    Hvol=Hgs(1)*Hgs(2)*Hgs(3)
    Miter = 0        ! Miter: Iteration counter set to zero
    itmg=img
    call set_imesh_oddeven(itmg)
    call init_mesh(itmg)
    call set_gridcoo
    call init_mesh_s
    call check_ng

    call init_updown
    call init_itype
    call init_sendrecv_matrix
    select case(iperiodic)
    case(0)
      if(MEO==2.or.MEO==3) call make_corr_pole
    end select
    call make_icoobox_bound
        
    call allocate_mat
    call set_icoo1d
    call allocate_sendrecv
    call init_persistent_requests

    if(iperiodic==3)then
      allocate (zpsi_tmp(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
                 1:iobnum,k_sta:k_end))
      allocate(k_rd(3,num_kpoints_rd),ksquare(num_kpoints_rd))
      call init_k_rd(k_rd,ksquare,1)
    end if

    allocate( Vpsl(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
    if(icalcforce==1)then
      allocate( Vpsl_atom(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),MI) )
    end if
    
    if(iperiodic==3)then
      call prep_poisson_fft
    end if

    if(iflag_ps.eq.0)then
      Vpsl=0d0
    else
      call read_pslfile
      call allocate_psl
      call init_ps
    end if

    if(iobnum >= 1)then
      select case(iperiodic)
      case(0)
        allocate( psi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum,k_sta:k_end) )
      case(3)
        allocate( ttpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
        allocate( zpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum,k_sta:k_end) )
      end select
    end if
    if(iswitch_orbital_mesh==1.or.iflag_subspace_diag==1)then
      select case(iperiodic)
      case(0)
        allocate( psi_mesh(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:itotMST,1) ) 
      case(3)
        allocate( zpsi_mesh(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:itotMST,num_kpoints_rd) )
      end select
    end if

    call init_wf_ns(1)

    select case(iperiodic)
    case(0)
      call Gram_Schmidt_ns
    case(3)
      call Gram_Schmidt_periodic
    end select

    allocate( rho(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )  
    allocate( rho_in(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:num_rho_stock+1) )  
    allocate( rho_out(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:num_rho_stock+1) ) 
    rho_in=0.d0
    rho_out=0.d0
                                
    if(ilsda == 1)then
      allocate( rho_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2) )  
      allocate( rho_s_in(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:num_rho_stock+1,2) )  
      allocate( rho_s_out(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:num_rho_stock+1,2) )  
      rho_s_in=0.d0
      rho_s_out=0.d0
    end if
    rho=0.d0 

    select case(iperiodic)
    case(0)
      call calc_density(psi)
    case(3)
      call calc_density(zpsi)
    end select

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
    allocate( esp(itotMST,num_kpoints_rd) )

    call exc_cor_ns

    call allgatherv_vlocal

    select case(iperiodic)
    case(0)
      call Total_Energy(psi)
    case(3)
      call Total_Energy_periodic_scf_esp(zpsi)
      do ik=k_sta,k_end
      do iob=1,iobnum
!$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          zpsi_tmp(ix,iy,iz,iob,ik)=zpsi(ix,iy,iz,iob,ik)
        end do
        end do
        end do
      end do
      end do
      call Total_Energy_periodic_scf(zpsi_tmp)
    end select
      
!------------------------------ Continue the previous calculation

  case(1,3)
    call IN_data

    call allocate_mat
    call set_icoo1d
    call allocate_sendrecv
    call init_persistent_requests

    if(iperiodic==3)then
      allocate (zpsi_tmp(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
                 1:iobnum,k_sta:k_end))
      allocate(k_rd(3,num_kpoints_rd),ksquare(num_kpoints_rd))
      call init_k_rd(k_rd,ksquare,1)
    end if

    if(iperiodic==3)then
      call prep_poisson_fft
    end if

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

elp3(103)=get_wtime()

if(comm_is_root(nproc_id_global)) then
  write(*,*) '-----------------------------------------------'
  select case(iperiodic)
  case(0)
    write(*,'(1x,"iter =",i6,5x,"Total Energy =",f19.8,5x,"Vh iteration =",i4)') Miter,Etot*2d0*Ry,iterVh
  case(3)
    write(*,'(1x,"iter =",i6,5x,"Total Energy =",f19.8)') Miter,Etot*2d0*Ry
  end select
  do ik=1,num_kpoints_rd
    if(ik<=3)then
      if(iperiodic==3) write(*,*) "k=",ik
      do p5=1,(itotMST+3)/4
        p1=4*(p5-1)+1
        p2=4*p5 ; if ( p2 > itotMST ) p2=itotMST
        write(*,'(1x,4(i5,f15.4,2x))') (iob,esp(iob,ik)*2d0*Ry,iob=p1,p2)
      end do
      if(iperiodic==3) write(*,*)
    end if
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
!$OMP parallel do private(iz,iy,ix)
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

  elp3(111)=get_wtime()

  select case(convergence)
    case('rho_dne')
      if(sum1<threshold) cycle DFT_Iteration
    case('norm_rho','norm_rho_dng')
      if(sum1<threshold_norm_rho) cycle DFT_Iteration
    case('norm_pot','norm_pot_dng')
      if(sum1<threshold_norm_pot) cycle DFT_Iteration
  end select 

  elp3(112)=get_wtime()
  elp3(122)=elp3(122)+elp3(112)-elp3(111)

  Miter=Miter+1

  if(temperature_k>=0.d0.and.Miter>iditer_notemperature) then
    if(iperiodic.eq.3) then
      call ne2mu_p
    else
      call ne2mu
    endif
  else
    call calc_occupation
  endif

  call copy_density

  if(iscf_order==1)then
   
    if( amin_routine == 'cg' .or.       &
   (amin_routine == 'cg-diis' .and. Miter <= iDiterYBCG) ) then
      elp3(181)=get_wtime()
      select case(iperiodic)
      case(0)
        select case(gscg)
        case('y')
          call sgscg(psi,iflag)
        case('n')
          call DTcg(psi,iflag)
        end select
      case(3)
        select case(gscg)
        case('y')
          call gscg_periodic(zpsi,iflag)
        case('n')
          call DTcg_periodic(zpsi,iflag)
        end select
      end select
      elp3(182)=get_wtime()
      elp3(183)=elp3(183)+elp3(182)-elp3(181)
    else if( amin_routine  == 'diis' .or. amin_routine == 'cg-diis' ) then
      elp3(181)=get_wtime()
      select case(iperiodic)
      case(0)
        call rmmdiis(psi)
      case(3)
        stop "rmmdiis method is not implemented for periodic systems."
      end select
      elp3(182)=get_wtime()
      elp3(184)=elp3(184)+elp3(182)-elp3(181)
    end if
  
    elp3(113)=get_wtime()
    elp3(123)=elp3(123)+elp3(113)-elp3(112)

    select case(iperiodic)
    case(0)
      call Gram_Schmidt_ns
    case(3)
      call Gram_Schmidt_periodic
    end select
  
    if(iflag_subspace_diag==1)then
      if(Miter>iDiter_nosubspace_diag)then
        select case(iperiodic)
        case(0)
          call subspace_diag
        case(3)
          call subspace_diag_periodic
        end select
      end if
    end if
  
    elp3(114)=get_wtime()
    elp3(124)=elp3(124)+elp3(114)-elp3(113)
    

    select case(iperiodic)
    case(0)
      call calc_density(psi)
    case(3)
      call calc_density(zpsi)
    end select

    select case(amixing)
      case ('simple')
        call simple_mixing(1.d0-rmixrate,rmixrate)
      case ('broyden')
        call buffer_broyden_ns(iter)
    end select
    
    elp3(115)=get_wtime()
    elp3(125)=elp3(125)+elp3(115)-elp3(114)
  
    if(imesh_s_all==1.or.(imesh_s_all==0.and.nproc_id_global<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
      call Hartree_ns
    end if
  
    elp3(116)=get_wtime()
    elp3(126)=elp3(126)+elp3(116)-elp3(115)
  
    if(imesh_s_all==1.or.(imesh_s_all==0.and.nproc_id_global<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
      call exc_cor_ns
    end if
   
    call allgatherv_vlocal
    
    elp3(117)=get_wtime()
    elp3(127)=elp3(127)+elp3(117)-elp3(116)
  
    select case(iperiodic)
    case(0)
      call Total_Energy(psi)
    case(3)
      call Total_Energy(zpsi)
      do ik=k_sta,k_end
      do iob=1,iobnum
!$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          zpsi_tmp(ix,iy,iz,iob,ik)=zpsi(ix,iy,iz,iob,ik)
        end do
        end do
        end do
      end do
      end do
      call Total_Energy_periodic_scf(zpsi_tmp)
    end select
  
    elp3(118)=get_wtime()
    elp3(128)=elp3(128)+elp3(118)-elp3(117)
    elp3(131)=get_wtime()
  
    elp3(132)=get_wtime()
    elp3(142)=elp3(142)+elp3(132)-elp3(131)
    
    elp3(118)=get_wtime()

    if(iperiodic==0)then  
      call change_order(psi)
    end if
  
  else if(iscf_order==2)then

    select case(iperiodic)
    case(0)
      call Gram_Schmidt_ns
    case(3)
      call Gram_Schmidt_periodic
    end select

    if(Miter>iDiter_nosubspace_diag)then
      select case(iperiodic)
      case(0)
        call subspace_diag
      case(3)
        call subspace_diag_periodic
      end select
    end if

    select case(iperiodic)
    case(0)
      call Gram_Schmidt_ns
    case(3)
      call Gram_Schmidt_periodic
    end select

    if( amin_routine == 'cg' .or. (amin_routine == 'cg-diis' .and. Miter <= iDiterYBCG) ) then
      select case(iperiodic)
      case(0)
        select case(gscg)
        case('y')
          call sgscg(psi,iflag)
        case('n')
          call DTcg(psi,iflag)
        end select
      case(3)
        select case(gscg)
        case('y')
          call gscg_periodic(zpsi,iflag)
        case('n')
          call DTcg_periodic(zpsi,iflag)
        end select
      end select
    else if( amin_routine == 'diis' .or. amin_routine == 'cg-diis' ) then
      select case(iperiodic)
      case(0)
        call rmmdiis(psi)
      case(3)
        stop "rmmdiis method is not implemented for periodic systems."
      end select
    end if

    select case(iperiodic)
    case(0)
      call Gram_Schmidt_ns
    case(3)
      call Gram_Schmidt_periodic
    end select

    select case(iperiodic)
    case(0)
      call calc_density(psi)
    case(3)
      call calc_density(zpsi)
    end select

    select case(amixing)
      case ('simple')
        call simple_mixing(1.d0-rmixrate,rmixrate)
      case ('broyden')
        call buffer_broyden_ns(iter)
    end select
    
    if(imesh_s_all==1.or.(imesh_s_all==0.and.nproc_id_global<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
      call Hartree_ns
    end if
  
    if(imesh_s_all==1.or.(imesh_s_all==0.and.nproc_id_global<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
      call exc_cor_ns
    end if
   
    call allgatherv_vlocal
    
    select case(iperiodic)
    case(0)
      call calc_density(psi)
    case(3)
      call calc_density(zpsi)
    end select

    select case(iperiodic)
    case(0)
      call Total_Energy(psi)
    case(3)
      call Total_Energy(zpsi)
      do ik=k_sta,k_end
      do iob=1,iobnum
!$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          zpsi_tmp(ix,iy,iz,iob,ik)=zpsi(ix,iy,iz,iob,ik)
        end do
        end do
        end do
      end do
      end do
      call Total_Energy_periodic_scf(zpsi_tmp)
    end select
  end if

  select case(convergence)
    case('rho_dne')
      sum0=0.d0
!$OMP parallel do reduction(+:sum0) private(iz,iy,ix)
      do iz=ng_sta(3),ng_end(3) 
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        sum0=sum0+abs(rho(ix,iy,iz)-rho_stock(ix,iy,iz,1))
      end do
      end do
      end do
      call comm_summation(sum0,sum1,nproc_group_h)
      if(ispin==0)then
        sum1=sum1*Hvol/(dble(ifMST(1))*2.d0)
      else if(ispin==1)then
        sum1=sum1*Hvol/dble(ifMST(1)+ifMST(2))
      end if
    case('norm_rho','norm_rho_dng')
      sum0=0.d0
!$OMP parallel do reduction(+:sum0) private(iz,iy,ix)
      do iz=ng_sta(3),ng_end(3) 
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        sum0=sum0+(rho(ix,iy,iz)-rho_stock(ix,iy,iz,1))**2
      end do
      end do
      end do
      call comm_summation(sum0,sum1,nproc_group_h)
      if(convergence=='norm_rho_dng')then
        sum1=sum1/dble(lg_num(1)*lg_num(2)*lg_num(3))
      end if
    case('norm_pot','norm_pot_dng')
      sum0=0.d0
!$OMP parallel do reduction(+:sum0) private(iz,iy,ix)
      do iz=ng_sta(3),ng_end(3) 
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        sum0=sum0+(Vlocal(ix,iy,iz,1)-Vlocal_stock(ix,iy,iz,1))**2
      end do
      end do
      end do
      call comm_summation(sum0,sum1,nproc_group_h)
      if(convergence=='norm_pot_dng')then
        sum1=sum1/dble(lg_num(1)*lg_num(2)*lg_num(3))
      end if
  end select 

  if(comm_is_root(nproc_id_global)) then
    write(*,*) '-----------------------------------------------'
    select case(iperiodic)
    case(0)
      if(iflag_diisjump == 1) then
        write(*,'("Diisjump occured. Steepest descent was used.")')
      end if
      write(*,'(1x,"iter =",i6,5x,"Total Energy =",f19.8,5x,"Vh iteration =",i4)') Miter,Etot*2d0*Ry,iterVh
    case(3)
      write(*,'(1x,"iter =",i6,5x,"Total Energy =",f19.8,5x)') Miter,Etot*2d0*Ry
    end select
    do ik=1,num_kpoints_rd
      if(ik<=3)then
        if(iperiodic==3) write(*,*) "k=",ik
        do p5=1,(itotMST+3)/4
          p1=4*(p5-1)+1
          p2=4*p5 ; if ( p2 > itotMST ) p2=itotMST
          write(*,'(1x,4(i5,f15.4,2x))') (iob,esp(iob,ik)*2d0*Ry,iob=p1,p2)
        end do
        if(iperiodic==3) write(*,*) 
      end if
    end do
    select case(convergence)
      case('rho_dne')
        write(*,'("iter and int_x|rho_i(x)-rho_i-1(x)|dx/nelec     = ",i6,e15.8)') Miter,sum1
      case('norm_rho')
        write(*,'("iter and ||rho_i(ix)-rho_i-1(ix)||**2              = ",i6,e15.8)') Miter,sum1/a_B**6
      case('norm_rho_dng')
        write(*,'("iter and ||rho_i(ix)-rho_i-1(ix)||**2/(# of grids) = ",i6,e15.8)') Miter,sum1/a_B**6
      case('norm_pot')
        write(*,'("iter and ||Vlocal_i(ix)-Vlocal_i-1(ix)||**2              = ",i6,e15.8)') Miter,     &
                                                                         sum1*(2.d0*Ry)**2/a_B**6
      case('norm_pot_dng')
        write(*,'("iter and ||Vlocal_i(ix)-Vlocal_i-1(ix)||**2/(# of grids) = ",i6,e15.8)') Miter,     &
                                                                         sum1*(2.d0*Ry)**2/a_B**6
    end select
  end if 
  rNebox1=0.d0 
!$OMP parallel do reduction(+:rNebox1) private(iz,iy,ix)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    rNebox1=rNebox1+rho(ix,iy,iz)
  end do
  end do
  end do
  call comm_summation(rNebox1,rNebox2,nproc_group_global)
  if(comm_is_root(nproc_id_global))then
    write(*,*) "Ne=",rNebox2*Hvol
  end if

  elp3(119)=get_wtime()
  elp3(129)=elp3(129)+elp3(119)-elp3(118)
  elp3(130)=elp3(130)+elp3(119)-elp3(111)

if(ilsda==0)then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    rho_stock(ix,iy,iz,1)=rho(ix,iy,iz)
    Vlocal_stock(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)
  end do
  end do
  end do
else if(ilsda==1)then
!$OMP parallel do private(iz,iy,ix)
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
elp3(104)=get_wtime()

deallocate(idiis_sd)

if(icalcforce==1) then
  if(iperiodic==0) then
    call calc_force
  elseif(iperiodic==3) then
    if(comm_is_root(nproc_id_global)) write(*,*) &
      "This version is not supporting force-calculation with spatial-grid parallelization."
    stop
  end if
  if(comm_is_root(nproc_id_global))then
    write(*,*) "===== force ====="
    do iatom=1,MI
      select case(unit_system)
      case('au','a.u.')
        write(*,'(i6,3e16.8)') iatom,(rforce(ix,iatom),ix=1,3)
      case('A_eV_fs')
        write(*,'(i6,3e16.8)') iatom,(rforce(ix,iatom)*2.d0*Ry/a_B,ix=1,3)
      end select
    end do
  end if
end if

if(iflag_stopt==1) then
  call structure_opt_check(MI,istopt,istopt_tranc,rforce)
  if(istopt_tranc/=1) call structure_opt(MI,istopt,rforce,Rion)
  if(comm_is_root(nproc_id_global))then
    write(*,*) "atomic coordinate"
    do iatom=1,MI
      write(*,'(a3,3f16.8,i3,a3)') AtomName(Kion(iatom)), (Rion(jj,iatom)*ulength_from_au,jj=1,3), &
                                   Kion(iatom),flag_geo_opt_atom(iatom)
    end do
  end if
  if(istopt_tranc==1) then
    call structure_opt_fin
    exit Multigrid_Iteration
  end if
end if

end do Multigrid_Iteration
if(istopt_tranc==1)then
  exit Structure_Optimization_Iteration
end if
end do Structure_Optimization_Iteration

!---------------------------------------- Output

call write_eigen

if(out_psi=='y') then
  call writepsi
end if

if(out_dns=='y') then
  call writedns
end if

if(out_dos=='y') then
  call calc_dos
end if

if(out_pdos=='y') then
  call calc_pdos
end if

if(OC==2)then
  call prep_ini
end if

if(out_elf=='y')then
  allocate(elf(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),      &
               lg_sta(3):lg_end(3)))
  call calcELF
  call writeelf
  deallocate(elf)
end if

! LDA data
! subroutines in scf_data.f90
if ( OC==1.or.OC==2.or.OC==3 ) then
  call OUT_data
end if
elp3(105)=get_wtime()

! LDA information

if(comm_is_root(nproc_id_global)) then
  open(1,file=LDA_info)

  write(1,*) "Total number of iteration = ", Miter
  write(1,*)
  select case (ilsda)
  case(0)
    write(1,*) "Number of states = ", nstate
    write(1,*) "Number of electrons = ", ifMST(1)*2
  case(1)
    write(1,*) "Number of states = ", (nstate_spin(is),is=1,2)
    write(1,*) "Number of electrons = ", (nelec_spin(is),is=1,2)
  end select
  write(1,*)
  write(1,*) "Total energy (eV) = ", Etot*2d0*Ry
  write(1,*) "1-particle energies (eV)"
  select case (ilsda)
  case(0)
    do p5=1,(nstate+3)/4
      p1=4*(p5-1)+1
      p2=4*p5 ; if ( p2 > nstate ) p2=nstate
      write(1,'(1x,4(i5,f15.4,2x))') (iob,esp(iob,1)*2d0*Ry,iob=p1,p2)
    end do
  case(1)
    do is=1,2
      select case(is)
      case(1)
        write(1,*) "for up-spin"
        do p5=1,(nstate_spin(is)+3)/4
          p1=4*(p5-1)+1
          p2=4*p5 ; if ( p2 > nstate_spin(1) ) p2=nstate_spin(1)
          write(1,'(1x,4(i5,f15.4,2x))') (iob,esp(iob,1)*2d0*Ry,iob=p1,p2)
        end do
      case(2)
        write(1,*) "for down-spin"
        do p5=1,(nstate_spin(is)+3)/4
          p1=4*(p5-1)+1+nstate_spin(1)
          p2=4*p5+nstate_spin(1) ; if ( p2 > nstate_spin(1)+nstate_spin(2) ) p2=nstate_spin(1)+nstate_spin(2)
          write(1,'(1x,4(i5,f15.4,2x))') (iob-nstate_spin(1),esp(iob,1)*2d0*Ry,iob=p1,p2)
        end do
      end select
    end do
  end select
  write(1,*)

  do ii=1,ntmg
    write(1,'(1x,a,3f14.8)') "Size of the box (A) = ", rLsize(:,ii)*a_B
  end do

  write(1,'(1x,a,3f14.8)')   "Grid spacing (A)    = ", (Hgs(jj)*a_B,jj=1,3)
  write(1,*)
  write(1,'(1x,"Number of atoms = ",i8)') MI
  do ik=1,MKI
    write(1,'(1x,"iZatom(",i3,")     = ",i8)') ik, iZatom(ik)
  end do
  write(1,*)
  write(1,*) "Ref. and max angular momentum",      &
        " and pseudo-core radius of PP (A)"
  do ikoa=1,MKI
     write(1,'(1x,"(",i3,")  "," Ref, Max, Rps =",2i4,f8.3)')      &
                              ikoa,Lref(ikoa),Mlps(ikoa),Rps(ikoa)*a_B
  end do

  close(1)

end if
elp3(106)=get_wtime()

if(comm_is_root(nproc_id_global))then
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

call finalize_xc(xc_func)

END subroutine Real_Space_DFT

!=======================================================================
!========================================= Grid generation and labelling

SUBROUTINE init_mesh
use salmon_parallel, only: nproc_id_global, nproc_size_global
use salmon_communication, only: comm_is_root
use inputoutput, only: iperiodic
use global_variables_scf
implicit none

real(8) :: rLsize1(3)

if(comm_is_root(nproc_id_global))      &
    print *,"----------------------------------- init_mesh"

rLsize1(:)=rLsize(:,img)
call setlg(lg_sta,lg_end,lg_num,ista_Mx_ori,iend_Mx_ori,inum_Mx_ori,    &
           Hgs,Nd,rLsize1,imesh_oddeven,iperiodic)
call check_fourier

allocate(ista_Mxin(3,0:nproc_size_global-1),iend_Mxin(3,0:nproc_size_global-1))
allocate(inum_Mxin(3,0:nproc_size_global-1))

call setmg(mg_sta,mg_end,mg_num,ista_Mxin,iend_Mxin,inum_Mxin,  &
           lg_sta,lg_num,nproc_size_global,nproc_id_global,nproc_Mxin,nproc_k,nproc_ob,isequential)

if(comm_is_root(nproc_id_global)) write(*,*) "Mx     =", iend_Mx_ori

return

END SUBROUTINE init_mesh

