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
module structure_opt_sub
  implicit none
  real(8),allocatable :: r1_opt(:),r2_opt(:)
  real(8),allocatable :: H_opt(:,:),H_opt_temp(:,:)
contains
  !=======================================================================
  !==============================================================initilize
  subroutine structure_opt_ini(iMI_opt)
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    implicit none
    integer,intent(in) :: iMI_opt
    allocate(r1_opt(3*iMI_opt),r2_opt(3*iMI_opt))
    allocate(H_opt(3*iMI_opt,3*iMI_opt),H_opt_temp(3*iMI_opt,3*iMI_opt))
    r1_opt(:)=0.0d0; r2_opt(:)=0.0d0;
    H_opt(:,:)=0.0d0; H_opt_temp(:,:)=0.0d0;
    if(comm_is_root(nproc_id_global))then
      write(*,*) "===== Grand State Optimization Start ====="
      write(*,*) "       (Quasi-Newton method using Force only)       "
    end if
  end subroutine structure_opt_ini
  !=======================================================================
  !======================================================convergence check
  subroutine structure_opt_check(iMI_opt,iopt,itranc,rforce_opt)
    use salmon_global, only: convrg_opt_fmax,unit_system,flag_geo_opt_atom
    use salmon_parallel, only: nproc_id_global,nproc_group_global
    use salmon_communication, only: comm_is_root,comm_bcast
    implicit none
    integer,intent(in) :: iMI_opt,iopt
    integer,intent(inout) :: itranc
    real(8),intent(in) :: rforce_opt(3,iMI_opt)
    real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0
    integer :: iatom,iatom_count
    real(8) :: fabs,fmax,fave
    fmax=0.0d0; fave= 0d0;
    iatom_count=0
    do iatom=1,iMI_opt
      if(flag_geo_opt_atom(iatom)=='y') then
        iatom_count=iatom_count+1
        fabs=rforce_opt(1,iatom)**2.0d0+rforce_opt(2,iatom)**2.0d0+rforce_opt(3,iatom)**2.0d0
        fave=fave+fabs
        if(fabs>=fmax) fmax=fabs
      end if
    enddo
    select case(unit_system)
    case('au','a.u.')
      fmax = sqrt(fmax)
      fave = sqrt(fave/iatom_count)
    case('A_eV_fs')
      fmax = sqrt(fmax)*2.d0*Ry/a_B
      fave = sqrt(fave/iatom_count)*2.d0*Ry/a_B
    end select
    if(comm_is_root(nproc_id_global))then
      write(*,*) " Max-force=",fmax, "  Mean-force=",fave
      write(*,*) "==================================================="
      write(*,*) "Quasi-Newton Optimization Step = ", iopt
      if(fmax<=convrg_opt_fmax) itranc=1;
    end if
    call comm_bcast(itranc,nproc_group_global)
  end subroutine structure_opt_check
  !=======================================================================
  !===========================================================optimization
  subroutine structure_opt(iMI_opt,iopt,rforce_opt,Rion_opt)
    use salmon_global, only: flag_geo_opt_atom
    use salmon_parallel, only: nproc_group_global
    use salmon_communication, only: comm_bcast
    implicit none
    integer,intent(in) :: iMI_opt,iopt
    real(8),intent(in) :: rforce_opt(3,iMI_opt)
    real(8),intent(inout) :: Rion_opt(3,iMI_opt)
    real(8), parameter :: alpha=1.0d0,theta_opt=1.0d0  !theta_opt=0.0d0:DFP,theta_opt=1.0d0:BFGS in Quasi_Newton method
    integer :: ii,ij,icount_opt,iatom
    real(8) :: const1_opt,const2_opt
    real(8) :: rforce_1d(3*iMI_opt),del_Rion_1d(3*iMI_opt),optmat_1d(3*iMI_opt)
    real(8) :: del_Rion(3,iMI_opt)
    real(8) :: optmat1_2d((3*iMI_opt),(3*iMI_opt)),optmat2_2d((3*iMI_opt),(3*iMI_opt)),optmat3_2d((3*iMI_opt),(3*iMI_opt))
    !transrate rforce to rforce_1d
    icount_opt=1
    do iatom=1,iMI_opt
      do ii=1,3
        rforce_1d(icount_opt)=rforce_opt(ii,iatom)
        icount_opt=icount_opt+1
      end do
    end do
    if(iopt==1)then
      !update H_opt
      do ii=1,(3*iMI_opt)
        do ij=1,(3*iMI_opt)
          if(ii==ij)then
            H_opt(ii,ij)=1.0d0
          else
            H_opt(ii,ij)=0.0d0
          end if
          H_opt_temp(ii,ij)=H_opt(ii,ij)
        end do
      end do
    else
      !update r2_opt
      r2_opt=-(rforce_1d-r2_opt)
      !prepare const and matrix
      call dgemm('n','n',1,1,(3*iMI_opt),1.0d0,r1_opt,1,r2_opt,(3*iMI_opt),0.0d0,const1_opt,1)
      call dgemm('n','n',(3*iMI_opt),1,(3*iMI_opt),1.0d0,H_opt,(3*iMI_opt),r2_opt,(3*iMI_opt),0.0d0,optmat_1d,(3*iMI_opt))
      call dgemm('n','n',1,1,(3*iMI_opt),1.0d0,r2_opt,1,optmat_1d,(3*iMI_opt),0.0d0,const2_opt,1)
      call dgemm('n','n',(3*iMI_opt),(3*iMI_opt),1,1.0d0,r1_opt,(3*iMI_opt),r1_opt,1,0.0d0,optmat1_2d,(3*iMI_opt))
      !update H_opt
      H_opt=H_opt_temp+((const1_opt+theta_opt*const2_opt)/(const1_opt**2.0d0))*optmat1_2d
      if(theta_opt==0.0d0)then
        !theta_opt=0.0d0:DFP
        call dgemm('n','n',(3*iMI_opt),(3*iMI_opt),1,1.0d0,optmat_1d,(3*iMI_opt),optmat_1d,1,0.0d0,optmat2_2d,(3*iMI_opt))
        H_opt=H_opt-(1.0d0/const2_opt)*optmat2_2d
      elseif(theta_opt==1.0d0)then
        !theta_opt=1.0d0:BFGS
        call dgemm('n','n',(3*iMI_opt),(3*iMI_opt),1,1.0d0,optmat_1d,(3*iMI_opt),r1_opt,1,0.0d0,optmat2_2d,(3*iMI_opt))
        call dgemm('n','n',(3*iMI_opt),(3*iMI_opt),1,1.0d0,r1_opt,(3*iMI_opt),optmat_1d,1,0.0d0,optmat3_2d,(3*iMI_opt))
        H_opt=H_opt-(theta_opt/const1_opt)*(optmat2_2d+optmat3_2d)
      endif
      !update H_opt_temp
      H_opt_temp=H_opt
    end if
    !update del_Rion_1d and del_Rion
    del_Rion_1d(:)=0.0d0;
    call dgemm('n','n',(3*iMI_opt),1,(3*iMI_opt),1.0d0,H_opt,(3*iMI_opt),rforce_1d,(3*iMI_opt),0.0d0,del_Rion_1d,(3*iMI_opt))
    do iatom=1,iMI_opt
      del_Rion(1:3,iatom)=del_Rion_1d((1+3*(iatom-1)):(3+3*(iatom-1)))
    end do
    !update r1_opt,r2_opt
    r1_opt=alpha*del_Rion_1d
    r2_opt=rforce_1d
    !update Rion
    do iatom=1,iMI_opt
      if(flag_geo_opt_atom(iatom)=='y') then
        Rion_opt(1:3,iatom)=Rion_opt(1:3,iatom)+alpha*del_Rion(1:3,iatom)
      end if
    end do
    call comm_bcast(Rion_opt,nproc_group_global)
  end subroutine structure_opt
  !=======================================================================
  !===============================================================finilize
  subroutine structure_opt_fin
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    implicit none
    deallocate(r1_opt,r2_opt)
    deallocate(H_opt,H_opt_temp)
    if(comm_is_root(nproc_id_global)) write(*,*) "Optimization Converged"
  end subroutine structure_opt_fin
end module structure_opt_sub