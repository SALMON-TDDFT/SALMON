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
module read_pslfile_sub

  use scf_data
  
  integer,parameter :: Nlps=3, Nlmps=16
  
  integer,allocatable :: Mlps0(:)
  
  !integer :: maxMr
  real(8) :: rmin_step
  real(8) :: rmaxRps
  
  real(8) :: rPC,rRC(0:Nlps)
  
  integer,allocatable :: Mr(:)
  real(8),allocatable :: step(:)
  
  real(8), allocatable :: upp_f(:,:,:)
  real(8), allocatable :: rhopp_f(:,:)
  real(8), allocatable :: vpp_f(:,:,:)
  
  real(8), allocatable :: rad_f(:,:)
  
  contains
  !==================================================================================================
  subroutine read_pslfile
    use salmon_communication, only: comm_is_root
    use salmon_pp, only: init_pp
    implicit none
    integer :: ak,i,ll
    
    allocate( Mlps0(MKI) )
    allocate( Mr(MKI) )
    allocate( step(MKI) )
    
    allocate( Zps(MKI) )
    allocate( Rps(MKI) )
    allocate( Mass(MKI) )
    
    call init_pp(pp,Nrmax,Lmax,flag_nlcc)
    call init_mps(ppg)
    call init_mps(ppg_all)
    
    allocate(upp_f(0:Nrmax,0:Nlps,MKI))
    allocate(rhopp_f(0:Nrmax,MKI))
    allocate(vpp_f(0:Nrmax,0:Nlps,MKI))
    allocate(rad_f(0:Nrmax,MKI) )
    
    call input_pp(pp,harray(1,1),harray(2,1),harray(3,1))
 
    Mass(1:MKI)=pp%rmass(1:MKI)
  
    do ak=1,MKI
      Mr(ak)=pp%mr(ak)
      Mlps0(ak)=pp%mlps(ak)
      do ll=0,Mlps0(ak)
        do i=0,Mr(ak)
          upp_f(i,ll,ak)=pp%upp_f(i,ll,ak)
          vpp_f(i,ll,ak)=pp%vpp_f(i,ll,ak)
        end do
      end do
      do i=1,Mr(ak)
        rad_f(i-1,ak)=pp%rad(i,ak)
      end do
    end do
  
    Zps(1:MKI)=pp%zps(1:MKI)
    Rps(1:MKI)=pp%rps(1:MKI)
    rmaxRps=maxval(Rps(1:MKI))
  
    step(1:MKI)=1.d8
    do ak=1,MKI
      select case(ps_format(ak))
      case('KY')
        step(ak)=pp%rad(2,ak)-pp%rad(1,ak)
      case('ABINIT')
        step(ak)=0.01d0/a_B
      case('ABINITFHI','FHI')
        step(ak)=pp%rad(2,ak)/pp%rad(1,ak)
      end select
    end do
    rmin_step=minval(step(1:MKI))
  
    maxlm=0
    do ak=1,MKI
      if(Mlps(ak)>maxlm) maxlm=Mlps(ak)
    end do
    maxlm=(maxlm+1)**2

    do ak=1,MKI
      select case(ps_format(ak))
      case('KY')
        ipsfileform(ak) = n_Yabana_Bertsch_psformat
      case('ABINIT')
        ipsfileform(ak) = n_ABINIT_psformat
      case('ABINITFHI')
        ipsfileform(ak) = n_ABINITFHI_psformat
      case('FHI')
        ipsfileform(ak) = n_FHI_psformat
      end select
    end do
    
    do ak=1,MKI
      select case( iZatom(ak) )
        case(1) ; Atomname(ak)='H'  
        case(2) ; Atomname(ak)='He' 
        case(3) ; Atomname(ak)='Li' 
        case(4) ; Atomname(ak)='Be' 
        case(5) ; Atomname(ak)='B'  
        case(6) ; Atomname(ak)='C'  
        case(7) ; Atomname(ak)='N'  
        case(8) ; Atomname(ak)='O'  
        case(9) ; Atomname(ak)='F'  
        case(10); Atomname(ak)='Ne' 
        case(11); Atomname(ak)='Na' 
        case(12); Atomname(ak)='Mg' 
        case(13); Atomname(ak)='Al' 
        case(14); Atomname(ak)='Si' 
        case(15); Atomname(ak)='P'  
        case(16); Atomname(ak)='S'  
        case(17); Atomname(ak)='Cl' 
        case(18); Atomname(ak)='Ar' 
        case(19); Atomname(ak)='K'  
        case(20); Atomname(ak)='Ca' 
        case(21); Atomname(ak)='Sc' 
        case(22); Atomname(ak)='Ti' 
        case(23); Atomname(ak)='V'  
        case(24); Atomname(ak)='Cr' 
        case(25); Atomname(ak)='Mn' 
        case(26); Atomname(ak)='Fe' 
        case(27); Atomname(ak)='Co' 
        case(28); Atomname(ak)='Ni' 
        case(29); Atomname(ak)='Cu' 
        case(30); Atomname(ak)='Zn' 
        case(31); Atomname(ak)='Ga' 
        case(32); Atomname(ak)='Ge' 
        case(33); Atomname(ak)='As' 
        case(34); Atomname(ak)='Se' 
        case(35); Atomname(ak)='Br' 
        case(36); Atomname(ak)='Kr' 
        case(47); Atomname(ak)='Ag' 
        case(79); Atomname(ak)='Au' 
      end select
    end do

  end subroutine read_pslfile
    
end module read_pslfile_sub
