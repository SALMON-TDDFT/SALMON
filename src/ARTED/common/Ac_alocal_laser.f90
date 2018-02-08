! This file is for "alocal_laser" option -- just trial version (by AY)
! This is an option for applying electric field only at specific domain or atoms.
! Usually this is not used but for something specific case, you may want to use
! for e.x. local excitation.
! Currently, space coordinate dependent vector potential A(r,t) is implemented
! for KS equation and total energy calculation. (NOT for current, force etc)
! It works for calculating dielectric function with "impulse" option, 
! as A(r,t) is only used in the first time step, doesn't affect on current so much.
!
! Should include A(r,t) in future:
!  - Energy (if you only see the initial and final state, not necessary to fix)
!  - Force from non-local pseudopotential term 
!  - Non-local pseudopotential term in KS equation
!  - Current
!
module Ac_alocal_laser
  implicit none

contains

  !-----------------------------
  subroutine init_Ac_alocal
  use Global_Variables
  use salmon_parallel, only: nproc_group_global, nproc_id_global
  use salmon_communication, only: comm_bcast, comm_is_root
  implicit none
  logical :: flag_use_gauss_for_weight_ion
  integer :: i,j,NI_alocal
  integer :: ix,iy,iz,ia,ia2,iatom(NI)
  real(8) :: x,y,z,r,tmpx,tmpy,tmpz,tmp1,radi(NI),sgm(NI),ave_weight_Ac_alocal
  character(10) :: ckeyword

  allocate(Ac_ext_al(-1:Nt+1,3),Ac_tot_al(-1:Nt+1,3))
  Ac_tot_al = Ac_tot
  Ac_ext_al = Ac_ext
  Ac_tot    = 0d0
  Ac_ext    = 0d0

  if(comm_is_root(nproc_id_global))then
     open(898,file='input_Ac_alocal.dat')

     read(898,*) ckeyword
     if(ckeyword=="atoms") then

        read(898,*) NI_alocal
        do i=1,NI_alocal
           read(898,*) iatom(i),radi(i),sgm(i)   !unit=[au]
        enddo

     else if(ckeyword=="grids") then
        write(*,*) "Error: keyword in input_Ac_alocal.dat is not supported now"
        stop
     else
        write(*,*) "Error: keyword in input_Ac_alocal.dat is bad"
        stop
     endif
     close(898)
  end if

  call comm_bcast(ckeyword, nproc_group_global)
  call comm_bcast(NI_alocal,nproc_group_global)

  allocate(weight_Ac_alocal(NL),weight_Ac_alocal_ion(NI),divA_al(NL))
  weight_Ac_alocal(:)    =0d0
  weight_Ac_alocal_ion(:)=0d0

  if(ckeyword=="atoms") then

     if(comm_is_root(nproc_id_global))then
        write(*,*) "# input_Ac_alocal.dat was read:"
        write(*,*) "    ",trim(ckeyword)
        do i=1,NI_alocal
           write(*,*) "  ",iatom(i),real(radi(i)),real(sgm(i))
        enddo
     endif

     call comm_bcast(iatom,nproc_group_global)
     call comm_bcast(radi, nproc_group_global)
     call comm_bcast(sgm,  nproc_group_global)

!$omp parallel
     do j=1,NI_alocal
        ia=iatom(j)
        do ix=-2,2
        do iy=-2,2
        do iz=-2,2
           tmpx = Rion(1,ia)+ix*aLx
           tmpy = Rion(2,ia)+iy*aLy
           tmpz = Rion(3,ia)+iz*aLz
!$omp do private(i,x,y,z,r,tmp1)
           do i=1,NL
              x=Lx(i)*Hx-tmpx
              y=Ly(i)*Hy-tmpy
              z=Lz(i)*Hz-tmpz
              r=sqrt(x*x+y*y+z*z)
              if(r<radi(j)) then
                 weight_Ac_alocal(i)=1d0
              else
                 tmp1 = exp(-(r-radi(j))**2d0/(sgm(j)*sgm(j)))
                 weight_Ac_alocal(i)= weight_Ac_alocal(i) + tmp1
              endif
           enddo
!$omp end do
           do ia2=1,NI
              x=Rion(1,ia2)-tmpx
              y=Rion(2,ia2)-tmpy
              z=Rion(3,ia2)-tmpz
              r=sqrt(x*x+y*y+z*z)
              if(r<radi(j)) then
                 weight_Ac_alocal_ion(ia2)=1d0
              else
                 tmp1 = exp(-(r-radi(j))**2d0/(sgm(j)*sgm(j)))
                 weight_Ac_alocal_ion(ia2)= weight_Ac_alocal_ion(ia2) + tmp1
              endif
           enddo
        enddo
        enddo
        enddo
     end do
!$omp end parallel

     do i=1,NL
        if(weight_Ac_alocal(i).gt.1d0) weight_Ac_alocal(i)=1d0
     enddo
     do ia=1,NI
        if(weight_Ac_alocal_ion(ia).gt.1d0) weight_Ac_alocal_ion(ia)=1d0
     enddo

!     flag_use_gauss_for_weight_ion=.true.
     flag_use_gauss_for_weight_ion=.false.
     if(.not. flag_use_gauss_for_weight_ion)then
        weight_Ac_alocal_ion(:)=0d0
        do j=1,NI_alocal
           ia=iatom(j)
           weight_Ac_alocal_ion(ia)=1d0
        enddo
     endif


     !log
     ave_weight_Ac_alocal = sum(weight_Ac_alocal(:))/dble(NL)
     if(comm_is_root(nproc_id_global)) then
       write(*,*)"  average spatial weight for local_laser=",real(ave_weight_Ac_alocal)
       write(*,*)"  weight on atoms for local_laser="
       do ia=1,NI
          write(*,'(4X,i6,f10.5)') ia,weight_Ac_alocal_ion(ia)
       enddo
     endif

     call write_weight_cube_alocal

  endif


  call prep_RT_Ac_alocal_laser(0)


  contains

    subroutine write_weight_cube_alocal
      implicit none
      integer :: i, j, ix, iy, iz, ifh
      real(8) :: r
      character(256) :: file_name_alocal

      ifh= 505
      file_name_alocal = "weight_alocal.cube"
      open(ifh,file=trim(file_name_alocal),status="unknown")
      
      write(ifh, '(A)') "# SALMON"
      write(ifh, '(A)') "# COMMENT"
      write(ifh, '(I5,3(F12.6))') NI, 0.00, 0.00, 0.00
      write(ifh, '(I5,3(F12.6))') NLx, Hx, 0.00, 0.00
      write(ifh, '(I5,3(F12.6))') NLy, 0.00, Hy, 0.00
      write(ifh, '(I5,3(F12.6))') NLz, 0.00, 0.00, Hz
      
      do i=1, NI
         write(ifh,'(I5,4(F12.6))') Zatom(Kion(i)), 0d0, (Rion(j,i),j=1,3)
      end do
  
      ! Gaussian .cube file (x-slowest index, z-fastest index)
      i=1
      do ix=0, NLx-1
      do iy=0, NLy-1
      do iz=0, NLz-1
        r = weight_Ac_alocal(Lxyz(ix,iy,iz))

        if(mod(i,6)==0) then
           write(ifh,10) r
        else
           write(ifh,10,advance='no') r
        endif
        i=i+1
      end do
      end do
      end do
      close(ifh)

10    format(ES12.4)
      return
    end subroutine write_weight_cube_alocal

  end subroutine init_Ac_alocal

  !---------------------------------------
  subroutine prep_RT_Ac_alocal_laser(it)
  use Global_Variables
  implicit none
  integer :: ik,it
  real(8),allocatable:: tmpAc2(:)

  Ac_al_amp(:) = Ac_tot_al(it,:)

  allocate(tmpAc2(NL))
  if(.not. allocated(Ac1x_al)) then
     allocate(Ac1x_al(NL),Ac1y_al(NL),Ac1z_al(NL),Ac2_al(NL,NK_s:NK_e))
  endif
  Ac1x_al(:) = Ac_al_amp(1)*weight_Ac_alocal(:)
  Ac1y_al(:) = Ac_al_amp(2)*weight_Ac_alocal(:)
  Ac1z_al(:) = Ac_al_amp(3)*weight_Ac_alocal(:)
  Ac2_al_amp = sum(Ac_al_amp(:)**2d0)
  tmpAc2(:)  = 0.5d0 * Ac2_al_amp * weight_Ac_alocal(:)**2d0
  do ik=NK_s,NK_e
     Ac2_al(:,ik) = tmpAc2(:) + sum(kAc0(ik,:)*Ac_al_amp(:))*weight_Ac_alocal(:)
  enddo
  nabt_al( 1: 4)= nabx(1:4)
  nabt_al( 5: 8)= naby(1:4)
  nabt_al( 9:12)= nabz(1:4)

  deallocate(tmpAc2)
  end subroutine prep_RT_Ac_alocal_laser
  !----------------------------------------------------------
  subroutine get_Eelemag_FionAc_alocal_laser(it)
    use Global_Variables
    implicit none
    integer :: i,it,ia
    real(8) :: Etot(3),E_tot_xyz(3,NL)

    !cancel out Eelemag from Eall added in tddft_sc subroutine 
    !but actually Eelemag=0 for local_laser
    Eall = Eall - Eelemag  

    !calc Eelemag
    Eelemag = 0d0
    Etot(:) = -(Ac_tot_al(it+1,:)-Ac_tot_al(it-1,:))/(2d0*dt)
    do i=1,NL
       E_tot_xyz(:,i)= Etot(:)*weight_Ac_alocal(i)
       Eelemag = Eelemag + sum(E_tot_xyz(:,i)**2)
    enddo
   !Eelemag=aLxyz*sum(E_tot(it,:)**2)/(8.d0*Pi)
    Eelemag = Eelemag/dble(NL)/(8.d0*Pi)*aLxyz
    Eall = Eall + Eelemag

    !calc force from local field
    do ia=1,NI
      FionAc(:,ia) = Zps(Kion(ia))*Etot(:)*weight_Ac_alocal_ion(ia)
    enddo

  end subroutine get_Eelemag_FionAc_alocal_laser
  !----------------------------------------------------------
  subroutine hpsi1_RT_stencil_add_Ac_alocal(B,Cx,Cy,Cz,D,E,F)
  use global_variables, only: NLx,NLy,NLz,zI
#ifndef ARTED_DOMAIN_POWER_OF_TWO
  use opt_variables, only: modx, mody, modz
#endif
  use opt_variables, only: PNLx,PNLy,PNLz
  implicit none
  real(8),   intent(in)   :: D(12)
  real(8),   intent(in)   :: Cx(0: NLz-1,0: NLy-1,0: NLx-1)
  real(8),   intent(in)   :: Cy(0: NLz-1,0: NLy-1,0: NLx-1)
  real(8),   intent(in)   :: Cz(0: NLz-1,0: NLy-1,0: NLx-1)
  real(8),   intent(in)   ::  B(0: NLz-1,0: NLy-1,0: NLx-1)
  complex(8),intent(in)   ::  E(0:PNLz-1,0:PNLy-1,0:PNLx-1)
  complex(8),intent(inout)::  F(0:PNLz-1,0:PNLy-1,0:PNLx-1)

  integer    :: ix,iy,iz
  complex(8) :: w1,w2,w3,w,v1,v2,v3,v

#ifdef ARTED_DOMAIN_POWER_OF_TWO
# define IDX(dt) iz,iy,and(ix+(dt)+NLx,NLx-1)
# define IDY(dt) iz,and(iy+(dt)+NLy,NLy-1),ix
# define IDZ(dt) and(iz+(dt)+NLz,NLz-1),iy,ix
#else
# define IDX(dt) iz,iy,modx(ix+(dt)+NLx)
# define IDY(dt) iz,mody(iy+(dt)+NLy),ix
# define IDZ(dt) modz(iz+(dt)+NLz),iy,ix
#endif

  do ix=0,NLx-1
  do iy=0,NLy-1
  do iz=0,NLz-1

    w1=(D( 9)*(E(IDZ(1))-E(IDZ(-1))) &
    &  +D(10)*(E(IDZ(2))-E(IDZ(-2))) &
    &  +D(11)*(E(IDZ(3))-E(IDZ(-3))) &
    &  +D(12)*(E(IDZ(4))-E(IDZ(-4))))

    w2=(D( 5)*(E(IDY(1))-E(IDY(-1))) &
    &  +D( 6)*(E(IDY(2))-E(IDY(-2))) &
    &  +D( 7)*(E(IDY(3))-E(IDY(-3))) &
    &  +D( 8)*(E(IDY(4))-E(IDY(-4))))

    w3=(D( 1)*(E(IDX(1))-E(IDX(-1))) &
    &  +D( 2)*(E(IDX(2))-E(IDX(-2))) &
    &  +D( 3)*(E(IDX(3))-E(IDX(-3))) &
    &  +D( 4)*(E(IDX(4))-E(IDX(-4))))

    v1=(D( 9)*(Cz(IDZ(1))-Cz(IDZ(-1))) &
    &  +D(10)*(Cz(IDZ(2))-Cz(IDZ(-2))) &
    &  +D(11)*(Cz(IDZ(3))-Cz(IDZ(-3))) &
    &  +D(12)*(Cz(IDZ(4))-Cz(IDZ(-4))))

    v2=(D( 5)*(Cy(IDY(1))-Cy(IDY(-1))) &
    &  +D( 6)*(Cy(IDY(2))-Cy(IDY(-2))) &
    &  +D( 7)*(Cy(IDY(3))-Cy(IDY(-3))) &
    &  +D( 8)*(Cy(IDY(4))-Cy(IDY(-4))))

    v3=(D( 1)*(Cx(IDX(1))-Cx(IDX(-1))) &
    &  +D( 2)*(Cx(IDX(2))-Cx(IDX(-2))) &
    &  +D( 3)*(Cx(IDX(3))-Cx(IDX(-3))) &
    &  +D( 4)*(Cx(IDX(4))-Cx(IDX(-4))))

    w = Cz(iz,iy,ix)*w1 + Cy(iz,iy,ix)*w2 + Cx(iz,iy,ix)*w3
    v = 0.5d0*(v1+v2+v3)*E(iz,iy,ix)
    F(iz,iy,ix) = F(iz,iy,ix) + B(iz,iy,ix)*E(iz,iy,ix) - zI*(w+v)
  end do
  end do
  end do


  end subroutine hpsi1_RT_stencil_add_Ac_alocal

  !---------------------------------------------------------------
  subroutine total_energy_stencil_add_Ac_alocal(B,Cx,Cy,Cz,D,E,F)
  use global_variables, only: NLx,NLy,NLz,zI
#ifndef ARTED_DOMAIN_POWER_OF_TWO
  use opt_variables, only: modx, mody, modz
#endif
  implicit none
  real(8),   intent(in)  :: D(12)
  real(8),   intent(in)  :: Cx(0:NLz-1,0:NLy-1,0:NLx-1)
  real(8),   intent(in)  :: Cy(0:NLz-1,0:NLy-1,0:NLx-1)
  real(8),   intent(in)  :: Cz(0:NLz-1,0:NLy-1,0:NLx-1)
  real(8),   intent(in)  ::  B(0:NLz-1,0:NLy-1,0:NLx-1)
  complex(8),intent(in)  ::  E(0:NLz-1,0:NLy-1,0:NLx-1)
  complex(8),intent(out) ::  F

  integer    :: ix,iy,iz
  complex(8) :: w1,w2,w3,w,v1,v2,v3,v,z

  F = 0d0

#ifdef ARTED_DOMAIN_POWER_OF_TWO
# define IDX(dt) iz,iy,and(ix+(dt)+NLx,NLx-1)
# define IDY(dt) iz,and(iy+(dt)+NLy,NLy-1),ix
# define IDZ(dt) and(iz+(dt)+NLz,NLz-1),iy,ix
#else
# define IDX(dt) iz,iy,modx(ix+(dt)+NLx)
# define IDY(dt) iz,mody(iy+(dt)+NLy),ix
# define IDZ(dt) modz(iz+(dt)+NLz),iy,ix
#endif

  do ix=0,NLx-1
  do iy=0,NLy-1
  do iz=0,NLz-1

    w1=(D( 9)*(E(IDZ(1))-E(IDZ(-1))) &
    &  +D(10)*(E(IDZ(2))-E(IDZ(-2))) &
    &  +D(11)*(E(IDZ(3))-E(IDZ(-3))) &
    &  +D(12)*(E(IDZ(4))-E(IDZ(-4))))

    w2=(D( 5)*(E(IDY(1))-E(IDY(-1))) &
    &  +D( 6)*(E(IDY(2))-E(IDY(-2))) &
    &  +D( 7)*(E(IDY(3))-E(IDY(-3))) &
    &  +D( 8)*(E(IDY(4))-E(IDY(-4))))

    w3=(D( 1)*(E(IDX(1))-E(IDX(-1))) &
    &  +D( 2)*(E(IDX(2))-E(IDX(-2))) &
    &  +D( 3)*(E(IDX(3))-E(IDX(-3))) &
    &  +D( 4)*(E(IDX(4))-E(IDX(-4))))

    v1=(D( 9)*(Cz(IDZ(1))-Cz(IDZ(-1))) &
    &  +D(10)*(Cz(IDZ(2))-Cz(IDZ(-2))) &
    &  +D(11)*(Cz(IDZ(3))-Cz(IDZ(-3))) &
    &  +D(12)*(Cz(IDZ(4))-Cz(IDZ(-4))))

    v2=(D( 5)*(Cy(IDY(1))-Cy(IDY(-1))) &
    &  +D( 6)*(Cy(IDY(2))-Cy(IDY(-2))) &
    &  +D( 7)*(Cy(IDY(3))-Cy(IDY(-3))) &
    &  +D( 8)*(Cy(IDY(4))-Cy(IDY(-4))))

    v3=(D( 1)*(Cx(IDX(1))-Cx(IDX(-1))) &
    &  +D( 2)*(Cx(IDX(2))-Cx(IDX(-2))) &
    &  +D( 3)*(Cx(IDX(3))-Cx(IDX(-3))) &
    &  +D( 4)*(Cx(IDX(4))-Cx(IDX(-4))))

    w = Cz(iz,iy,ix)*w1 + Cy(iz,iy,ix)*w2 + Cx(iz,iy,ix)*w3
    v = 0.5*(v1+v2+v3)*E(iz,iy,ix)

    z = B(iz,iy,ix)*E(iz,iy,ix) - zI*(w+v)
    F = F + conjg(E(iz,iy,ix)) * z

  end do
  end do
  end do

  end subroutine total_energy_stencil_add_Ac_alocal
end module Ac_alocal_laser
