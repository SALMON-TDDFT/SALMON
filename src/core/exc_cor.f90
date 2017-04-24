subroutine core_exc_cor(xc, ispin, cval, nx, ny, nz, dv, rho, rho_s, tau, rj, grho, rlrho, rho_nlcc, vxc, vxc_s, exc, tot_exc)
  implicit none
  character(20), intent(in) :: xc
  integer, intent(in) :: ispin ! ispin=0(nonmag) =1(mag)
  real(8), intent(in) :: cval
  integer, intent(in) :: nx, ny, nz
  real(8), intent(in) :: dv
  real(8), intent(in) :: rho(nx, ny, nz)
  real(8), intent(in) :: rho_s(nx, ny, nz, 2)
  real(8), intent(in) :: tau(nx, ny, nz)
  real(8), intent(in) :: rj(nx, ny, nz, 3)
  real(8), intent(in) :: grho(nx, ny, nz, 3)
  real(8), intent(in) :: rlrho(nx, ny, nz)
  real(8), intent(in) :: rho_nlcc(nx, ny, nz)
  real(8), intent(out) :: vxc(nx, ny, nz) ! lda
  real(8), intent(out) :: vxc_s(nx, ny, nz, 2) !lsda
  real(8), intent(out) :: exc(nx, ny, nz)
  real(8), intent(out) :: tot_exc ! sum(exc)
  
  select case(xc)
  case("pz")
    if (ispin == 0) then
      call core_exc_cor_pz() !lda
    else
      call core_exc_cor_pz_spin() !lsda
    endif
  case("pbe")
    call core_exc_cor_pbe()
  case("lda")
    call core_exc_cor_lda()
  case("tbmbj")
    call core_exc_cor_tbmbj()
  case("tpss")
    call core_exc_cor_tpss()
  end select
  
  return
  
contains
  
  
  Subroutine core_exc_cor_pz()
    implicit none
    real(8),parameter :: Pi=3.141592653589793d0
    real(8),parameter :: gammaU=-0.1423d0,beta1U=1.0529d0
    real(8),parameter :: beta2U=0.3334d0,AU=0.0311d0,BU=-0.048d0
    real(8),parameter :: CU=0.002d0,DU=-0.0116d0
    real(8),parameter :: const = 0.75d0*(3d0/(2d0*pi))**(2d0/3d0)
    real(8) :: trho,rs,rssq,rsln
    real(8) :: de_xc_drho, e_xc
    integer :: ix, iy, iz
    
    tot_exc = 0.0
    do iz = 1, nz
      do iy = 1, ny
        do ix = 1, nx
          trho = rho(ix, iy, iz) + 1d-10
          rs=(3d0/(4*Pi*trho))**(1d0/3d0)
          e_xc=-const/rs
          de_xc_drho=e_xc/(3*trho)
          if (rs>1d0) then
            rssq=sqrt(rs)
            e_xc=e_xc+gammaU/(1d0+beta1U*rssq+beta2U*rs)
            de_xc_drho=de_xc_drho+gammaU*(0.5d0*beta1U*rssq+beta2U*rs)/(3*trho)/(1+beta1U*rssq+beta2U*rs)**2
          else
            rsln=log(rs)
            e_xc=e_xc+AU*rsln+BU+CU*rs*rsln+DU*rs
            de_xc_drho=de_xc_drho-rs/(3*trho)*(AU/rs+CU*(rsln+1)+DU)
          endif
          
          exc(ix,iy,iz)=e_xc*trho
          vxc(ix,iy,iz)=e_xc+trho*de_xc_drho
          
          tot_exc = tot_exc + exc(ix, iy, iz) * dv
        end do
      end do
    end do
      
    return
  End Subroutine core_exc_cor_pz

  

  
  
  
  
  
  

  subroutine core_exc_cor_pz_spin()
    implicit none
    return 
  end subroutine core_exc_cor_pz_spin

  subroutine core_exc_cor_pbe()
    implicit none
    return 
  end subroutine core_exc_cor_pbe
  
  subroutine core_exc_cor_lda()
    implicit none
    return 
  end subroutine core_exc_cor_lda
  
  subroutine core_exc_cor_tbmbj()
    implicit none
    return 
  end subroutine core_exc_cor_tbmbj
  
  subroutine core_exc_cor_tpss()
    implicit none
    return 
  end subroutine core_exc_cor_tpss
  
end subroutine

