subroutine core_exc_cor(xc, cval, nx, ny, nz, rho, rho_s, tau, rj, grho, rlrho, rhn_nlcc, vxc, exc)
  implicit none
  character(20), intent(in) :: xc
  real(8), intent(in) :: cval
  integer, intent(in) :: nx, ny, nz
  real(8), intent(in) :: rho(nx, ny, nz)
  real(8), intent(in) :: rho_s(nx, ny, nz, 2)
  real(8), intent(in) :: tau(nx, ny, nz)
  real(8), intent(in) :: rj(nx, ny, nz, 3)
  real(8), intent(in) :: grho(nx, ny, nz, 3)
  real(8), intent(in) :: rlrho(nx, ny, nz)
  real(8), intent(out) :: vxc(nx, ny, nz)
  real(8), intent(out) :: vxc_s(nx, ny, nz, 2)
  real(8), intent(out) :: exc
  
  select case(xc)
  case("pz")
    call core_exc_cor_pz()
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
  
  subroutine core_exc_cor_pz()
    
  end subroutine core_exc_cor_pz

  subroutine core_exc_cor_pbe()
    
  end subroutine core_exc_cor_pbe
  
  subroutine core_exc_cor_lda()
    
  end subroutine core_exc_cor_lda
  
  subroutine core_exc_cor_tbmbj()
    
  end subroutine core_exc_cor_tbmbj
  
  subroutine core_exc_cor_tpss()
    
  end subroutine core_exc_cor_tpss
  
end subroutine
