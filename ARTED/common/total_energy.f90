!
!  Copyright 2016 ARTED developers
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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
subroutine Total_Energy_omp(Rion_update,GS_RT,ixy_m)
  use Global_Variables, only: zu_t,zu_m,zu_GS,NB,NBoccmax,calc_mode_gs,calc_mode_rt
  use communication
  implicit none
  integer,intent(in) :: GS_RT
  logical,intent(in) :: Rion_update
  integer,intent(in),optional :: ixy_m

  select case(GS_RT)
    case(calc_mode_gs)
      call impl(Rion_update,zu_GS,NB)
    case(calc_mode_rt)
      if (present(ixy_m)) then
        call impl(Rion_update,zu_m(:,:,:,ixy_m),NBoccmax)
      else
        call impl(Rion_update,zu_t,NBoccmax)
      end if
    case default
      call err_finalize('total_energy_omp: GS_RT flag')
  end select

contains
  subroutine impl(Rion_update,zutmp,zu_NB)
    use Global_Variables
    use Opt_Variables
    use timer
    implicit none
    logical,intent(in)       :: Rion_update
    integer,intent(in)       :: zu_NB
    complex(8),intent(inout) :: zutmp(0:NL-1,zu_NB,NK_s:NK_e)

    integer      :: ik,ib,ia,ix,iy,iz,n,ilma,j,i
    real(8)      :: rab(3),rab2,G2,Gd,kr
    complex(8)   :: uVpsi
    real(8)      :: Ekin_l,Enl_l,Eh_l,Eion_l,sum_tmp(5),sum_result(5)
    real(8)      :: Eion_tmp1,Eion_tmp2,Eloc_l1,Eloc_l2
    complex(8)   :: tpsum
    real(8)      :: nabt(12)
    real(8)      :: lap0_2
    integer      :: thr_id,omp_get_thread_num

#if defined(__KNC__) || defined(__AVX512F__)
# define MEM_ALIGNED 64
#else
# define MEM_ALIGNED 32
#endif
!dir$ attributes align:MEM_ALIGNED :: nabt

    call timer_begin(LOG_TOTAL_ENERGY)

    !ion-ion
    if (Rion_update) then
      thr_id=0
      Eion_tmp1=0.d0
      Eion_l=0.d0
!$omp parallel private(ia,ix,iy,iz,ib,rab,rab2,n,G2,thr_id)
!$    thr_id=omp_get_thread_num()

!$omp do reduction(+:Eion_tmp1) collapse(5)
      do ia=1,NI
      do ix=-NEwald,NEwald
      do iy=-NEwald,NEwald
      do iz=-NEwald,NEwald
      do ib=1,NI
        if (ix**2+iy**2+iz**2 == 0 .and. ia == ib) then
          cycle
        end if
        rab(1)=Rion(1,ia)-ix*aLx-Rion(1,ib)
        rab(2)=Rion(2,ia)-iy*aLy-Rion(2,ib)
        rab(3)=Rion(3,ia)-iz*aLz-Rion(3,ib)
        rab2=sum(rab(:)**2)
        Eion_tmp1=Eion_tmp1 + 0.5d0*Zps(Kion(ia))*Zps(Kion(ib))*erfc(sqrt(aEwald*rab2))/sqrt(rab2)
      end do
      end do
      end do
      end do
      end do
!$omp end do

!$omp do reduction(+:Eion_l)
      do n=NG_s,NG_e
        if(n == nGzero) cycle
        G2=Gx(n)**2+Gy(n)**2+Gz(n)**2
        Eion_l=Eion_l+aLxyz*(4*Pi/G2)*(abs(rhoion_G(n))**2*exp(-G2/(4*aEwald))*0.5d0)
      end do
!$omp end do
!$omp end parallel

      Eion_tmp1=Eion_tmp1-Pi*sum(Zps(Kion(:)))**2/(2*aEwald*aLxyz)-sqrt(aEwald/Pi)*sum(Zps(Kion(:))**2)
    end if

    lap0_2=-(lapx(0)+lapy(0)+lapz(0))*0.5d0

    Eall=0.d0
    Ekin_l=0.d0
    Eh_l=0.d0
    Eloc_l1=0.d0
    Eloc_l2=0.d0
    Enl_l=0.d0

!$omp parallel private(thr_id)
!$  thr_id=omp_get_thread_num()

!$omp do collapse(2) private(ia,n,Gd) reduction(+:Eloc_l2)
    do ia=1,NI
    do n=NG_s,NG_e
      Gd=Gx(n)*Rion(1,ia)+Gy(n)*Rion(2,ia)+Gz(n)*Rion(3,ia)
      Eloc_l2=Eloc_l2+conjg(rhoe_G(n))*dVloc_G(n,Kion(ia))*exp(-zI*Gd)
    end do
    end do
!$omp end do

!$omp do private(n,G2) reduction(+:Eh_l,Eloc_l1)
    do n=NG_s,NG_e
      if(n == nGzero) cycle
      G2=Gx(n)**2+Gy(n)**2+Gz(n)**2
      Eh_l=Eh_l+aLxyz*(4*Pi/G2)*(abs(rhoe_G(n))**2*0.5d0)
      Eloc_l1=Eloc_l1+aLxyz*(4*Pi/G2)*(-rhoe_G(n)*conjg(rhoion_G(n)))
    end do
!$omp end do

!$omp do private(ia,j,i,ix,iy,iz,kr) collapse(2)
    do ik=NK_s,NK_e
    do ia=1,NI
    do j=1,Mps(ia)
      i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
      kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
      ekr_omp(j,ia,ik)=exp(zI*kr)
    end do
    end do
    end do
!$omp end do

!$omp do private(ik,ib,nabt,tpsum,i,j,ilma,uVpsi,ia) &
!$omp   &reduction(+:Ekin_l,Enl_l) &
!$omp   &collapse(2)
    do ik=NK_s,NK_e
    do ib=1,NBoccmax

      nabt( 1: 4)=kAc(ik,1)*nabx(1:4)
      nabt( 5: 8)=kAc(ik,2)*naby(1:4)
      nabt( 9:12)=kAc(ik,3)*nabz(1:4)

      call total_energy_stencil(lap0_2,lapt,nabt,zutmp(:,ib,ik),tpsum);
      Ekin_l=Ekin_l+occ(ib,ik)*tpsum*Hxyz+occ(ib,ik)*sum(kAc(ik,:)**2)/2.d0

!dir$ vector aligned
      do ilma=1,Nlma
        ia=a_tbl(ilma)
        uVpsi=0.d0
        do j=1,Mps(ia)
          uVpsi=uVpsi+uV(j,ilma)*ekr_omp(j,ia,ik)*zutmp(zJxyz(j,ia),ib,ik)
        enddo
        uVpsi=uVpsi*Hxyz
        Enl_l=Enl_l+occ(ib,ik)*abs(uVpsi)**2*iuV(ilma)
      end do
    end do
    end do
!$omp end do nowait
!$omp end parallel
    call timer_end(LOG_TOTAL_ENERGY)

    call timer_begin(LOG_ALLREDUCE)
    !summarize
    if (Rion_update) then
      call comm_summation(Eion_l,Eion_tmp2,proc_group(2))
      Eion=Eion_tmp1+Eion_tmp2
    end if

    sum_tmp(1) = Ekin_l
    sum_tmp(2) = Eh_l
    sum_tmp(3) = Enl_l
    sum_tmp(4) = Eloc_l1
    sum_tmp(5) = Eloc_l2
    call comm_summation(sum_tmp,sum_result,5,proc_group(2))
    Ekin = sum_result(1)
    Eh   = sum_result(2)
    Enl  = sum_result(3)
    Eloc = sum_result(4) + sum_result(5)

    !Exchange correlation
    Exc=sum(Eexc)*Hxyz

    Eall=Ekin+Eloc+Enl+Exc+Eh+Eion
    call timer_end(LOG_ALLREDUCE)

  end subroutine
end subroutine Total_Energy_omp
