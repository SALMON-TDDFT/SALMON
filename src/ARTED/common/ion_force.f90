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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine Ion_Force_omp(Rion_update,GS_RT,ixy_m)
  use Global_Variables, only: zu_t,zu_m,zu_GS,NB,NBoccmax,calc_mode_gs,calc_mode_rt,Rion,Rion_m,force,force_m
  implicit none
  integer,intent(in) :: GS_RT
  logical,intent(in) :: Rion_update
  integer,intent(in),optional :: ixy_m

  select case(GS_RT)
    case(calc_mode_gs)
      call impl(Rion_update,zu_GS,NB)
    case(calc_mode_rt)
      if (present(ixy_m)) then
        Rion(:,:) = Rion_m(:,:,ixy_m) 
        call impl(Rion_update,zu_m(:,:,:,ixy_m),NBoccmax)
        force_m(:,:,ixy_m) = force(:,:)
      else
        call impl(Rion_update,zu_t,NBoccmax)
      end if
    case default
      call err_finalize('ion_force_omp: gs_rt flag')
  end select

contains
  subroutine impl(Rion_update,zutmp,zu_NB)
    use Global_Variables
    use salmon_parallel, only: nproc_group_tdks,get_thread_id
    use salmon_communication, only: comm_summation, comm_is_root
    use timer
    use salmon_math
    use projector
    use opt_variables, only: NUMBER_THREADS_POW2
    implicit none
    logical,intent(in)       :: Rion_update
    integer,intent(in)       :: zu_NB
    complex(8),intent(inout) :: zutmp(NL,zu_NB,NK_s:NK_e)

    integer      :: ia,ib,ilma,ik,ix,iy,iz,n,j,i
    integer      :: ix_s,ix_e,iy_s,iy_e,iz_s,iz_e
    integer,allocatable:: idx(:),idy(:),idz(:)
    real(8)      :: rab(3),rab2,Gvec(3),G2,Gd
    real(8)      :: ftmp_l(3,NI)
    real(8)      :: FionR(3,NI), FionG(3,NI),nabt_wrk(4,3)
    complex(8)   :: uVpsi,duVpsi(3)
    complex(8)   :: dzudrzu(3),dzuekrdr(3)

    complex(8), allocatable :: dzudr(:,:,:,:)

    integer :: tid
    real(8) :: ftmp_t(3,NI,0:NUMBER_THREADS_POW2-1)

    allocate(dzudr(3,NL,zu_NB,NK_s:NK_e))
    ftmp_t(:,:,:) = 0.d0


    !flag_use_grad_wf_on_force is given in Gloval_Variable
    !if .true.   use gradient of wave-function (better accuracy)
    !if .false., old way which use gradient of potential is used (less accurate)

    call timer_begin(LOG_ION_FORCE)

    !ion-ion interaction (point charge Ewald)
    if (Rion_update) then

      !ion-ion: Ewald short interaction term (real space)
!$omp parallel private(tid)
      tid = get_thread_id()
      ftmp_t(:,:,tid)=0.d0
!$omp do private(ia, ik,ix,iy,iz,ib,rab,rab2) collapse(5)
      do ia=1,NI
      do ix=-NEwald,NEwald
      do iy=-NEwald,NEwald
      do iz=-NEwald,NEwald
      do ib=1,NI
        ik=Kion(ia)
        if(ix**2+iy**2+iz**2 == 0 .and. ia == ib) cycle
        rab(1)=Rion(1,ia)-ix*aLx-Rion(1,ib)
        rab(2)=Rion(2,ia)-iy*aLy-Rion(2,ib)
        rab(3)=Rion(3,ia)-iz*aLz-Rion(3,ib)
        rab2=sum(rab(:)**2)
        ftmp_t(:,ia,tid)=ftmp_t(:,ia,tid)&
             &-Zps(Kion(ia))*Zps(Kion(ib))*rab(:)/sqrt(rab2)*(-erfc_salmon(sqrt(aEwald*rab2))/rab2&
             &-2*sqrt(aEwald/(rab2*Pi))*exp(-aEwald*rab2))
      end do
      end do
      end do
      end do
      end do
!$omp end do
      call mreduce_omp(ftmp_l,ftmp_t,3*NI,tid)
!$omp end parallel
      FionR = ftmp_l


      !ion-ion: Ewald long interaction term (wave-number space)
!$omp parallel private(ia,tid)
      tid = get_thread_id()
      ftmp_t(:,:,tid)=0.d0
      do ia=1,NI
!$omp do private(ik,n,Gvec,G2,Gd)
      do n=NG_s,NG_e
        if(n == nGzero) cycle
        ik=Kion(ia)
        Gvec(1)=Gx(n); Gvec(2)=Gy(n); Gvec(3)=Gz(n)
        G2=sum(Gvec(:)**2)
        Gd=sum(Gvec(:)*Rion(:,ia))
        ftmp_t(:,ia,tid) = ftmp_t(:,ia,tid) &
        &   + Gvec(:)*(4*Pi/G2)*exp(-G2/(4*aEwald))*Zps(ik) &
        &     *zI*0.5d0*(conjg(rhoion_G(n))*exp(-zI*Gd)-rhoion_G(n)*exp(zI*Gd))
      end do
!$omp end do
      end do
      call mreduce_omp(ftmp_l,ftmp_t,3*NI,tid)
!$omp end parallel

      call comm_summation(ftmp_l,FionG,3*NI,nproc_group_tdks)
      Fion = FionR + FionG

    end if


    call update_projector_omp(kac)


    ! ion-electron 
    if(flag_use_grad_wf_on_force)then

    ! Use gradient of wave-func for calculating force on ions

    !(prepare gradient of w.f.)
    ix_s = 0  ;  ix_e = NLx-1
    iy_s = 0  ;  iy_e = NLy-1
    iz_s = 0  ;  iz_e = NLz-1

    allocate(idx(ix_s-4:ix_e+4),idy(iy_s-4:iy_e+4),idz(iz_s-4:iz_e+4))
    do i=ix_s-4,ix_e+4
      idx(i) = mod(NLx+i,NLx)
    end do
    do i=iy_s-4,iy_e+4
      idy(i) = mod(NLy+i,NLy)
    end do
    do i=iz_s-4,iz_e+4
      idz(i) = mod(NLz+i,NLz)
    end do

    nabt_wrk(1:4,1) = nabx(1:4)
    nabt_wrk(1:4,2) = naby(1:4)
    nabt_wrk(1:4,3) = nabz(1:4)

!$omp parallel do collapse(2) default(none) &
!$omp          private(ik,ib) &
!$omp          shared(zutmp,dzudr,idx,idy,idz,nabt_wrk) &
!$omp          firstprivate(ix_s,ix_e,iy_s,iy_e,iz_s,iz_e,nk_s,nk_e,nboccmax)
    do ik=NK_s,NK_e
    do ib=1,NBoccmax
       call stencil_C_zu(zutmp(:,ib,ik) &
       &                ,dzudr(:,:,ib,ik) &
       &                ,ix_s,ix_e,iy_s,iy_e,iz_s,iz_e &
       &                ,idx,idy,idz,nabt_wrk)
    enddo
    enddo
!$omp end parallel do

    !Force from Vlocal with wave-func gradient --
!$omp parallel private(ia,tid)
    tid = get_thread_id()
    ftmp_t(:,:,tid) = 0.d0
!$omp do private(ik,ib,i,ia,dzudrzu) collapse(3)
    do ik=NK_s,NK_e
    do ib=1,NBoccmax
    do i=1,NL

       dzudrzu(:)=conjg(dzudr(:,i,ib,ik))*zutmp(i,ib,ik)

       do ia=1,NI
          ftmp_t(:,ia,tid) = ftmp_t(:,ia,tid) &
          &  -2d0* dble(dzudrzu(:)*Vpsl_ia(i,ia))*occ(ib,ik)*Hxyz
       enddo

    enddo
    enddo
    enddo
!$omp end do
    call mreduce_omp(ftmp_l,ftmp_t,3*NI,tid)
!$omp end parallel

    call timer_begin(LOG_ALLREDUCE)
    call comm_summation(ftmp_l,Floc,3*NI,nproc_group_tdks)
    call timer_end(LOG_ALLREDUCE)


    !Non-Local pseudopotential term using gradient of w.f.
!$omp parallel private(ia,tid)
    tid = get_thread_id()
    ftmp_t(:,:,tid) = 0.d0
!$omp do private(ik,j,i,ib,ilma,uVpsi,duVpsi,dzuekrdr) collapse(2)
    do ik=NK_s,NK_e
    do ib=1,NBoccmax
       do ilma=1,Nlma
          ia=a_tbl(ilma)
           uVpsi   =0.d0
          duVpsi(:)=0.d0
          do j=1,Mps(ia)
             i=Jxyz(j,ia)
            uVpsi      =uVpsi + uV(j,ilma)*ekr_omp(j,ia,ik)*zutmp(i,ib,ik)
            dzuekrdr(:)=(dzudr(:,i,ib,ik)+zI*kAc(ik,:)*zutmp(i,ib,ik))*ekr_omp(j,ia,ik)
            duVpsi(:)  =duVpsi(:) + conjg(dzuekrdr(:))*uV(j,ilma)
          end do
          uVpsi    =uVpsi    *Hxyz
          duVpsi(:)=duVpsi(:)*Hxyz
          ftmp_t(:,ia,tid) = ftmp_t(:,ia,tid) &
          &                  -2d0* dble(uVpsi*duVpsi(:))*iuV(ilma)*occ(ib,ik)
       end do
    end do
    end do
!$omp end do
    call mreduce_omp(ftmp_l,ftmp_t,3*NI,tid)
!$omp end parallel

    call timer_begin(LOG_ALLREDUCE)
    call comm_summation(ftmp_l,Fnl,3*NI,nproc_group_tdks)
    call timer_end(LOG_ALLREDUCE)


    else

    ! Use gradient of potential for calculating force on ions
    ! (older method, less accurate for larger grid size)

!$omp parallel private(ia,tid)
    tid = get_thread_id()
    ftmp_t(:,:,tid) = 0.d0
    do ia=1,NI
!$omp do private(ik,n,Gvec,G2,Gd)
    do n=NG_s,NG_e
      if(n == nGzero) cycle
      ik=Kion(ia)
      Gvec(1)=Gx(n); Gvec(2)=Gy(n); Gvec(3)=Gz(n)
      G2=sum(Gvec(:)**2)
      Gd=sum(Gvec(:)*Rion(:,ia))
      ftmp_t(:,ia,tid) = ftmp_t(:,ia,tid) &
      &    + zI*Gvec(:)*(4*Pi/G2)*Zps(ik)*exp(zI*Gd)*rhoe_G(n) &
      &    + conjg(rhoe_G(n))*dVloc_G(n,ik)*zI*Gvec(:)*exp(-zI*Gd)
    end do
!$omp end do
    end do
    call mreduce_omp(ftmp_l,ftmp_t,3*NI,tid)
!$omp end parallel

    call comm_summation(ftmp_l,Floc,3*NI,nproc_group_tdks)

    !non-local pseudopotential term
!$omp parallel private(ia,tid)
    tid = get_thread_id()
    ftmp_t(:,:,tid) = 0.d0
!$omp do private(ik,j,i,ib,ilma,uVpsi,duVpsi) collapse(2)
    do ik=NK_s,NK_e
    do ib=1,NBoccmax
      do ilma=1,Nlma
        ia=a_tbl(ilma)
        uVpsi=0.d0
        duVpsi(:)=0.d0
        do j=1,Mps(ia)
          i=Jxyz(j,ia)
          uVpsi=uVpsi+uV(j,ilma)*ekr_omp(j,ia,ik)*zutmp(i,ib,ik)
          duVpsi(:)=duVpsi(:)+duV(j,ilma,:)*ekr_omp(j,ia,ik)*zutmp(i,ib,ik)
        end do
        uVpsi=uVpsi*Hxyz; duVpsi(:)=duVpsi(:)*Hxyz
        ftmp_t(:,ia,tid)=ftmp_t(:,ia,tid)+(conjg(uVpsi)*duVpsi(:)+uVpsi*conjg(duVpsi(:)))*iuV(ilma)*occ(ib,ik)
      end do
    end do
    end do
!$omp end do
    call mreduce_omp(ftmp_l,ftmp_t,3*NI,tid)
!$omp end parallel

    call timer_begin(LOG_ALLREDUCE)
    call comm_summation(ftmp_l,fnl,3*NI,nproc_group_tdks)
    call timer_end(LOG_ALLREDUCE)

    endif ! flag_use_grad_wf_on_force


    call timer_begin(LOG_ALLREDUCE)
    force = Fion + Floc + Fnl
    call timer_end(LOG_ALLREDUCE)

    call timer_end(LOG_ION_FORCE)
  end subroutine

!OpenMP manual reduction routine
  subroutine mreduce_omp(vout,vtmp,vsize,tid)
    use opt_variables, only: NUMBER_THREADS_POW2
    implicit none
    integer,intent(in)  :: vsize,tid
    real(8),intent(out) :: vout(vsize)
    real(8)             :: vtmp(vsize,0:NUMBER_THREADS_POW2-1)

    integer :: i

    i = NUMBER_THREADS_POW2/2
    do while(i > 0)
      if(tid < i) then
        vtmp(:,tid) = vtmp(:,tid) + vtmp(:,tid + i)
      end if
      i = i/2
!$omp barrier
    end do

    if (tid == 0) then
      vout(:) = vtmp(:,0)
    end if
  end subroutine

!Gradient of wave function (du/dr) with nine points formura
# define DX(dt) iz,iy,idx(ix+(dt))
# define DY(dt) iz,idy(iy+(dt)),ix
# define DZ(dt) idz(iz+(dt)),iy,ix
  subroutine stencil_C_zu(zu0,dzu0dr &
  &                      ,ix_s,ix_e,iy_s,iy_e,iz_s,iz_e &
  &                      ,idx,idy,idz,nabt)
  implicit none
  integer   ,intent(in)  :: ix_s,ix_e,iy_s,iy_e,iz_s,iz_e
  integer   ,intent(in)  :: idx(ix_s-4:ix_e+4),idy(iy_s-4:iy_e+4),idz(iz_s-4:iz_e+4)
  real(8)   ,intent(in)  :: nabt(4,3)
  complex(8),intent(in)  :: zu0(iz_s:iz_e,iy_s:iy_e,ix_s:ix_e)
  complex(8),intent(out) :: dzu0dr(3,iz_s:iz_e,iy_s:iy_e,ix_s:ix_e)
  !
  integer :: iz,iy,ix
  complex(8) :: w(3)

  do ix=ix_s,ix_e
  do iy=iy_s,iy_e
!dir$ vector nontemporal(dzu0dr)
  do iz=iz_s,iz_e

    w(1) =  nabt(1,1)*(zu0(DX(1)) - zu0(DX(-1))) &
           +nabt(2,1)*(zu0(DX(2)) - zu0(DX(-2))) &
           +nabt(3,1)*(zu0(DX(3)) - zu0(DX(-3))) &
           +nabt(4,1)*(zu0(DX(4)) - zu0(DX(-4)))

    w(2) =  nabt(1,2)*(zu0(DY(1)) - zu0(DY(-1))) &
           +nabt(2,2)*(zu0(DY(2)) - zu0(DY(-2))) &
           +nabt(3,2)*(zu0(DY(3)) - zu0(DY(-3))) &
           +nabt(4,2)*(zu0(DY(4)) - zu0(DY(-4)))

    w(3) =  nabt(1,3)*(zu0(DZ(1)) - zu0(DZ(-1))) &
           +nabt(2,3)*(zu0(DZ(2)) - zu0(DZ(-2))) &
           +nabt(3,3)*(zu0(DZ(3)) - zu0(DZ(-3))) &
           +nabt(4,3)*(zu0(DZ(4)) - zu0(DZ(-4)))

    dzu0dr(:,iz,iy,ix) = w(:)
  end do
  end do
  end do

  return
  end subroutine stencil_C_zu


end subroutine Ion_Force_omp
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
