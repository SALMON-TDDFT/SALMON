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
!This file is "prep_ps.f90"
!This file conatain one soubroutine.
!SUBROUTINE prep_ps_periodic(property)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
Subroutine prep_ps_periodic(property)
  use Global_Variables
  use salmon_parallel, only: nproc_id_global, nproc_group_tdks
  use salmon_communication, only: comm_summation, comm_is_root
  use salmon_math
  use opt_variables, only: zJxyz,zKxyz,init_for_padding, nprojector,idx_proj,idx_lma,pseudo_start_idx, init_projector
#ifdef ARTED_LBLK
  use opt_variables, only: t4ppt_nlma,t4ppt_nlma,t4ppt_i2vi,t4ppt_vi2i,t4ppt_ilma,t4ppt_j,  opt_vars_init_t4ppt
#endif
  implicit none
  character(17) :: property
  logical :: flag_alloc1, flag_alloc2
  integer :: ik,i,a,j,ix,iy,iz,lma,l,m,lm,ir,intr
  integer :: PNLx,PNLy,PNLz,narray
  real(8) :: x,y,z,r
  real(8) :: ratio1,ratio2,rc


  !(Local pseudopotential in G-space (radial part))
  if(property == 'initial') then

    allocate(lma_tbl((Lmax+1)**2,NI))

    call calc_vloc(pp,dVloc_G,Gx,Gy,Gz,NG,NG_s,NG_e,ngzero)

    !(this save is used only for opt and md options)
    save_dVloc_G(:,:)=dVloc_G(:,:)
     
  else

    dVloc_G(:,:)=save_dVloc_G(:,:)

    if(property == 'update_all') then
#define SAFE_DEALLOCATE(var) if(allocated(var)) deallocate(var)
#ifdef ARTED_LBLK
       SAFE_DEALLOCATE(t4ppt_nlma)
       SAFE_DEALLOCATE(t4ppt_i2vi)
       SAFE_DEALLOCATE(t4ppt_vi2i)
       SAFE_DEALLOCATE(t4ppt_ilma)
       SAFE_DEALLOCATE(t4ppt_j)
#endif
       SAFE_DEALLOCATE(nprojector)
       SAFE_DEALLOCATE(idx_proj)
       SAFE_DEALLOCATE(idx_lma)
       SAFE_DEALLOCATE(pseudo_start_idx)
    endif
  endif

  call calc_vpsl(pp,rhoion_G,Vpsl_ia,Vpsl,  &
                 dVloc_G,nGzero,Gx,Gy,Gz,NG,NG_s,NG_e,NL,aLxyz,Lx,Ly,Lz,Hx,Hy,Hz)
 
  !(Non-Local pseudopotential)
  if(property /= 'update_wo_realloc') then

  if (comm_is_root(nproc_id_global) .and. property=='initial') then
    write(*,*) ''
    write(*,*) '============nonlocal grid data=============='
  endif

  call calc_mps(pp,mps,nps,alx,aly,alz,lx,ly,lz,nl,hx,hy,hz)

  if (comm_is_root(nproc_id_global) .and. property == 'initial') then
     do a=1,NI
        write(*,*) 'a =',a,'Mps(a) =',Mps(a)
     end do
  endif

  !(allocate/deallocate with Nps)
  if(property == 'initial') then
     flag_alloc1=.true.
  else if(property == 'update_all') then
     narray=ubound(Jxyz,1)
     if(Nps.ne.narray)then
        deallocate(Jxyz,Jxx,Jyy,Jzz,zJxyz)
        deallocate(ekr,ekr_omp)
#ifdef ARTED_STENCIL_PADDING
        deallocate(zKxyz)
#endif
        flag_alloc1=.true.
     else
        flag_alloc1=.false.
     endif
  endif
  if(flag_alloc1)then
     allocate(Jxyz(Nps,NI),Jxx(Nps,NI),Jyy(Nps,NI),Jzz(Nps,NI),zJxyz(Nps,NI))
     allocate(ekr_omp(Nps,NI,NK_s:NK_e),ekr(Nps,NI))
#ifdef ARTED_STENCIL_PADDING
     allocate(zKxyz(Nps,NI))
#endif
  endif

  call calc_jxyz(pp,Jxyz,Jxx,Jyy,Jzz,Nps,aLx,aLy,aLz,Lx,Ly,Lz,NL,Hx,Hy,Hz)

  if(property == 'update_all') then
     zJxyz(1:Nps,1:NI) = Jxyz(1:Nps,1:NI) - 1

#ifdef ARTED_STENCIL_PADDING
     !call init_for_padding
     PNLx = NLx
     PNLy = NLy + 1
     PNLz = NLz

!$omp parallel
!$omp do private(a,j,i)
    do a=1,NI
    do j=1,Mps(a)
       i=Jxyz(j,a)
       zKxyz(j,a)=Lx(i)*PNLy*PNLz + Ly(i)*PNLz + Lz(i)
    enddo
    enddo
!$omp end do
!$omp end parallel
#endif

  endif

  lma=0
  do a=1,NI
    ik=Kion(a)
    do l=0,Mlps(ik)
      if(inorm(l,ik)==0) cycle
      do m=-l,l
        lma=lma+1
      enddo
    enddo
  enddo
  Nlma=lma


  !(allocate/deallocate with Nlma)
  if(property == 'initial') then
     flag_alloc2=.true.
  else if(property == 'update_all') then
     narray=ubound(a_tbl,1)
     if(Nlma.ne.narray .or. flag_alloc1)then
        deallocate(a_tbl,uV,duV,iuV,zproj)
        flag_alloc2=.true.
     else
        flag_alloc2=.false.
     endif
  endif
  if(flag_alloc2)then
     allocate(a_tbl(Nlma),uV(Nps,Nlma),iuV(Nlma),duV(Nps,Nlma,3))
     allocate(zproj(Nps,Nlma,NK_s:NK_e))
  endif

  lma=0
  do a=1,NI
    ik=Kion(a)
    lm=0
    do l=0,Mlps(ik)
      if(inorm(l,ik)==0) cycle
      do m=-l,l
        lm=lm+1
        lma=lma+1
        a_tbl(lma)=a
        lma_tbl(lm,a)=lma
      enddo
    enddo
  enddo

  endif  !for /= 'update_wo_realloc'


  if(property /= 'update_wo_realloc') then

  if(property == 'initial') then
     allocate( save_udVtbl_a(Nrmax,0:Lmax,NI) )
     allocate( save_udVtbl_b(Nrmax,0:Lmax,NI) )
     allocate( save_udVtbl_c(Nrmax,0:Lmax,NI) )
     allocate( save_udVtbl_d(Nrmax,0:Lmax,NI) )
  endif

  end if

  call calc_uv(pp,save_udVtbl_a,save_udVtbl_b,save_udVtbl_c,save_udVtbl_d,uV,duV, &
               Jxyz,Jxx,Jyy,Jzz,Mps,Nps,nlma,Lx,Ly,Lz,NL,Hx,Hy,Hz,aLx,aLy,aLz,  &
               lma_tbl,flag_use_grad_wf_on_force,property)

  do a=1,natom

    if(property /= 'update_wo_realloc') then
    lm=0
    do l=0,Mlps(ik)
      if(inorm(l,ik)==0) cycle
      do m=-l,l
        lm=lm+1
        iuV(lma_tbl(lm,a))=inorm(l,ik)
      enddo
    enddo
    endif

  enddo

! nonlinear core-correction
  rho_nlcc = 0d0
  tau_nlcc = 0d0
  if(flag_nlcc)then
    if(comm_is_root(nproc_id_global))write(*,"(A)")"Preparation: Non-linear core correction"
    do a=1,NI
      ik=Kion(a)
      rc = 15d0 ! maximum
      do i=1,Nrmax
        if(rho_nlcc_tbl(i,ik) + tau_nlcc_tbl(i,ik) < 1d-6)then
          rc = rad(i,ik)
          exit
        end if
        if(i == Nrmax) stop "no-cut-off"
      end do
      
      do ix=-2,2; do iy=-2,2; do iz=-2,2
        do i=1,NL
          x=Lx(i)*Hx-(Rion(1,a)+ix*aLx)
          y=Ly(i)*Hy-(Rion(2,a)+iy*aLy)
          z=Lz(i)*Hz-(Rion(3,a)+iz*aLz)
          r=sqrt(x**2+y**2+z**2)
          if(r > rc)cycle
          
          do ir=1,NRmax
            if(rad(ir,ik).gt.r) exit
          enddo
          intr=ir-1
          if (intr.lt.0.or.intr.ge.NRmax)stop 'bad intr at prep_ps'
          ratio1=(r-rad(intr,ik))/(rad(intr+1,ik)-rad(intr,ik))
          ratio2=1-ratio1
          rho_nlcc(i) = rho_nlcc(i) &
            +ratio1*rho_nlcc_tbl(intr+1,ik)+ratio2*rho_nlcc_tbl(intr,ik)
          tau_nlcc(i) = tau_nlcc(i) &
            +ratio1*tau_nlcc_tbl(intr+1,ik)+ratio2*tau_nlcc_tbl(intr,ik)
          
        enddo
        
      end do; end do; end do
    end do
  end if


  if(property == 'update_all') then
#ifdef ARTED_STENCIL_PADDING
    call init_projector(zKxyz)
#else
    call init_projector(zJxyz)
#endif
#ifdef ARTED_LBLK
     call opt_vars_init_t4ppt
#endif
  endif

  return
End Subroutine prep_ps_periodic
