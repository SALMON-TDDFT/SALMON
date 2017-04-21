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
subroutine prep_backup_values(is_backup)
  use global_variables
  use timer,           only: timer_reentrance_read, timer_reentrance_write
  use opt_variables,   only: opt_vars_initialize_p1, opt_vars_initialize_p2
  use backup_routines, only: backup_value
  use misc_routines,   only: get_wtime
  use communication
  implicit none
  logical, intent(in) :: is_backup

  integer, parameter :: iounit = 500
  character(256)     :: reent_filename, dump_filename
  real(8) :: beg_time, end_time
  integer :: gNt


  call comm_sync_all; beg_time = get_wtime()

  if (is_backup) then
    ! backup
    if (comm_is_root()) then
      ! create entry input file.
      write(reent_filename,'(A,A,I6.6,A)') trim(SYSname),'_re_',iter_now,'.dat'
      write(dump_filename,'(A,A,A,I6.6)') 'dump_',trim(SYSname),'_',iter_now
      open(iounit, file=gen_filename(reent_filename))
      write(iounit,"(A)") '&group_function'
      write(iounit,"(A,A,A)") "cfunction='",trim(calc_mode),"'"
      write(iounit,"(A)") '/'
      write(iounit,"(A)") '&control'
      write(iounit,"(A)") "entrance_option = 'reentrance'  "
      write(iounit,"(A,e26.16e3)") "Time_shutdown=",Time_shutdown
      write(iounit,"(A)") '/'
      write(iounit,"(A)") '&reentrance'
      write(iounit,"(A)") "'"//trim(directory)//"'  !directory"
      write(iounit,"(A)") "'"//trim(dump_filename)//"'  !dump_filename"
      if(calc_mode == calc_mode_ms)then
         write(iounit,'(I7,A)')Nt,'  ! Nt: If you want continuous execution, please change the value.'
      end if
      write(iounit,"(A)") '/'
      close(iounit)
    end if
    call comm_bcast(dump_filename, proc_group(1))
    open(iounit, status='replace', form='unformatted', file=gen_filename(dump_filename, procid(1)))
  else
    ! restore
    if(comm_is_root()) then
      open(iounit, file='.reenetrance.tmp', status='old')
      read(iounit, *)directory
      read(iounit, *)dump_filename
      if(calc_mode == calc_mode_ms) read(iounit, *)gNt
      close(iounit)

    end if
    call comm_bcast(directory, proc_group(1))
    call comm_bcast(dump_filename, proc_group(1))
    if(calc_mode == calc_mode_ms) call comm_bcast(gNt, proc_group(1))

    write (process_directory,'(A,A,I5.5,A)') trim(directory),'/work_p',procid(1),'/'
    open(iounit, status='old', form='unformatted', file=gen_filename(dump_filename, procid(1)))
  end if

!============= backup values ==============
#define BACKUP(TARGET_VAR) call backup_value(is_backup, iounit, TARGET_VAR)

  BACKUP(iter_now)
  BACKUP(entrance_iter)

! grid
  BACKUP(NLx)
  BACKUP(NLy)
  BACKUP(NLz)
  BACKUP(Nd)
  BACKUP(NL)
  BACKUP(NG)
  BACKUP(NKx)
  BACKUP(NKy)
  BACKUP(NKz)
  BACKUP(NK)
  BACKUP(Sym)
  BACKUP(nGzero)
  BACKUP(NKxyz)
  BACKUP(aL)
  BACKUP(ax)
  BACKUP(ay)
  BACKUP(az)
  BACKUP(aLx)
  BACKUP(aLy)
  BACKUP(aLz)
  BACKUP(aLxyz)
  BACKUP(bLx)
  BACKUP(bLy)
  BACKUP(bLz)
  BACKUP(Hx)
  BACKUP(Hy)
  BACKUP(Hz)
  BACKUP(Hxyz)
  BACKUP(Lx)
  BACKUP(Ly)
  BACKUP(Lz)
  BACKUP(Lxyz)
  BACKUP(ifdx)
  BACKUP(ifdy)
  BACKUP(ifdz)
  BACKUP(Gx)
  BACKUP(Gy)
  BACKUP(Gz)
  BACKUP(lap)
  BACKUP(nab)
  BACKUP(lapx)
  BACKUP(lapy)
  BACKUP(lapz)
  BACKUP(nabx)
  BACKUP(naby)
  BACKUP(nabz)

! pseudopotential
  BACKUP(ps_type)
  BACKUP(ps_format)
  BACKUP(PSmask_option)
  BACKUP(alpha_mask)
  BACKUP(gamma_mask)
  BACKUP(eta_mask)
  BACKUP(Nps)
  BACKUP(Nlma)
  BACKUP(Mps)
  BACKUP(Jxyz)
  BACKUP(Jxx)
  BACKUP(Jyy)
  BACKUP(Jzz)
  BACKUP(Mlps)
  BACKUP(Lref)
  BACKUP(Zps)
  BACKUP(NRloc)
  BACKUP(NRps)
  BACKUP(inorm)
  BACKUP(iuV)
  BACKUP(a_tbl)
  BACKUP(rad)
  BACKUP(Rps)
  BACKUP(vloctbl)
  BACKUP(udVtbl)
  BACKUP(radnl)
  BACKUP(Rloc)
  BACKUP(uV)
  BACKUP(duV)
  BACKUP(anorm)
  BACKUP(dvloctbl)
  BACKUP(dudVtbl)

! material
  BACKUP(NI)
  BACKUP(NE)
  BACKUP(NB)
  BACKUP(NBoccmax)
  BACKUP(Ne_tot)
  BACKUP(Zatom)
  BACKUP(Kion)
  BACKUP(Rion)
  BACKUP(Mass)
  BACKUP(Rion_eq)
  BACKUP(dRion)
  BACKUP(occ)
  BACKUP(wk)

! physical quantities
  BACKUP(Eall)
  BACKUP(Eall0)
  BACKUP(jav(3))
  BACKUP(Tion)
  BACKUP(Ekin)
  BACKUP(Eloc)
  BACKUP(Enl)
  BACKUP(Eh)
  BACKUP(Exc)
  BACKUP(Eion)
  BACKUP(Eelemag)
  BACKUP(javt)
  BACKUP(Vpsl)
  BACKUP(Vh)
  BACKUP(Vexc)
  BACKUP(Eexc)
  BACKUP(Vloc)
  BACKUP(Vloc_GS)
  BACKUP(Vloc_t)
  BACKUP(Vloc_new)
  BACKUP(Vloc_old)
  BACKUP(tmass)
  BACKUP(tjr)
  BACKUP(tjr2)
  BACKUP(tmass_t)
  BACKUP(tjr_t)
  BACKUP(tjr2_t)
  BACKUP(dVloc_G)
  BACKUP(rho)
  BACKUP(rho_gs)
  BACKUP(rhoe_G)
  BACKUP(rhoion_G)
  BACKUP(force)
  BACKUP(esp)
  BACKUP(force_ion)
  BACKUP(Floc)
  BACKUP(Fnl)
  BACKUP(Fion)
  BACKUP(ovlp_occ_l)
  BACKUP(ovlp_occ)
  BACKUP(Nelec)
  BACKUP(NBocc)
  BACKUP(esp_vb_min)
  BACKUP(esp_vb_max)
  BACKUP(esp_cb_min)
  BACKUP(esp_cb_max)

! Nonlinear core correction
  BACKUP(flag_nlcc)
  BACKUP(rho_nlcc)
  BACKUP(tau_nlcc)

! wave functions, work array
  BACKUP(zu_t)
  BACKUP(zu_GS)
  BACKUP(zu_GS0)
  BACKUP(esp_var)

! variables for 4-times loop in Fourier transportation
  BACKUP(nxyz)
  BACKUP(rho_3D)
  BACKUP(Vh_3D)
  BACKUP(rhoe_G_temp)
  BACKUP(rhoe_G_3D)
  BACKUP(f1)
  BACKUP(f2)
  BACKUP(f3)
  BACKUP(f4)
  BACKUP(eGx)
  BACKUP(eGy)
  BACKUP(eGz)
  BACKUP(eGxc)
  BACKUP(eGyc)
  BACKUP(eGzc)

! Bloch momentum,laser pulse, electric field
  BACKUP(AE_shape)
  BACKUP(f0_1)
  BACKUP(IWcm2_1)
  BACKUP(tpulsefs_1)
  BACKUP(omegaev_1)
  BACKUP(omega_1)
  BACKUP(tpulse_1)
  BACKUP(Epdir_1(3))
  BACKUP(phi_CEP_1)
  BACKUP(f0_2)
  BACKUP(IWcm2_2)
  BACKUP(tpulsefs_2)
  BACKUP(omegaev_2)
  BACKUP(omega_2)
  BACKUP(tpulse_2)
  BACKUP(Epdir_2(3))
  BACKUP(phi_CEP_2)
  BACKUP(T1_T2fs)
  BACKUP(T1_T2)
  BACKUP(E_ext)
  BACKUP(E_ind)
  BACKUP(E_tot)
  BACKUP(kAc)
  BACKUP(kAc0)
  BACKUP(kAc_new)
  BACKUP(Ac_ext)
  BACKUP(Ac_ind)
  BACKUP(Ac_tot)

! control parameters
  BACKUP(NEwald)
  BACKUP(aEwald)
  BACKUP(Ncg)
  BACKUP(dt)
  BACKUP(dAc)
  BACKUP(domega)
  BACKUP(Nscf)
  BACKUP(Nt)
  BACKUP(Nomega)
  BACKUP(Nmemory_MB)
  BACKUP(alpha_MB)
  BACKUP(NFSset_start)
  BACKUP(NFSset_every)

! file names, flags, etc
  BACKUP(SYSname)
  BACKUP(directory)
  BACKUP(file_GS)
  BACKUP(file_RT)
  BACKUP(file_epst)
  BACKUP(file_epse)
  BACKUP(file_force_dR)
  BACKUP(file_j_ac)
  BACKUP(file_DoS)
  BACKUP(file_band)
  BACKUP(file_dns)
  BACKUP(file_ovlp)
  BACKUP(file_nex)
  BACKUP(file_kw)
  BACKUP(file_energy_transfer)
  BACKUP(file_ac_vac)
  BACKUP(file_ac_vac_back)
  BACKUP(file_ac_m)
  BACKUP(file_ac)

  BACKUP(ext_field)
  BACKUP(Longi_Trans)
  BACKUP(FSset_option)
  BACKUP(MD_option)
  BACKUP(AD_RHO)
  BACKUP(functional)
  BACKUP(cval)
  BACKUP(propagator)

  BACKUP(NK_ave)
  BACKUP(NG_ave)
  BACKUP(NK_s)
  BACKUP(NK_e)
  BACKUP(NG_s)
  BACKUP(NG_e)
  BACKUP(NK_remainder)
  BACKUP(NG_remainder)

  BACKUP(position_option)

  BACKUP(ekr)

  BACKUP(ekr_omp)
  BACKUP(tpsi_omp)
  BACKUP(ttpsi_omp)
  BACKUP(htpsi_omp)
  BACKUP(xk_omp)
  BACKUP(hxk_omp)
  BACKUP(gk_omp)
  BACKUP(pk_omp)
  BACKUP(pko_omp)
  BACKUP(txk_omp)
  BACKUP(NKB)
  BACKUP(ik_table)
  BACKUP(ib_table)

  BACKUP(tau_s_l_omp)
  BACKUP(j_s_l_omp)

! sym
  BACKUP(crystal_structure)
  BACKUP(itable_sym)
  BACKUP(rho_l)
  BACKUP(rho_tmp1)
  BACKUP(rho_tmp2)

! Finite temperature
  BACKUP(KbTev)

! multi scale
  BACKUP(FDTDdim)
  BACKUP(TwoD_shape)
  BACKUP(NXY_s)
  BACKUP(NXY_e)
  BACKUP(NKsplit)
  BACKUP(NXYsplit)
  BACKUP(macRANK)
  BACKUP(kRANK)

  BACKUP(HX_m)
  BACKUP(HY_m)
  BACKUP(NX_m)
  BACKUP(NXvacL_m)
  BACKUP(NXvacR_m)
  BACKUP(NY_m)
  BACKUP(NYvacT_m)
  BACKUP(NYvacB_m)
  BACKUP(Ac_m)
  BACKUP(Ac_new_m)
  BACKUP(Ac_old_m)
  BACKUP(Elec)
  BACKUP(Bmag)
  BACKUP(j_m)
  BACKUP(jmatter_m)
  BACKUP(jmatter_m_l)
  BACKUP(g)
  BACKUP(bcon)
  BACKUP(NX_table)
  BACKUP(NY_table)
  BACKUP(BC_my)

  BACKUP(zu_m)
  BACKUP(Vh_m)
  BACKUP(Vexc_m)
  BACKUP(Eexc_m)
  BACKUP(Vloc_m)
  BACKUP(Vloc_old_m)
  BACKUP(rho_m)
  BACKUP(energy_joule)
  BACKUP(energy_elec_Matter_l)
  BACKUP(energy_elec_Matter)
  BACKUP(energy_elec)
  BACKUP(energy_elemag)
  BACKUP(energy_total)
  BACKUP(excited_electron_l)
  BACKUP(excited_electron)

  BACKUP(data_out)
  BACKUP(data_local_Ac)
  BACKUP(data_local_jm)
  BACKUP(data_vac_Ac)
  BACKUP(Nstep_write)


  BACKUP(backup_frequency)
  BACKUP(need_backup)

  if (is_backup) then
    call timer_reentrance_write(iounit)
  else
    call timer_reentrance_read(iounit)
  end if
!============= backup values ==============
  close(iounit)

! initialize
  if (.not. is_backup) then
    call comm_set_level2_group(macRANK, kRANK)
    call opt_vars_initialize_p1
    call opt_vars_initialize_p2

    if(calc_mode == calc_mode_ms)then
  ! continuous execution (is available only multi-scale mode)
      if (gNt /= Nt) then
        if (comm_is_root()) then
          print '(A,I7,A,I7)', '*** [Nt is updated] start continuous execution:',Nt,' to ',gNt
        end if
        Nt                 = gNt
        Ndata_out          = Nt / Nstep_write
        Ndata_out_per_proc = Ndata_out / nprocs(1)
        call resize_arrays
      end if
    end if
  end if

  call comm_sync_all; end_time = get_wtime()

  if (comm_is_root()) then
    if (is_backup) then
      write(*,*) 'Backup time =',end_time-beg_time,' sec'
    else
      write(*,*) 'Restore time =',end_time-beg_time,' sec'
    end if
  end if

contains
  function gen_filename(base, procid)
    use global_variables, only: directory, process_directory
    implicit none
    character(*), intent(in)      :: base
    integer, intent(in), optional :: procid
    character(128)                :: gen_filename
    character(32) :: cMyrank
    if (present(procid)) then
      write(cMyrank,'(I5.5)') procid
      gen_filename = trim(process_directory)//trim(base)//'.p'//trim(cMyrank)
    else
      gen_filename = trim(directory)//trim(base)
    end if
  end function

  ! TODO: An array resizing subroutine should be provided.
  subroutine resize_arrays
    implicit none
    real(8), allocatable :: tmp2(:,:),tmp3(:,:,:),tmp4(:,:,:,:)
    integer :: mt

    ! javt
    mt = min(Nt, ubound(javt, 1))
    allocate(tmp2(0:Nt,3))
    tmp2(:,:) = 0.d0
    tmp2(0:mt,:) = javt(0:mt,:)
    deallocate(javt)
    allocate(javt(0:Nt,3))
    javt(:,:) = tmp2(:,:)
    deallocate(tmp2)

    ! E_ext
    mt = min(Nt, ubound(E_ext, 1))
    allocate(tmp2(0:Nt,3))
    tmp2(:,:) = 0.d0
    tmp2(0:mt,:) = E_ext(0:mt,:)
    deallocate(E_ext)
    allocate(E_ext(0:Nt,3))
    E_ext(:,:) = tmp2(:,:)
    deallocate(tmp2)

    ! E_ind
    mt = min(Nt, ubound(E_ind, 1))
    allocate(tmp2(0:Nt,3))
    tmp2(:,:) = 0.d0
    tmp2(0:mt,:) = E_ind(0:mt,:)
    deallocate(E_ind)
    allocate(E_ind(0:Nt,3))
    E_ind(:,:) = tmp2(:,:)
    deallocate(tmp2)

    ! E_tot
    mt = min(Nt, ubound(E_tot, 1))
    allocate(tmp2(0:Nt,3))
    tmp2(:,:) = 0.d0
    tmp2(0:mt,:) = E_tot(0:mt,:)
    deallocate(E_tot)
    allocate(E_tot(0:Nt,3))
    E_tot(:,:) = tmp2(:,:)
    deallocate(tmp2)

    ! Ac_ext
    mt = min(Nt+1, ubound(Ac_ind, 1))
    allocate(tmp2(-1:Nt+1,3))
    tmp2(:,:) = 0.d0
    tmp2(-1:mt,:) = Ac_ext(-1:mt,:)
    deallocate(Ac_ext)
    allocate(Ac_ext(-1:Nt+1,3))
    Ac_ext(:,:) = tmp2(:,:)
    deallocate(tmp2)

    ! Ac_ind
    mt = min(Nt+1, ubound(Ac_ind, 1))
    allocate(tmp2(-1:Nt+1,3))
    tmp2(:,:) = 0.d0
    tmp2(-1:mt,:) = Ac_ind(-1:mt,:)
    deallocate(Ac_ind)
    allocate(Ac_ind(-1:Nt+1,3))
    Ac_ind(:,:) = tmp2(:,:)
    deallocate(tmp2)

    ! Ac_tot
    mt = min(Nt+1, ubound(Ac_tot, 1))
    allocate(tmp2(-1:Nt+1,3))
    tmp2(:,:) = 0.d0
    tmp2(-1:mt,:) = Ac_tot(-1:mt,:)
    deallocate(Ac_tot)
    allocate(Ac_tot(-1:Nt+1,3))
    Ac_tot(:,:) = tmp2(:,:)
    deallocate(tmp2)

    ! data_local_Ac
    mt = min(Nt, ubound(data_local_Ac,3))
    allocate(tmp3(3,NXY_s:NXY_e,0:Nt))
    tmp3(:,:,:) = 0.d0
    tmp3(:,:,0:mt) = data_local_Ac(:,:,0:mt)
    deallocate(data_local_Ac)
    allocate(data_local_Ac(3,NXY_s:NXY_e,0:Nt))
    data_local_Ac(:,:,:) = tmp3(:,:,:)
    deallocate(tmp3)

    ! data_local_jm
    mt = min(Nt, ubound(data_local_jm,3))
    allocate(tmp3(3,NXY_s:NXY_e,0:Nt))
    tmp3(:,:,:) = 0.d0
    tmp3(:,:,0:mt) = data_local_jm(:,:,0:mt)
    deallocate(data_local_jm)
    allocate(data_local_jm(3,NXY_s:NXY_e,0:Nt))
    data_local_jm(:,:,:) = tmp3(:,:,:)
    deallocate(tmp3)

    ! data_vac_Ac
    mt = min(Nt, ubound(data_vac_Ac,3))
    allocate(tmp3(3,2,0:Nt))
    tmp3(:,:,:) = 0.d0
    tmp3(:,:,0:mt) = data_vac_Ac(:,:,0:mt)
    deallocate(data_vac_Ac)
    allocate(data_vac_Ac(3,2,0:Nt))
    data_vac_Ac(:,:,:) = tmp3
    deallocate(tmp3)

    ! data_out
    mt = min(Ndata_out_per_proc, ubound(data_out, 4))
    allocate(tmp4(16,NXvacL_m:NXvacR_m,NY_m+1,0:Ndata_out_per_proc))
    tmp4(:,:,:,:) = 0.d0
    tmp4(:,:,:,0:mt) = data_out(:,:,:,0:mt)
    deallocate(data_out)
    allocate(data_out(16,NXvacL_m:NXvacR_m,NY_m+1,0:Ndata_out_per_proc))
    data_out(:,:,:,:) = tmp4(:,:,:,:)
    deallocate(tmp4)
  end subroutine
end subroutine

subroutine prep_Reentrance_Read
  implicit none
  call prep_backup_values(.FALSE.)
end subroutine

subroutine prep_Reentrance_Write
  use global_variables, only: iter_now, entrance_iter
  implicit none
  entrance_iter = iter_now
  call prep_backup_values(.TRUE.)
end subroutine
