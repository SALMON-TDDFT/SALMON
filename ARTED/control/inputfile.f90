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
module inputfile
  implicit none
  
  integer, parameter :: fh_namelist = 901
  integer, parameter :: fh_atomic_spiecies = 902
  integer, parameter :: fh_atomic_positions = 903
  integer, parameter :: fh_reentrance = 904

  integer :: inml_control
  integer :: inml_system
  integer :: inml_incident
  integer :: inml_propagation
  integer :: inml_rgrid
  integer :: inml_kgrid
  integer :: inml_tstep
  integer :: inml_electrons
  integer :: inml_pseudo
  integer :: inml_response
  integer :: inml_multiscale
  
  
  
contains



  subroutine extract_stdin()
    use communication, only: comm_is_root, comm_sync_all
    implicit none
    
    integer :: cur = fh_namelist
    integer :: ret = 0
    character(100) :: buff, text
    
    if (comm_is_root()) then
      open(fh_namelist, file='.namelist.tmp', status='replace')
      open(fh_atomic_spiecies, file='.atomic_spiecies.tmp', status='replace')
      open(fh_atomic_positions, file='.atomic_positions.tmp', status='replace')
      open(fh_reentrance, file='.reenetrance.tmp', status='replace')
      
      do while (.true.)
        read(*, '(a)', iostat=ret) buff
        if (ret < 0) then
          exit
        else
          text = trim(adjustl(buff))
          ! Comment lines
          if (text(1:1) == '!') cycle
          ! Beginning of 'atomic_species' part
          if (text == '&atomic_spiecies') then
            cur = fh_atomic_spiecies
            cycle
          end if
          ! Beginning of 'atomic_positions' part
          if (text == '&atomic_positions') then
            cur = fh_atomic_positions
            cycle
          end if
          ! Beginning of 'atomic_species' part
          if (text == '&reentrance') then
            cur = fh_reentrance
            cycle
          end if
          ! End of 'atomic_(spiecies|positions)' part
          if ((text == '/') .and. (cur /= fh_namelist)) then
            cur = fh_namelist
            cycle
          end if
          
          write(cur, '(a)') text
        end if
      end do
      close(fh_namelist)
      close(fh_atomic_positions)
      close(fh_atomic_spiecies)
      close(fh_reentrance)
    end if
    call comm_sync_all()
    return
  end subroutine extract_stdin

  
  
  subroutine set_default_param()
    use Global_Variables
    implicit none
    ! control
    entrance_option = 'new'
    Time_shutdown = 86400
    backup_frequency = 0
    entrance_iter = 0
    SYSname = ''
    directory = './'
    ! system
    functional = ''
    cval = 1
    aL = 0
    ax = 0
    ay = 0
    az = 0
    Sym = -1
    crystal_structure = ''
    NB = 0
    Nelec = 0
    ext_field = ''
    MD_option = 'N'
    AD_RHO = 'N'
    NE = 0
    NI = 0
    ! rgrid
    Nd = 4
    NLx = 0
    NLy = 0
    NLz = 0
    ! kgrid
    NKx = 0
    NKy = 0
    NKz = 0
    file_kw = ''
    ! tstep
    Nt = -1
    dt = 0
    ! pseudo
    ps_format = ''
    PSmask_option = 'n'
    alpha_mask = 0.8d0
    gamma_mask = 1.8d0
    eta_mask = 15
    ! electrons
    NEwald = 4
    aEwald = 0.5d0
    KbTev = 0
    Ncg = 5
    Nmemory_MB = 8
    alpha_MB = 0.75d0
    FSset_option = 'N'
    NFSset_start = 75
    NFSset_every = 25
    Nscf = 0
    ! incident
    Longi_Trans = 'Tr'
    dAc = 0
    AE_shape = ''
    IWcm2_1 = 0
    tpulsefs_1 = 0
    omegaev_1 = 0
    phi_CEP_1 = 0
    Epdir_1 = (/0,0,0/)
    IWcm2_2 = 0
    tpulsefs_2 = 0
    omegaev_2 = 0
    phi_CEP_2 = 0
    Epdir_2 = (/0,0,0/)
    T1_T2fs = 0
    ! response
    Nomega = -1
    domega = 0
    ! propagation
    propagator = 'default'
    ! multiscale
    FDTDdim = ""
    TwoD_shape = ""
    NX_m = 0
    NY_m = 0
    HX_m = 0
    HY_m = 0
    NKsplit = 0
    NXYsplit = 0
    NXvacL_m = 0
    NXvacR_m = 0
    return
  end subroutine set_default_param

  
  subroutine read_namelist()
    use communication, only: comm_is_root, comm_bcast, comm_sync_all, proc_group
    use Global_Variables
    implicit none
    
    namelist/control/ &
            & entrance_option, &
            & Time_shutdown, &
            & backup_frequency, &
            & entrance_iter, &
            & SYSname, &
            & directory
    namelist/system/ &
            & functional, &
            & cval, &
            & aL, &
            & ax, &
            & ay, &
            & az, &
            & Sym, &
            & crystal_structure, &
            & NB, &
            & Nelec, &
            & ext_field, &
            & MD_option, &
            & AD_RHO, &
            & NE, &
            & NI
    namelist/rgrid/ &
            & Nd, &
            & NLx, &
            & NLy, &
            & NLz
    namelist/kgrid/ &
            & NKx, &
            & NKy, &
            & NKz, &
            & file_kw
    namelist/tstep/ &
            & Nt, &
            & dt
    namelist/pseudo/ &
            & ps_format, &
            & PSmask_option, &
            & alpha_mask, &
            & gamma_mask, &
            & eta_mask
    namelist/electrons/ &
            & NEwald, &
            & aEwald, &
            & KbTev, &
            & Ncg, &
            & Nmemory_MB, &
            & alpha_MB, &
            & FSset_option, &
            & NFSset_start, &
            & NFSset_every, &
            & Nscf
    namelist/incident/ &
            & Longi_Trans, &
            & dAc, &
            & AE_shape, &
            & IWcm2_1, &
            & tpulsefs_1, &
            & omegaev_1, &
            & phi_CEP_1, &
            & Epdir_1, &
            & IWcm2_2, &
            & tpulsefs_2, &
            & omegaev_2, &
            & phi_CEP_2, &
            & Epdir_2, &
            & T1_T2fs
    namelist/propagation/ &
            & propagator
    namelist/response/ &
            & Nomega, &
            & domega
    namelist/multiscale/ &
            & FDTDdim, &
            & TwoD_shape, &
            & NX_m, &
            & NY_m, &
            & HX_m, &
            & HY_m, &
            & NKsplit, &
            & NXYsplit, &
            & NXvacL_m, &
            & NXvacR_m
      
    if (comm_is_root()) then
      call set_default_param()
      
      open(fh_namelist, file='.namelist.tmp', status='old')
      read(fh_namelist, nml=control, iostat=inml_control)
      rewind(fh_namelist)
      read(fh_namelist, nml=system, iostat=inml_system)
      rewind(fh_namelist)
      read(fh_namelist, nml=propagation, iostat=inml_propagation)
      rewind(fh_namelist)
      read(fh_namelist, nml=incident, iostat=inml_incident)
      rewind(fh_namelist)
      read(fh_namelist, nml=rgrid, iostat=inml_rgrid)
      rewind(fh_namelist)
      read(fh_namelist, nml=kgrid, iostat=inml_kgrid)
      rewind(fh_namelist)
      read(fh_namelist, nml=tstep, iostat=inml_tstep)
      rewind(fh_namelist)
      read(fh_namelist, nml=electrons, iostat=inml_electrons)
      rewind(fh_namelist)
      read(fh_namelist, nml=pseudo, iostat=inml_pseudo)
      rewind(fh_namelist)
      read(fh_namelist, nml=response, iostat=inml_response)
      rewind(fh_namelist)
      read(fh_namelist, nml=multiscale, iostat=inml_multiscale)
      close(fh_namelist)
    end if
    
    call comm_bcast(entrance_option, proc_group(1))
    call comm_bcast(Time_shutdown, proc_group(1))
    call comm_bcast(backup_frequency,proc_group(1))
    call comm_bcast(entrance_iter, proc_group(1))
    call comm_bcast(SYSname, proc_group(1))
    call comm_bcast(directory, proc_group(1))
    call comm_bcast(functional, proc_group(1))
    call comm_bcast(cval, proc_group(1))
    call comm_bcast(propagator,proc_group(1))
    call comm_bcast(aL, proc_group(1))
    call comm_bcast(ax, proc_group(1))
    call comm_bcast(ay, proc_group(1))
    call comm_bcast(az, proc_group(1))
    call comm_bcast(Sym, proc_group(1))
    call comm_bcast(crystal_structure, proc_group(1))
    call comm_bcast(NB, proc_group(1))
    call comm_bcast(Nelec, proc_group(1))
    call comm_bcast(ext_field, proc_group(1))
    call comm_bcast(MD_option, proc_group(1))
    call comm_bcast(AD_RHO, proc_group(1))
    call comm_bcast(NE, proc_group(1))
    call comm_bcast(NI, proc_group(1))
    call comm_bcast(Nd, proc_group(1))
    call comm_bcast(NLx, proc_group(1))
    call comm_bcast(NLy, proc_group(1))
    call comm_bcast(NLz, proc_group(1))
    call comm_bcast(NKx, proc_group(1))
    call comm_bcast(NKy, proc_group(1))
    call comm_bcast(NKz, proc_group(1))
    call comm_bcast(file_kw, proc_group(1))
    call comm_bcast(Nt, proc_group(1))
    call comm_bcast(dt, proc_group(1))
    call comm_bcast(ps_format, proc_group(1))
    call comm_bcast(PSmask_option, proc_group(1))
    call comm_bcast(alpha_mask, proc_group(1))
    call comm_bcast(gamma_mask, proc_group(1))
    call comm_bcast(eta_mask, proc_group(1))
    call comm_bcast(NEwald, proc_group(1))
    call comm_bcast(aEwald, proc_group(1))
    call comm_bcast(KbTev, proc_group(1))
    call comm_bcast(Ncg, proc_group(1))
    call comm_bcast(Nmemory_MB, proc_group(1))
    call comm_bcast(alpha_MB, proc_group(1))
    call comm_bcast(FSset_option, proc_group(1))
    call comm_bcast(NFSset_start, proc_group(1))
    call comm_bcast(NFSset_every, proc_group(1))
    call comm_bcast(Nscf, proc_group(1))
    call comm_bcast(Longi_Trans, proc_group(1))
    call comm_bcast(dAc, proc_group(1))
    call comm_bcast(AE_shape, proc_group(1))
    call comm_bcast(IWcm2_1, proc_group(1))
    call comm_bcast(tpulsefs_1, proc_group(1))
    call comm_bcast(omegaev_1, proc_group(1))
    call comm_bcast(phi_CEP_1, proc_group(1))
    call comm_bcast(Epdir_1, proc_group(1))
    call comm_bcast(IWcm2_2, proc_group(1))
    call comm_bcast(tpulsefs_2, proc_group(1))
    call comm_bcast(omegaev_2, proc_group(1))
    call comm_bcast(phi_CEP_2, proc_group(1))
    call comm_bcast(Epdir_2, proc_group(1))
    call comm_bcast(T1_T2fs, proc_group(1))
    call comm_bcast(Nomega, proc_group(1))
    call comm_bcast(domega, proc_group(1))
    call comm_bcast(FDTDdim, proc_group(1))
    call comm_bcast(TwoD_shape, proc_group(1))
    call comm_bcast(NX_m, proc_group(1))
    call comm_bcast(NY_m, proc_group(1))
    call comm_bcast(HX_m, proc_group(1))
    call comm_bcast(HY_m, proc_group(1))
    call comm_bcast(NKsplit, proc_group(1))
    call comm_bcast(NXYsplit, proc_group(1))
    call comm_bcast(NXvacL_m, proc_group(1))
    call comm_bcast(NXvacR_m, proc_group(1))
    call comm_sync_all()
    return
  end subroutine read_namelist


  subroutine read_atomic_spiecies()
    use communication, only: comm_is_root, comm_bcast, comm_sync_all, proc_group
    use Global_Variables, only: NE, Zatom, Lref
    implicit none
    integer :: i, index, Zatom_tmp, Lref_tmp
    
    allocate(Zatom(NE), Lref(NE))
    
    if (comm_is_root()) then
      open(fh_atomic_spiecies, file='.atomic_spiecies.tmp', status='old')
      do i=1, NE
        read(fh_atomic_spiecies, *) index, zatom_tmp, lref_tmp
        if (i == index) then
          Zatom(i) = Zatom_tmp
          Lref(i) = Lref_tmp
        else
          call err_finalize('atomic_spiecies is not ordered')
        end if
      end do
      close(fh_atomic_positions)
    end if
    call comm_bcast(Zatom, proc_group(1))
    call comm_bcast(Lref, proc_group(1))
    call comm_sync_all()
    return
  end subroutine


  subroutine read_atomic_positions()
    use communication, only: comm_is_root, comm_bcast, comm_sync_all, proc_group
    use Global_Variables, only: NI, Rion, Kion
    implicit none
    integer :: i, index, Kion_tmp
    real(8) :: Rion_tmp(3)
    
    allocate(Rion(3,NI), Kion(NI))
    
    if (comm_is_root()) then
      open(fh_atomic_positions, file='.atomic_positions.tmp', status='old')
      do i=1, NI
        read(fh_atomic_positions, *) index, Rion_tmp, Kion_tmp
        if (i == index) then
          Rion(:,i) = Rion_tmp
          Kion(i) = Kion_tmp
        else
          call err_finalize('atomic_positions is not ordered')
        end if
      end do
      close(fh_atomic_positions)
    end if
    call comm_bcast(Rion, proc_group(1))
    call comm_bcast(Kion, proc_group(1))
    call comm_sync_all()
    return
  end subroutine read_atomic_positions





  subroutine read_input()
    use global_variables, only: entrance_option
    implicit none
    
    call extract_stdin()
    call read_namelist()
    if(entrance_option == 'reentrance')return
    call read_atomic_spiecies()
    call read_atomic_positions()
    
    return
  end subroutine read_input
  
  
  
  subroutine dump_inputdata()
    use communication, only: comm_is_root, comm_sync_all
    use Global_Variables
    implicit none
    integer :: i
  
    if (comm_is_root()) then
      print '("#namelist: ",A,", status=",I1)', 'control', inml_control
      print '("#",4X,A,"=",A)', 'entrance_option', entrance_option
      print '("#",4X,A,"=",ES12.5)', 'Time_shutdown', Time_shutdown
      print '("#",4X,A,"=",I1)', 'entrance_iter', entrance_iter
      print '("#",4X,A,"=",A)', 'SYSname', SYSname
      print '("#",4X,A,"=",A)', 'directory', directory
      print '("#namelist: ",A,", status=",I1)', 'system', inml_system
      print '("#",4X,A,"=",A)', 'functional', functional
      print '("#",4X,A,"=",ES12.5)', 'cval', cval
      print '("#",4X,A,"=",ES12.5)', 'aL', aL
      print '("#",4X,A,"=",ES12.5)', 'ax', ax
      print '("#",4X,A,"=",ES12.5)', 'ay', ay
      print '("#",4X,A,"=",ES12.5)', 'az', az
      print '("#",4X,A,"=",I1)', 'Sym', Sym
      print '("#namelist: ",A,", status=",I1)', 'propagation', inml_propagation
      print '("#",4X,A,"=",A)', 'propagator', propagator
      print '("#",4X,A,"=",A)', 'crystal_structure', crystal_structure
      print '("#",4X,A,"=",I1)', 'NB', NB
      print '("#",4X,A,"=",I1)', 'Nelec', Nelec
      print '("#",4X,A,"=",A)', 'ext_field', ext_field
      print '("#",4X,A,"=",A)', 'MD_option', MD_option
      print '("#",4X,A,"=",A)', 'AD_RHO', AD_RHO
      print '("#",4X,A,"=",I1)', 'NE', NE
      print '("#",4X,A,"=",I1)', 'NI', NI
      print '("#namelist: ",A,", status=",I1)', 'rgrid', inml_rgrid
      print '("#",4X,A,"=",I1)', 'Nd', Nd
      print '("#",4X,A,"=",I1)', 'NLx', NLx
      print '("#",4X,A,"=",I1)', 'NLy', NLy
      print '("#",4X,A,"=",I1)', 'NLz', NLz
      print '("#namelist: ",A,", status=",I1)', 'kgrid', inml_kgrid
      print '("#",4X,A,"=",I1)', 'NKx', NKx
      print '("#",4X,A,"=",I1)', 'NKy', NKy
      print '("#",4X,A,"=",I1)', 'NKz', NKz
      print '("#",4X,A,"=",A)', 'file_kw', file_kw
      print '("#namelist: ",A,", status=",I1)', 'tstep', inml_tstep
      print '("#",4X,A,"=",I1)', 'Nt', Nt
      print '("#",4X,A,"=",ES12.5)', 'dt', dt
      print '("#namelist: ",A,", status=",I1)', 'pseudo', inml_pseudo
      print '("#",4X,A,"=",A)', 'ps_format', ps_format
      print '("#",4X,A,"=",A)', 'PSmask_option', PSmask_option
      print '("#",4X,A,"=",ES12.5)', 'alpha_mask', alpha_mask
      print '("#",4X,A,"=",ES12.5)', 'gamma_mask', gamma_mask
      print '("#",4X,A,"=",ES12.5)', 'eta_mask', eta_mask
      print '("#namelist: ",A,", status=",I1)', 'electrons', inml_electrons
      print '("#",4X,A,"=",I1)', 'NEwald', NEwald
      print '("#",4X,A,"=",ES12.5)', 'aEwald', aEwald
      print '("#",4X,A,"=",ES12.5)', 'KbTev', KbTev
      print '("#",4X,A,"=",I1)', 'Ncg', Ncg
      print '("#",4X,A,"=",I1)', 'Nmemory_MB', Nmemory_MB
      print '("#",4X,A,"=",ES12.5)', 'alpha_MB', alpha_MB
      print '("#",4X,A,"=",A)', 'FSset_option', FSset_option
      print '("#",4X,A,"=",I1)', 'NFSset_start', NFSset_start
      print '("#",4X,A,"=",I1)', 'NFSset_every', NFSset_every
      print '("#",4X,A,"=",I1)', 'Nscf', Nscf
      print '("#namelist: ",A,", status=",I1)', 'incident', inml_incident
      print '("#",4X,A,"=",A)', 'Longi_Trans', Longi_Trans
      print '("#",4X,A,"=",ES12.5)', 'dAc', dAc
      print '("#",4X,A,"=",A)', 'AE_shape', AE_shape
      print '("#",4X,A,"=",ES12.5)', 'IWcm2_1', IWcm2_1
      print '("#",4X,A,"=",ES12.5)', 'tpulsefs_1', tpulsefs_1
      print '("#",4X,A,"=",ES12.5)', 'omegaev_1', omegaev_1
      print '("#",4X,A,"=",ES12.5)', 'phi_CEP_1', phi_CEP_1
      print '("#",4X,A,"=",ES12.5,",",ES12.5,",",ES12.5)', 'Epdir_1', Epdir_1
      print '("#",4X,A,"=",ES12.5)', 'IWcm2_2', IWcm2_2
      print '("#",4X,A,"=",ES12.5)', 'tpulsefs_2', tpulsefs_2
      print '("#",4X,A,"=",ES12.5)', 'omegaev_2', omegaev_2
      print '("#",4X,A,"=",ES12.5)', 'phi_CEP_2', phi_CEP_2
      print '("#",4X,A,"=",ES12.5,",",ES12.5,",",ES12.5)', 'Epdir_2', Epdir_2
      print '("#",4X,A,"=",ES12.5)', 'T1_T2fs', T1_T2fs
      print '("#namelist: ",A,", status=",I1)', 'response', inml_response
      print '("#",4X,A,"=",I1)', 'Nomega', Nomega
      print '("#",4X,A,"=",ES12.5)', 'domega', domega
      print '("#namelist: ",A,", status=",I1)', 'multiscale', inml_multiscale
      print '("#",4X,A,"=",A)', 'FDTDdim', FDTDdim
      print '("#",4X,A,"=",A)', 'TwoD_shape', TwoD_shape
      print '("#",4X,A,"=",I1)', 'NX_m', NX_m
      print '("#",4X,A,"=",I1)', 'NY_m', NY_m
      print '("#",4X,A,"=",ES12.5)', 'HX_m', HX_m
      print '("#",4X,A,"=",ES12.5)', 'HY_m', HY_m
      print '("#",4X,A,"=",I1)', 'NKsplit', NKsplit
      print '("#",4X,A,"=",I1)', 'NXYsplit', NXYsplit
      print '("#",4X,A,"=",I1)', 'NXvacL_m', NXvacL_m
      print '("#",4X,A,"=",I1)', 'NXvacR_m', NXvacR_m
      
      print *, "#section: atomic_positions"
      do i=1, NE
        print '("#",4X,I1,",Zatom=",I1,",Lref=",I1)', i, Zatom(i), Lref(i)
      end do
      
      print *, "#section: atomic_positions"
      do i=1, NI
        print '("#",4X,I1,",Rion=",3ES12.5,",Kion=",I1)', i, Rion(:,i), Kion(i)
      end do
    end if
    call comm_sync_all()
    return
  end subroutine dump_inputdata
  
  
  
end module inputfile
