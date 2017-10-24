program main
  use salmon_global
  use salmon_parallel
  use inputoutput
  implicit none

  call setup_parallel
  if (nproc_id_global == 0) then
    call print_software_version
  endif

  call read_input

  select case(iperiodic)
  case(0)
    call gceed
  case(3)
    call arted
  case default
    stop 'invalid iperiodic'
  end select

  call end_parallel

contains
  subroutine print_software_version
    implicit none
    include 'versionf.h'
    print '(A)',       '##############################################################################'
    print '(A)',       '# SALMON: Scalable Ab-initio Light-Matter simulator for Optics and Nanoscience'
    if (GIT_FOUND) then
      print '(A)',       '#'
      print '(A,A,A,A)', '#   [Git revision] ', GIT_COMMIT_HASH, ' in ', GIT_BRANCH
    endif
    print '(A)',       '##############################################################################'
  end subroutine
end program main
