!
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!
!     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2011, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!         UNIVERSITY OF TSUKUBA
!         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!         E-MAIL: daisuke@cs.tsukuba.ac.jp
!
!
!     HEADER FILE FOR PARAMETERS
!
!     FORTRAN77 SOURCE PROGRAM
!
!     WRITTEN BY DAISUKE TAKAHASHI
!
! The maximum supported number of processors is 65536.
      PARAMETER (MAXNPU=65536)
! The maximum supported 2-D transform length is 65536.
      PARAMETER (NDA2=65536)
! The maximum supported 3-D transform length is 4096.
      PARAMETER (NDA3=4096)
! The parameter NBLK is a blocking parameter.
      PARAMETER (NBLK=16)
! The parameter NP is a padding parameter to avoid cache conflicts in
! the FFT routines.
      PARAMETER (NP=8)
! Size of L2 cache
      PARAMETER (L2SIZE=2097152)
