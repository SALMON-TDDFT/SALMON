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
subroutine print_optimize_message
  implicit none
  print *, 'Preprocessor: '

#ifdef ARTED_USE_TLOG
  print *, '  ARTED_USE_TLOG'
#endif
#ifdef ARTED_USE_PAPI
  print *, '  ARTED_USE_PAPI'
#endif
#ifdef ARTED_CURRENT_PREPROCESSING
  print *, '  ARTED_CURRENT_PREPROCESSING'
#endif
#ifdef ARTED_STENCIL_ORIGIN
  print *, '  ARTED_STENCIL_ORIGIN'
#endif
#ifdef ARTED_STENCIL_OPTIMIZED
  print *, '  ARTED_STENCIL_OPTIMIZED'
#endif
#ifdef ARTED_STENCIL_PADDING
  print *, '  ARTED_STENCIL_PADDING'
#endif
#ifdef ARTED_STENCIL_ENABLE_LOOP_BLOCKING
  print *, '  ARTED_STENCIL_ENABLE_LOOP_BLOCKING'
#endif
#ifdef ARTED_DOMAIN_POWER_OF_TWO
  print *, '  ARTED_DOMAIN_POWER_OF_TWO'
#endif
#ifdef ARTED_EXPLICIT_VECTORIZATION
  print *, '  ARTED_EXPLICIT_VECTORIZATION'
#endif
#ifdef ARTED_ENABLE_SOFTWARE_PREFETCH
  print *, '  ARTED_ENABLE_SOFTWARE_PREFETCH'
#endif
end subroutine
