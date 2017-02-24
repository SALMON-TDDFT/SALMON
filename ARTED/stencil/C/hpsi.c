/*
 *  Copyright 2016 ARTED developers
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */
#include <complex.h>
#include "./interop.h"

extern int NLx, NLy, NLz;
extern int PNLx, PNLy, PNLz;
extern int *modx, *mody, *modz;

#ifdef ARTED_STENCIL_LOOP_BLOCKING
extern int BX, BY;
#endif

void hpsi1_rt_stencil_( double         const* restrict A_
                      , double         const           B[NLx][NLy][NLz]
                      , double         const* restrict C
                      , double         const* restrict D
                      , double complex const           E[restrict PNLx][PNLy][PNLz]
                      , double complex                 F[restrict PNLx][PNLy][PNLz]
)
{
  const double A = *A_;

  int ix, iy, iz;
  double complex v, w;
  double complex t[8];

#ifdef ARTED_STENCIL_LOOP_BLOCKING
  int bx, by;
#endif

  __assume_aligned(B, MEM_ALIGNED);
  __assume_aligned(E, MEM_ALIGNED);
  __assume_aligned(F, MEM_ALIGNED);

#undef IDX
#undef IDY
#undef IDZ

#ifdef ARTED_DOMAIN_POWER_OF_TWO
# ifdef __INTEL_COMPILER
  __assume(NLx % VECTOR_SIZE == 0);
  __assume(NLy % VECTOR_SIZE == 0);
  __assume(NLz % VECTOR_SIZE == 0);
# endif
# define IDX(dt) (ix+(dt)+NLx)&(NLx-1)][iy][iz
# define IDY(dt) ix][(iy+(dt)+NLy)&(NLy-1)][iz
# define IDZ(dt) ix][iy][(iz+(dt)+NLz)&(NLz-1)
#else
# define IDX(dt) modx[ix+(dt)+NLx]][iy][iz
# define IDY(dt) ix][mody[iy+(dt)+NLy]][iz
# define IDZ(dt) ix][iy][modz[iz+(dt)+NLz]
#endif

#ifdef ARTED_STENCIL_LOOP_BLOCKING
  for(bx = 0 ; bx < NLx ; bx += BX)
  for(by = 0 ; by < NLy ; by += BY)
  for(ix = bx ; ix < MIN(bx+BX,NLx) ; ++ix)
  for(iy = by ; iy < MIN(by+BY,NLy) ; ++iy)
#else
  for(ix = 0 ; ix < NLx ; ++ix)
  for(iy = 0 ; iy < NLy ; ++iy)
#endif
  {
#pragma simd
#pragma vector nontemporal(F)
    for(iz = 0 ; iz < NLz ; ++iz)
    {
      __assume_aligned(&B[ix][iy][0], MEM_ALIGNED);
      __assume_aligned(&E[ix][iy][0], MEM_ALIGNED);
      __assume_aligned(&F[ix][iy][0], MEM_ALIGNED);

      /* z-dimension (unit stride) */
      {
        t[0] = E[IDZ( 1)];
        t[1] = E[IDZ( 2)];
        t[2] = E[IDZ( 3)];
        t[3] = E[IDZ( 4)];
        t[4] = E[IDZ(-1)];
        t[5] = E[IDZ(-2)];
        t[6] = E[IDZ(-3)];
        t[7] = E[IDZ(-4)];

        v = (C[ 8] * (t[0] + t[4])
            +C[ 9] * (t[1] + t[5])
            +C[10] * (t[2] + t[6])
            +C[11] * (t[3] + t[7]));
        w = (D[ 8] * (t[0] - t[4])
            +D[ 9] * (t[1] - t[5])
            +D[10] * (t[2] - t[6])
            +D[11] * (t[3] - t[7]));
      }
      /* y-dimension (NLz stride) */
      {
        t[0] = E[IDY( 1)];
        t[1] = E[IDY( 2)];
        t[2] = E[IDY( 3)];
        t[3] = E[IDY( 4)];
        t[4] = E[IDY(-1)];
        t[5] = E[IDY(-2)];
        t[6] = E[IDY(-3)];
        t[7] = E[IDY(-4)];

        v += (C[ 4] * (t[0] + t[4])
             +C[ 5] * (t[1] + t[5])
             +C[ 6] * (t[2] + t[6])
             +C[ 7] * (t[3] + t[7]));
        w += (D[ 4] * (t[0] - t[4])
             +D[ 5] * (t[1] - t[5])
             +D[ 6] * (t[2] - t[6])
             +D[ 7] * (t[3] - t[7]));
      }
      /* x-dimension (NLy*NLz stride)  */
      {
        t[0] = E[IDX( 1)];
        t[1] = E[IDX( 2)];
        t[2] = E[IDX( 3)];
        t[3] = E[IDX( 4)];
        t[4] = E[IDX(-1)];
        t[5] = E[IDX(-2)];
        t[6] = E[IDX(-3)];
        t[7] = E[IDX(-4)];

        v += (C[ 0] * (t[0] + t[4])
             +C[ 1] * (t[1] + t[5])
             +C[ 2] * (t[2] + t[6])
             +C[ 3] * (t[3] + t[7]));
        w += (D[ 0] * (t[0] - t[4])
             +D[ 1] * (t[1] - t[5])
             +D[ 2] * (t[2] - t[6])
             +D[ 3] * (t[3] - t[7]));
      }

      F[ix][iy][iz] = B[ix][iy][iz] * E[ix][iy][iz] + A * E[ix][iy][iz] - 0.5 * v - _Complex_I * w;
    }
  }
}
