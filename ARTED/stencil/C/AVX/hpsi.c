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

/* Hand-Code Vector processing for Intel Ivy-Bridge (AVX) */

#define ENABLE_STENCIL_CODE_WITH_PADDING

#include <complex.h>
#include "../interop.h"

extern int PNLx, PNLy, PNLz;
extern int NLx, NLy, NLz;

#ifndef ARTED_DOMAIN_POWER_OF_TWO
extern int *modx, *mody, *modz;
#endif

#ifdef ARTED_STENCIL_LOOP_BLOCKING
extern int BX, BY;
#endif

void hpsi1_rt_stencil_( double         const* restrict A_
                      , double         const           B[restrict NLx][NLy][NLz]
                      , double         const           C[restrict 12]
                      , double         const           D[restrict 12]
                      , double complex const           E[restrict PNLx][PNLy][PNLz]
                      , double complex                 F[restrict PNLx][PNLy][PNLz]
)
{
  const  double         A = *A_;
  double         const* b;
  double complex const* e;
  double complex      * f;

  int ix, iy, iz, n;

#ifdef ARTED_STENCIL_LOOP_BLOCKING
  int bx, by;
#endif

  const __m256d at   = _mm256_set1_pd(A);
  const __m256d HALF = _mm256_set1_pd(-0.5);
  const __m256i INV  = _mm256_set_epi64x(1LL << 63, 0, 1LL << 63, 0);

#pragma novector
#ifdef ARTED_STENCIL_LOOP_BLOCKING
  for(bx = 0 ; bx < NLx ; bx += BX)
  for(by = 0 ; by < NLy ; by += BY)
  for(ix = bx ; ix < MIN(bx+BX,NLx) ; ++ix) {
  for(iy = by ; iy < MIN(by+BY,NLy) ; ++iy)
#else
  for(ix = 0 ; ix < NLx ; ++ix) {
  for(iy = 0 ; iy < NLy ; ++iy)
#endif
  {
    b = &B[ix][iy][0];
    e = &E[ix][iy][0];
    f = &F[ix][iy][0];

    for(iz = 0 ; iz < NLz ; iz += 2)
    {
      __m256d bt;
      __m256d tt = _mm256_setzero_pd();
      __m256d ut = _mm256_setzero_pd();
      __m256d ez = _mm256_load_pd((double *)(e + iz));

#define STENCIL_CALC(MM,PP,CC,DD) \
      v2 = _mm256_mul_pd(_mm256_broadcast_sd(DD),  _mm256_sub_pd(PP, MM)); \
      v1 = _mm256_mul_pd(_mm256_broadcast_sd(CC),  _mm256_add_pd(PP, MM)); \
      ut = _mm256_add_pd(ut, v2); \
      tt = _mm256_add_pd(tt, v1);

#define STENCIL(IDM,IDP,CC,DD) \
      m = _mm256_load_pd((double *) (e + iz - IDM));     \
      p = _mm256_load_pd((double *) (e + iz - IDP));     \
      STENCIL_CALC(m,p,CC,DD) \

      __m256d m, p;
      __m256d v1, v2, v3, v4, v5, v6;

      /* x-dimension (NLy*NLz stride)  */
      {
        STENCIL(IDX(-1), IDX(1), C+0, D+0);
        STENCIL(IDX(-2), IDX(2), C+1, D+1);
        STENCIL(IDX(-3), IDX(3), C+2, D+2);
        STENCIL(IDX(-4), IDX(4), C+3, D+3);
      }

      /* y-dimension (NLz stride) */
      {
        STENCIL(IDY(-1), IDY(1), C+4, D+4);
        STENCIL(IDY(-2), IDY(2), C+5, D+5);
        STENCIL(IDY(-3), IDY(3), C+6, D+6);
        STENCIL(IDY(-4), IDY(4), C+7, D+7);
      }

      /* z-dimension (unit stride) */
      {
        __m256d z0,z1,z2,z3,z4,z5,z6,z7;
#ifdef ARTED_DOMAIN_POWER_OF_TWO
        z1 = _mm256_load_pd((double *)(e + ((iz - 2 + NLz) & (NLz - 1))));
        z2 = _mm256_load_pd((double *)(e + ((iz + 2 + NLz) & (NLz - 1))));
        z0 = _mm256_load_pd((double *)(e + ((iz - 4 + NLz) & (NLz - 1))));
        z3 = _mm256_load_pd((double *)(e + ((iz + 4 + NLz) & (NLz - 1))));
#else
        z1 = _mm256_load_pd((double *)(e + modz[iz - 2 + NLz]));
        z2 = _mm256_load_pd((double *)(e + modz[iz + 2 + NLz]));
        z0 = _mm256_load_pd((double *)(e + modz[iz - 4 + NLz]));
        z3 = _mm256_load_pd((double *)(e + modz[iz + 4 + NLz]));
#endif
        z6 = _mm256_permute2f128_pd(z0, z1, 0x21);
        z4 = _mm256_permute2f128_pd(z1, ez, 0x21);
        z5 = _mm256_permute2f128_pd(ez, z2, 0x21);
        z7 = _mm256_permute2f128_pd(z2, z3, 0x21);

        STENCIL_CALC(z1, z2, C+ 9, D+ 9);
        STENCIL_CALC(z0, z3, C+11, D+11);
        STENCIL_CALC(z4, z5, C+ 8, D+ 8);
        STENCIL_CALC(z6, z7, C+10, D+10);
      }

      bt = dcast_to_dcmplx(b + iz);
      v6 = _mm256_shuffle_pd(ut, ut, 0x05);
      v4 = _mm256_add_pd(at, bt);
      v5 = (__m256d) _mm256_xor_si256((__m256i) v6, INV);
      v2 = _mm256_mul_pd(HALF, tt);
      v3 = _mm256_mul_pd(ez, v4);
      v1 = _mm256_add_pd(v2, v3);
      v6 = _mm256_add_pd(v1, v5);

      _mm256_stream_pd((double *) (f + iz), v6);
    } /* NLz */
  } /* NLy */
  } /* NLx */
}
