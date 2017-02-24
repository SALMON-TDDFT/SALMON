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

#include <complex.h>
#include "../interop.h"

extern int NLx, NLy, NLz;

#ifndef ARTED_DOMAIN_POWER_OF_TWO
extern int *modx, *mody, *modz;
#endif

#ifdef ARTED_STENCIL_LOOP_BLOCKING
extern int BX, BY;
#endif

void current_stencil_( double         const           C[restrict 12]
                     , double complex const           E[restrict NLx][NLy][NLz]
                     , double              * restrict F
                     , double              * restrict G
                     , double              * restrict H
) {
  double complex const* e;

  int ix, iy, iz, n;

#ifdef ARTED_STENCIL_LOOP_BLOCKING
  int bx, by;
#endif

  const __m256i CONJ = _mm256_set_epi64x(1LL << 63, 0, 1LL << 63, 0);

  __m256d tx = _mm256_setzero_pd();
  __m256d ty = _mm256_setzero_pd();
  __m256d tz = _mm256_setzero_pd();

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
    e = &E[ix][iy][0];

    for(iz = 0 ; iz < NLz ; iz += 2)
    {
#define STENCIL_CALC(PP,CC) \
      v2 = _mm256_mul_pd(_mm256_broadcast_sd(CC), PP); \
      v1 = _mm256_add_pd(v1, v2);

#define STENCIL(IDP,CC) \
      p = _mm256_load_pd((double *) (e + iz - IDP));     \
      STENCIL_CALC(p,CC) \

      __m256d m, p;
      __m256d v1, v2, v3, v4;

      // conj(e[iz])
      __m256d ez = _mm256_load_pd((double *)(e + iz));
      __m256d w = (__m256d) _mm256_xor_si256((__m256i) ez, CONJ);

      /* x-dimension (NLy*NLz stride)  */
      {
        v1 = _mm256_setzero_pd();
        STENCIL(IDX(1), C+0);
        STENCIL(IDX(2), C+1);
        STENCIL(IDX(3), C+2);
        STENCIL(IDX(4), C+3);
        v3 = _mm256_shuffle_pd(v1, v1, 0x05);
        v4 = _mm256_mul_pd(w, v3);
        tx = _mm256_add_pd(v4, tx);
      }

      /* y-dimension (NLz stride) */
      {
        v1 = _mm256_setzero_pd();
        STENCIL(IDY(1), C+4);
        STENCIL(IDY(2), C+5);
        STENCIL(IDY(3), C+6);
        STENCIL(IDY(4), C+7);
        v3 = _mm256_shuffle_pd(v1, v1, 0x05);
        v4 = _mm256_mul_pd(w, v3);
        ty = _mm256_add_pd(v4, ty);
      }

      /* z-dimension (unit stride) */
      {
        __m256d z2,z3,z5,z7;
#ifdef ARTED_DOMAIN_POWER_OF_TWO
        z2 = _mm256_load_pd((double *)(e + ((iz + 2 + NLz) & (NLz - 1))));
        z3 = _mm256_load_pd((double *)(e + ((iz + 4 + NLz) & (NLz - 1))));
#else
        z2 = _mm256_load_pd((double *)(e + modz[iz + 2 + NLz]));
        z3 = _mm256_load_pd((double *)(e + modz[iz + 4 + NLz]));
#endif
        z5 = _mm256_permute2f128_pd(ez, z2, 0x21);
        z7 = _mm256_permute2f128_pd(z2, z3, 0x21);

        v1 = _mm256_setzero_pd();
        STENCIL_CALC(z2, C+ 9);
        STENCIL_CALC(z3, C+11);
        STENCIL_CALC(z5, C+ 8);
        STENCIL_CALC(z7, C+10);
        v3 = _mm256_shuffle_pd(v1, v1, 0x05);
        v4 = _mm256_mul_pd(w, v3);
        tz = _mm256_add_pd(v4, tz);
      }
    }  /* NLz */
  } /* NLy */
  } /* NLx */

  const __m256d two = _mm256_set1_pd(2);

  tx = _mm256_mul_pd(tx, two);
  ty = _mm256_mul_pd(ty, two);
  tz = _mm256_mul_pd(tz, two);

  __declspec(align(32)) double r[4];

  _mm256_store_pd(r, tx);
  *F = r[0] + r[1] + r[2] + r[3];

  _mm256_store_pd(r, ty);
  *G = r[0] + r[1] + r[2] + r[3];

  _mm256_store_pd(r, tz);
  *H = r[0] + r[1] + r[2] + r[3];
}
