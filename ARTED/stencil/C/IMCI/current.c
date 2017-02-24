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

/* Hand-Code Vector processing for Knights Corner */

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

  __m512d CONJ = _mm512_set_pd(-1, 1, -1, 1, -1, 1, -1, 1);

  __m512d tt[3];
  for(n = 0 ; n < 3 ; ++n)
    tt[n] = _mm512_setzero_pd();

  __m512i nly = _mm512_set1_epi32(NLy);
  __m512i nlz = _mm512_set1_epi32(NLz);
  __m512i nyx = _mm512_mask_blend_epi32(0xF0F0, _mm512_set1_epi32(NLy    ), _mm512_set1_epi32(NLx    ));
#ifdef ARTED_DOMAIN_POWER_OF_TWO
  __m512i myx = _mm512_mask_blend_epi32(0xF0F0, _mm512_set1_epi32(NLy - 1), _mm512_set1_epi32(NLx - 1));
#endif

  __declspec(align(64)) int yx_table[16];
  __m512i  yx_org = _mm512_setr_epi32(1, 2, 3, 4, 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0);
  __m512i *yx     = (__m512i*) yx_table;

  __m512i dnyx = _mm512_add_epi32(nyx, yx_org);
  __m512i nlyz = _mm512_mask_mullo_epi32(nlz, 0xF0F0, nlz, nly);

#pragma noprefetch
#pragma novector
#ifdef ARTED_STENCIL_LOOP_BLOCKING
  for(bx = 0 ; bx < NLx ; bx += BX)
  for(by = 0 ; by < NLy ; by += BY)
  for(ix = bx ; ix < MIN(bx+BX,NLx) ; ++ix)
#else
  for(ix = 0 ; ix < NLx ; ++ix)
#endif
  {
    __m512i tix = _mm512_set1_epi32(ix);
#ifndef ARTED_DOMAIN_POWER_OF_TWO
    __m512i xp = _mm512_loadu_prefetch_epi32(modx + (ix + 1 + NLx));
            xp = _mm512_permute4f128_epi32(xp, _MM_PERM_CDAB);
#endif
#ifdef ARTED_STENCIL_LOOP_BLOCKING
  for(iy = by ; iy < MIN(by+BY,NLy) ; ++iy)
#else
  for(iy = 0 ; iy < NLy ; ++iy)
#endif
  {
    __m512i tiy = _mm512_set1_epi32(iy);
#ifndef ARTED_DOMAIN_POWER_OF_TWO
    __m512i yp  = _mm512_loadu_prefetch_epi32(mody + (iy + 1 + NLy));
    __m512i uyx = _mm512_mask_blend_epi32(0xF0F0, yp, xp);
#endif
    __m512i tyx = _mm512_mask_blend_epi32(0xF0F0, tiy, tix);

    e = &E[ix][iy][0];

    for(iz = 0 ; iz < NLz ; iz += 4)
    {
      __m512i tiz = _mm512_set1_epi32(iz);
#ifdef ARTED_DOMAIN_POWER_OF_TWO
      __m512i mm  = _mm512_sub_epi32(tyx, _mm512_and_epi32(_mm512_add_epi32(dnyx, tyx), myx));
#else
      __m512i mm  = _mm512_sub_epi32(tyx, uyx);
#endif
      *yx = _mm512_sub_epi32(tiz, _mm512_mullo_epi32(mm, nlyz));

      __m512d wm[4];
      __m512d wp[4];
      __m512d v1, v2, v3;

      // conj(e[iz])
      __m512d ez = _mm512_load_prefetch_pd(e + iz);
      __m512d w  = _mm512_mul_pd(ez, CONJ);

      /* z-dimension (unit stride) */
      {
        __m512i z0, z2;
#ifdef ARTED_DOMAIN_POWER_OF_TWO
        z2 = _mm512_load_prefetch_epi64(e + ((iz + 4 + NLz) & (NLz - 1)));
#else
        z2 = _mm512_load_prefetch_epi64(e + modz[iz + 4 + NLz]);
#endif

        wp[0] = (__m512d) _mm512_alignr_epi32(z2, (__m512i) ez,  4);
        wp[1] = (__m512d) _mm512_alignr_epi32(z2, (__m512i) ez,  8);
        wp[2] = (__m512d) _mm512_alignr_epi32(z2, (__m512i) ez, 12);
        wp[3] = (__m512d) z2;

        v2 = _mm512_setzero_pd();
#pragma unroll(4)
        for(n = 0 ; n < 4 ; ++n) {
          v2 = _mm512_fmadd_pd(_mm512_set1_pd(C[n+8]), wp[n], v2);
        }
        v1    = (__m512d) _mm512_shuffle_epi32((__m512i) v2, _MM_PERM_BADC);
        tt[2] = _mm512_fmadd_pd(w, v1, tt[2]);
      }

      /* y-dimension (NLz stride) */
      {
        wp[0] = _mm512_load_prefetch_pd(e + yx_table[0]);
        wp[1] = _mm512_load_prefetch_pd(e + yx_table[1]);
        wp[2] = _mm512_load_prefetch_pd(e + yx_table[2]);
        wp[3] = _mm512_load_prefetch_pd(e + yx_table[3]);

        v2 = _mm512_setzero_pd();
#pragma unroll(4)
        for(n = 0 ; n < 4 ; ++n) {
          v2 = _mm512_fmadd_pd(_mm512_set1_pd(C[n+4]), wp[n], v2);
        }
        v1    = (__m512d) _mm512_shuffle_epi32((__m512i) v2, _MM_PERM_BADC);
        tt[1] = _mm512_fmadd_pd(w, v1, tt[1]);
      }

      /* x-dimension (NLy*NLz stride)  */
      {
        wp[0] = _mm512_load_prefetch_pd(e + yx_table[4]);
        wp[1] = _mm512_load_prefetch_pd(e + yx_table[5]);
        wp[2] = _mm512_load_prefetch_pd(e + yx_table[6]);
        wp[3] = _mm512_load_prefetch_pd(e + yx_table[7]);

        v2 = _mm512_setzero_pd();
#pragma unroll(4)
        for(n = 0 ; n < 4 ; ++n) {
          v2 = _mm512_fmadd_pd(_mm512_set1_pd(C[n]), wp[n], v2);
        }
        v1    = (__m512d) _mm512_shuffle_epi32((__m512i) v2, _MM_PERM_BADC);
        tt[0] = _mm512_fmadd_pd(w, v1, tt[0]);
      }
    }  /* NLz */
  } /* NLy */
  } /* NLx */

  const __m512d two = _mm512_set1_pd(2);

  tt[0] = _mm512_mul_pd(tt[0], two);
  tt[1] = _mm512_mul_pd(tt[1], two);
  tt[2] = _mm512_mul_pd(tt[2], two);

  *F = _mm512_reduce_add_pd(tt[0]);
  *G = _mm512_reduce_add_pd(tt[1]);
  *H = _mm512_reduce_add_pd(tt[2]);
}
