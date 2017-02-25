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

#define TUNING_COMPLEX_MUL
#define TUNING_Z_MEM_LOAD
#define TUNING_INDEX_CALC

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

  __m512d at   = _mm512_set1_pd(A);
  __m512d HALF = _mm512_set1_pd(-0.5);
#ifdef TUNING_COMPLEX_MUL
  __m512i INV  = _mm512_set4_epi64(1LL << 63, 0, 1LL << 63, 0);
#else
  __m512d ZI   = _mm512_set_pd(-1, 0, -1, 0, -1, 0, -1, 0);
#endif

  __declspec(align(64)) double G[12];
  for(n = 0 ; n < 12 ; ++n)
    G[n] = C[n] * -0.5;

  __m512i nly = _mm512_set1_epi32(PNLy);
  __m512i nlz = _mm512_set1_epi32(PNLz);
#ifdef ARTED_DOMAIN_POWER_OF_TWO
  __m512i myx = _mm512_mask_blend_epi32(0xFF00, _mm512_set1_epi32(NLy - 1), _mm512_set1_epi32(NLx - 1));
  __m512i nyx = _mm512_mask_blend_epi32(0xFF00, _mm512_set1_epi32(NLy    ), _mm512_set1_epi32(NLx    ));
#else
  __m512i nyx = _mm512_mask_blend_epi32(0xFF00, _mm512_set1_epi32(PNLy   ), _mm512_set1_epi32(PNLx   ));
#endif

  __declspec(align(64)) int yx_table[16];
  __m512i  yx_org = _mm512_setr_epi32(-4, -3, -2, -1, 1, 2, 3, 4, -4, -3, -2, -1, 1, 2, 3, 4);
  __m512i *yx     = (__m512i*) yx_table;

  __m512i dnyx = _mm512_add_epi32(nyx, yx_org);
  __m512i nlyz = _mm512_mask_mullo_epi32(nlz, 0xFF00, nlz, nly);

#ifndef TUNING_Z_MEM_LOAD
  __m512i mlz = _mm512_set1_epi32(NLz - 1);
  __m512i zm_org = _mm512_setr_epi32(-4, -3, -2, -1, 0, 1, 2, 0,
                                      1,  2,  3,  4, 5, 6, 7, 0);
  __m512i ze     = _mm512_setzero_epi32();
#endif

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
#ifdef TUNING_INDEX_CALC
    __m512i mxm = _mm512_loadu_prefetch_epi32(modx + (ix - 4 + NLx));
    __m512i mxp = _mm512_alignr_epi32(mxm, mxm, 1);
    __m512i xmp = _mm512_mask_blend_epi32(0xF0F0, mxm, mxp);
            xmp = _mm512_permute4f128_epi32(xmp, _MM_PERM_BADC);
#endif
#endif
#ifdef ARTED_STENCIL_LOOP_BLOCKING
  for(iy = by ; iy < MIN(by+BY,NLy) ; ++iy)
#else
  for(iy = 0 ; iy < NLy ; ++iy)
#endif
  {
    __m512i tiy = _mm512_set1_epi32(iy);
#ifndef ARTED_DOMAIN_POWER_OF_TWO
#ifdef TUNING_INDEX_CALC
    __m512i mym = _mm512_loadu_prefetch_epi32(mody + (iy - 4 + NLy));
    __m512i myp = _mm512_alignr_epi32(mym, mym, 1);
    __m512i ymp = _mm512_mask_blend_epi32(0xF0F0, mym, myp);
    __m512i uyx = _mm512_mask_blend_epi32(0xFF00, ymp, xmp);
#endif
#endif
    __m512i tyx = _mm512_mask_blend_epi32(0xFF00, tiy, tix);

    b = &B[ix][iy][0];
    e = &E[ix][iy][0];
    f = &F[ix][iy][0];

    for(iz = 0 ; iz < NLz ; iz += 4)
    {
      __m512i tiz = _mm512_set1_epi32(iz);
#ifdef ARTED_DOMAIN_POWER_OF_TWO
      __m512i mm  = _mm512_sub_epi32(tyx, _mm512_and_epi32(_mm512_add_epi32(dnyx, tyx), myx));
#else
#ifdef TUNING_INDEX_CALC
      __m512i mm  = _mm512_sub_epi32(tyx, uyx);
#else
      __m512i mm  = _mm512_sub_epi32(tyx, _mm512_rem_epi32(_mm512_add_epi32(dnyx, tyx), nyx));
#endif
#endif
      *yx = _mm512_sub_epi32(tiz, _mm512_mullo_epi32(mm, nlyz));

      __m512d ez = _mm512_load_prefetch_pd(e + iz);
      __m512d tt = _mm512_setzero_pd();
      __m512d ut = _mm512_setzero_pd();

      __m512d wm[4];
      __m512d wp[4];
      __m512d bt, v0, v1, v2, v3, v4;

      /* z-dimension (unit stride) */
      {
#ifdef TUNING_Z_MEM_LOAD
        __m512i z0, z2;
#ifdef ARTED_DOMAIN_POWER_OF_TWO
        z0 = _mm512_load_prefetch_epi64(e + ((iz - 4 + NLz) & (NLz - 1)));
        z2 = _mm512_load_prefetch_epi64(e + ((iz + 4 + NLz) & (NLz - 1)));
#else
#ifdef TUNING_INDEX_CALC
        z0 = _mm512_load_prefetch_epi64(e + modz[iz - 4 + NLz]);
        z2 = _mm512_load_prefetch_epi64(e + modz[iz + 4 + NLz]);
#else
        z0 = _mm512_load_prefetch_epi64(e + (iz - 4 + NLz) % NLz);
        z2 = _mm512_load_prefetch_epi64(e + (iz + 4 + NLz) % NLz);
#endif
#endif
        wm[3] = (__m512d) z0;
        wm[2] = (__m512d) _mm512_alignr_epi32((__m512i) ez, z0,  4);
        wm[1] = (__m512d) _mm512_alignr_epi32((__m512i) ez, z0,  8);
        wm[0] = (__m512d) _mm512_alignr_epi32((__m512i) ez, z0, 12);
        wp[0] = (__m512d) _mm512_alignr_epi32(z2, (__m512i) ez,  4);
        wp[1] = (__m512d) _mm512_alignr_epi32(z2, (__m512i) ez,  8);
        wp[2] = (__m512d) _mm512_alignr_epi32(z2, (__m512i) ez, 12);
        wp[3] = (__m512d) z2;
#else
        __m512i mz = _mm512_and_epi32(_mm512_add_epi32(_mm512_add_epi32(tiz, nlz), zm_org), mlz);

        wm[3] = dcomplex_gather( e, mz );
        wm[2] = dcomplex_gather( e, _mm512_alignr_epi32(ze, mz,  1) );
        wm[1] = dcomplex_gather( e, _mm512_alignr_epi32(ze, mz,  2) );
        wm[0] = dcomplex_gather( e, _mm512_alignr_epi32(ze, mz,  3) );
        wp[0] = dcomplex_gather( e, _mm512_alignr_epi32(ze, mz,  8) );
        wp[1] = dcomplex_gather( e, _mm512_alignr_epi32(ze, mz,  9) );
        wp[2] = dcomplex_gather( e, _mm512_alignr_epi32(ze, mz, 10) );
        wp[3] = dcomplex_gather( e, _mm512_alignr_epi32(ze, mz, 11) );
#endif

#pragma unroll(4)
        for(n = 0 ; n < 4 ; ++n) {
          v4 = _mm512_sub_pd(wp[n], wm[n]);
          v3 = _mm512_add_pd(wp[n], wm[n]);
          ut = _mm512_fmadd_pd(_mm512_set1_pd(D[n+8]), v4, ut);
          tt = _mm512_fmadd_pd(_mm512_set1_pd(G[n+8]), v3, tt);
        }
      }

      /* y-dimension (NLz stride) */
      {
        wm[3] = _mm512_load_prefetch_pd(e + yx_table[0]);
        wm[2] = _mm512_load_prefetch_pd(e + yx_table[1]);
        wm[1] = _mm512_load_prefetch_pd(e + yx_table[2]);
        wm[0] = _mm512_load_prefetch_pd(e + yx_table[3]);
        wp[0] = _mm512_load_prefetch_pd(e + yx_table[4]);
        wp[1] = _mm512_load_prefetch_pd(e + yx_table[5]);
        wp[2] = _mm512_load_prefetch_pd(e + yx_table[6]);
        wp[3] = _mm512_load_prefetch_pd(e + yx_table[7]);

#pragma unroll(4)
        for(n = 0 ; n < 4 ; ++n) {
          v4 = _mm512_sub_pd(wp[n], wm[n]);
          v3 = _mm512_add_pd(wp[n], wm[n]);
          ut = _mm512_fmadd_pd(_mm512_set1_pd(D[n+4]), v4, ut);
          tt = _mm512_fmadd_pd(_mm512_set1_pd(G[n+4]), v3, tt);
        }
      }

      /* x-dimension (NLy*NLz stride)  */
      {
        wm[3] = _mm512_load_prefetch_pd(e + yx_table[ 8]);
        wm[2] = _mm512_load_prefetch_pd(e + yx_table[ 9]);
        wm[1] = _mm512_load_prefetch_pd(e + yx_table[10]);
        wm[0] = _mm512_load_prefetch_pd(e + yx_table[11]);
        wp[0] = _mm512_load_prefetch_pd(e + yx_table[12]);
        wp[1] = _mm512_load_prefetch_pd(e + yx_table[13]);
        wp[2] = _mm512_load_prefetch_pd(e + yx_table[14]);
        wp[3] = _mm512_load_prefetch_pd(e + yx_table[15]);

#pragma unroll(4)
        for(n = 0 ; n < 4 ; ++n) {
          v4 = _mm512_sub_pd(wp[n], wm[n]);
          v3 = _mm512_add_pd(wp[n], wm[n]);
          ut = _mm512_fmadd_pd(_mm512_set1_pd(D[n]), v4, ut);
          tt = _mm512_fmadd_pd(_mm512_set1_pd(G[n]), v3, tt);
        }
      }

      bt = dcast_to_dcmplx(b + iz);
      v2 = _mm512_fmadd_pd(at, ez, tt);
#ifdef TUNING_COMPLEX_MUL
      v4 = (__m512d) _mm512_shuffle_epi32((__m512i) ut, _MM_PERM_BADC);
      v3 = (__m512d) _mm512_xor_si512((__m512i) v4, INV);
#else
      v3 = dcomplex_mul(ut, ZI);
#endif
      v1 = _mm512_add_pd(v2, v3);
      v0 = _mm512_fmadd_pd(bt, ez, v1);

      _mm512_storenrngo_pd(&f[iz], v0);
    } /* NLz */
  } /* NLy */
  } /* NLx */
}
