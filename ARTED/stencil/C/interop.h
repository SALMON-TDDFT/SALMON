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
#ifndef ARTED_INTEROP
#define ARTED_INTEROP

/* Stencil computation code with C supports Intel compiler only. */

/* currently Xeon CPUs and Xeon Phi are 64B cacheline. */
#ifndef CACHELINE_SIZE
# define CACHELINE_SIZE 64
#endif

#ifdef ARTED_ENABLE_SOFTWARE_PREFETCH
# define PREFETCH_L1(p, distance)  _mm_prefetch(((char const*)p) + distance, _MM_HINT_T0)
# define BUSY_PREFETCH_L1(p)       PREFETCH_L1(p, CACHELINE_SIZE)
#else
# define PREFETCH_L1(p, distance)
# define BUSY_PREFETCH_L1(p)
#endif

#if defined(__KNC__) || defined(__AVX512F__)
# define MEM_ALIGNED 64
# define VECTOR_SIZE 4
#else
# define MEM_ALIGNED 32
# define VECTOR_SIZE 2
#endif

#ifdef ARTED_DOMAIN_POWER_OF_TWO
# ifdef ENABLE_STENCIL_CODE_WITH_PADDING
#   define IDX(npt) ((ix - ((ix + (npt) + NLx) & (NLx-1))) * PNLy * PNLz)
#   define IDY(npt) ((iy - ((iy + (npt) + NLy) & (NLy-1))) * PNLz)
#   define IDZ(npt) ((iz - ((iz + (npt) + NLz) & (NLz-1))))
# else
#   define IDX(npt) ((ix - ((ix + (npt) + NLx) & (NLx-1))) * NLy * NLz)
#   define IDY(npt) ((iy - ((iy + (npt) + NLy) & (NLy-1))) * NLz)
#   define IDZ(npt) ((iz - ((iz + (npt) + NLz) & (NLz-1))))
# endif
#else
# ifdef ENABLE_STENCIL_CODE_WITH_PADDING
#   define IDX(npt) ((ix - modx[ix + (npt) + NLx]) * PNLy * PNLz)
#   define IDY(npt) ((iy - mody[iy + (npt) + NLy]) * PNLz)
#   define IDZ(npt) ((iz - modz[iz + (npt) + NLz]))
# else
#   define IDX(npt) ((ix - modx[ix + (npt) + NLx]) * NLy * NLz)
#   define IDY(npt) ((iy - mody[iy + (npt) + NLy]) * NLz)
#   define IDZ(npt) ((iz - modz[iz + (npt) + NLz]))
# endif
#endif /* ARTED_DOMAIN_POWER_OF_TWO */

#define MIN(n,m) (n < m ? n : m)

#include <immintrin.h>

#ifdef ARTED_STENCIL_LOOP_BLOCKING
#define BX   opt_variables_mp_stencil_blocking_x_
#define BY   opt_variables_mp_stencil_blocking_y_
#endif

#define NL   global_variables_mp_nl_
#define NLx  global_variables_mp_nlx_
#define NLy  global_variables_mp_nly_
#define NLz  global_variables_mp_nlz_
#define PNLx opt_variables_mp_pnlx_
#define PNLy opt_variables_mp_pnly_
#define PNLz opt_variables_mp_pnlz_

#define modx opt_variables_mp_modx_
#define mody opt_variables_mp_mody_
#define modz opt_variables_mp_modz_


#if defined(__AVX512F__) || (__KNC__)
inline
__m512i _mm512_load_prefetch_epi64(void const* c) {
  __m512i a = _mm512_load_epi64(c);
  BUSY_PREFETCH_L1(c);
  return a;
}

inline
__m512d _mm512_load_prefetch_pd(void const* c) {
  __m512d a = _mm512_load_pd(c);
  BUSY_PREFETCH_L1(c);
  return a;
}

inline
__m512i dcomplex_get_index(__m512i idx)
{
  const __m512i perm = _mm512_set_epi32(7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 0, 0);
  const __m512i one  = _mm512_set4_epi32(1, 0, 1, 0);
  __m512i x = _mm512_permutevar_epi32(perm, idx);
  __m512i y = _mm512_slli_epi32(x, 1);  /* x * 2 */
  return      _mm512_xor_si512(y, one); /* y + 1 or y + 0 */
}

inline
__m512d dcomplex_gather(void const* m, __m512i idx)
{
  const __m512i gidx = dcomplex_get_index(idx);
  return _mm512_i32loextgather_pd(gidx, m, _MM_UPCONV_PD_NONE, 8, _MM_HINT_NONE);
}
#endif

#if defined(__AVX512F__)
# include "./imci2avx512f.h"
#elif defined(__KNC__)
/* Knights Corner */
inline
__m512i _mm512_loadu_prefetch_epi32(int const* v) {
  __m512i w;
  w = _mm512_loadunpacklo_epi32(w, v + 0);
  w = _mm512_loadunpackhi_epi32(w, v + 16);
  BUSY_PREFETCH_L1(v);
  return w;
}

inline
__m512d dcast_to_dcmplx(double const *v) {
  const __m512i perm = _mm512_set_epi32(7, 6, 7, 6, 5, 4, 5, 4, 3, 2, 3, 2, 1, 0, 1, 0);
  __m512d w = _mm512_loadunpacklo_pd(_mm512_setzero_pd(), v);
  BUSY_PREFETCH_L1(v);
  return (__m512d) _mm512_permutevar_epi32(perm, (__m512i) w);
}

inline
__m512d dcomplex_mul(__m512d a, __m512d b) {
  __m512d ze = _mm512_setzero_pd();
  __m512d s0 = _mm512_swizzle_pd(a, _MM_SWIZ_REG_CDAB);  /* s0 = [a.i a.r] */
  __m512d re = _mm512_mask_blend_pd(0xAA, a, s0);        /* re = [a.r a.r] */
  __m512d im = _mm512_mask_sub_pd(a, 0x55, ze, s0);      /* im = [-a.i a.i] */
  __m512d t0 = _mm512_mul_pd(re, b);                     /* t0 = [a.r*b.r a.r*b.i] */
  __m512d s1 = _mm512_swizzle_pd(b, _MM_SWIZ_REG_CDAB);  /* s1 = [b.i b.r] */
  return       _mm512_fmadd_pd(im, s1, t0);              /* [-a.i*b.i+a.r*b.r a.i*b.r+a.r*b.i] */
}
#elif defined(__AVX__)
/* Sandy-Bridge or higher processors */
inline
__m256d dcast_to_dcmplx(double const *v) {
  __m256d a = _mm256_loadu_pd(v);
  __m256d b = _mm256_permute2f128_pd(a, a, 0x0);
  return _mm256_shuffle_pd(b, b, 0xC);
}
#endif

#endif /* ARTED_INTEROP */
