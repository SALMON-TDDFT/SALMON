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

#ifndef ARTED_IMCI_TO_AVX512F_HEADER
#define ARTED_IMCI_TO_AVX512F_HEADER

inline
__m512d dcast_to_dcmplx(double const *v) {
  __m512d w = _mm512_maskz_expandloadu_pd(0x55, v);
  BUSY_PREFETCH_L1(v);
  return _mm512_movedup_pd(w);
}

inline
__m512i _mm512_loadu_prefetch_epi32(int const* v) {
  __m512i w = _mm512_loadu_si512(v);
  BUSY_PREFETCH_L1(v);
  return w;
}

inline
__m512d dcomplex_mul(__m512d a, __m512d b) {
  __m512d ar = _mm512_movedup_pd(a);
  __m512d ai = _mm512_unpackhi_pd(a, a);
  __m512d bs = _mm512_shuffle_pd(b, b, 0x55);
          ai = _mm512_mul_pd(ai, bs);
  return       _mm512_fmaddsub_pd(b, ar, ai);
}

#define _mm512_permute4f128_epi32(a, p) _mm512_shuffle_i32x4(a, a, p)
#define _mm512_storenrngo_pd            _mm512_stream_pd

#endif
