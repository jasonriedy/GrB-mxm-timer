/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010-2011,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#if !defined(COMPAT_HEADER_)
#define COMPAT_HEADER_

#define _FILE_OFFSET_BITS 64
#define _THREAD_SAFE
#if !defined(_XOPEN_SOURCE)
#define _XOPEN_SOURCE 600
#endif
#define _XOPEN_SOURCE_EXTENDED
//#define _DEFAULT_SOURCE

#if __STDC_VERSION__ >= 199901L
#include <inttypes.h>
#elif defined(__MTA__)
#include <stdint.h>
#define PRId64 "d"
#define SCNd64 "d"
#else
#warning "Defining long as int64_t."
typedef long int64_t;
typedef unsigned uint32_fast_t;
#define PRId64 "ld"
#define SCNd64 "ld"
#if !defined(restrict)
#define restrict
#endif
#endif

#define parfor for

#if !defined(_OPENMP)
#define OMP(x)
#if defined(__GNUC__)
static int omp_get_thread_num (void) __attribute__((unused));
static int omp_get_num_threads (void) __attribute__((unused));
int omp_get_thread_num (void) { return 0; }
int omp_get_num_threads (void) { return 1; }
#else
static int omp_get_thread_num (void) { return 0; }
static int omp_get_num_threads (void) { return 1; }
#endif
#else
#define OMP_(x) _Pragma(#x)
#define OMP(x) OMP_(omp x)
#include <omp.h>
#undef parfor
#define parfor OMP(parallel for) for
#endif

#if defined(__cilk)
#include <cilk/cilk.h>
#undef parfor
#define parfor cilk_for
#endif

#if defined(__GNUC__)
#define CONST_FN_ATTR __attribute__((const))
#define UNUSED_FN_ATTR __attribute__((unused))
#else
#define CONST_FN_ATTR
#define UNUSED_FN_ATTR
#endif

#endif /* COMPAT_HEADER_ */
