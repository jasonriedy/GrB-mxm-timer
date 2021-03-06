#define _POSIX_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <errno.h>

#if defined(__CILK)
#include <cilk/cilk.h>
#define parfor cilk_for
#elif defined(_OPENMP)
#define parfor _Pragma(omp parallel for) for
#else
#define parfor for
#endif
#include <GraphBLAS.h>

#include "cmdline.h"

#include "globals.h"
#include "generator.h" // for make_edge
#include "prng.h" // for sample_roots
#include "hooks.h"

int verbose = 0;

#if defined(NDEBUG)
#define DEBUG_PRINT(...)
#else
#define DEBUG_PRINT(...) do { fprintf (stderr, __VA_ARGS__); } while (0)
#endif

#define VERBOSELVL_PRINT(lvl, ...) do { if (verbose >= lvl) { fprintf (stderr, __VA_ARGS__); fflush (stderr); } } while (0)
#define VERBOSE_PRINT(...) VERBOSELVL_PRINT(1, __VA_ARGS__)

static GrB_Info
make_A (GrB_Matrix *A, const GrB_Index NV, const GrB_Index NE, const GrB_Index NE_chunk_size)
{
     GrB_Info info = GrB_SUCCESS;

     GrB_Matrix tmpA;
     GrB_Index *I = NULL, *J = NULL;
     uint64_t *V = NULL;
     const GrB_Index nchunks = (NE + NE_chunk_size - 1)/NE_chunk_size;

     DEBUG_PRINT("new A\n");
     info = GrB_Matrix_new (&tmpA, GrB_UINT64, NV, NV);
     if (info != GrB_SUCCESS) return info;

     DEBUG_PRINT("new Atmp\n");
     info = GrB_Matrix_new (A, GrB_UINT64, NV, NV);
     if (info != GrB_SUCCESS) goto done;

     I = malloc (NE_chunk_size * sizeof (*I));
     J = malloc (NE_chunk_size * sizeof (*J));
     V = malloc (NE_chunk_size * sizeof (*V));
     if (!I || !J || !V) { info = GrB_OUT_OF_MEMORY; goto done; }

     for (GrB_Index ck = 0; ck < nchunks; ++ck) {
          GrB_Index ngen = NE_chunk_size;
          if (ck * NE_chunk_size + ngen > NE)
               ngen = NE - ck * NE_chunk_size;
          VERBOSELVL_PRINT(2, "  chunk %ld/%ld  %ld %ld\n", (long)ck+1, (long)nchunks, (long)NE, (long)ngen);
          edge_list_64 ((int64_t*)I, (int64_t*)J, V, ck*NE_chunk_size, ngen);
          info = GrB_Matrix_build (tmpA, I, J, V, ngen, GrB_FIRST_UINT64);
          if (info != GrB_SUCCESS) goto done;
          info = GrB_eWiseAdd (*A, GrB_NULL, GrB_NULL, GrB_FIRST_UINT64, *A, tmpA, GrB_DESC_R);
          if (info != GrB_SUCCESS) goto done;
          GrB_Matrix_clear (tmpA);
     }
     
done:
     if (V) free(V);
     if (J) free(J);
     if (I) free(I);
     GrB_free (&tmpA);
     return info;
}

static GrB_Info 
make_B (GrB_Matrix *B, GrB_Index NV, GrB_Index B_ncols, GrB_Index B_used_ncols, int B_nents_per_col)
{
     // Assume B_ncols << NV, and B_nents_per_col is small (e.g. 1 to ten).
     // So everything in one pass and sequentially here.

     GrB_Info info;

     const int64_t nroot = B_used_ncols * B_nents_per_col;
     GrB_Index I[nroot], J[nroot];
     uint64_t V[nroot];

     info = GrB_Matrix_new (B, GrB_UINT64, NV, B_ncols);
     if (info != GrB_SUCCESS) return info;

     // key doesn't really matter, but must be reproducible.
     sample_roots ((int64_t*)I, nroot, NV * B_used_ncols * B_nents_per_col);
     for (GrB_Index k = 0; k < nroot; ++k) {
          V[k] = 1;
          J[k] = k / B_nents_per_col;
     }

     info = GrB_Matrix_build (*B, I, J, V, nroot, GrB_FIRST_UINT64);
     if (info != GrB_SUCCESS) GrB_free (B);

     return info;
}

static GrB_Info
timed_loop (GrB_Matrix B, GrB_Matrix A, const int nhop)
{
     GrB_Info info;
     for (int k = 0; k < nhop; ++k) {
          info = GrB_mxm (B, GrB_NULL, GrB_NULL, GxB_MIN_FIRST_FP64, A, B, GrB_DESC_R);
          if (info != GrB_SUCCESS) return info;
     }
     GrB_wait();
     return info;
}

static struct gengetopt_args_info args;

int
main (int argc, char **argv) 
{
     VERBOSE_PRINT("LAUNCHED\n");
     if (0 != cmdline_parser (argc, argv, &args))
          exit (1);
     if (NULL != getenv ("VERBOSE")) {
          long lvl = strtol (getenv ("VERBOSE"), NULL, 10);
          if (lvl > 1) verbose = lvl;
     }
     if (args.verbose_given)
          verbose = args.verbose_arg;

     // Parse the khops arg
     int n_khops = 0;
     long *khops = malloc (strlen (args.khops_arg) * sizeof (*khops));
     if (!khops) {
          perror ("Cannot malloc khops");
          exit (1);
     }
     DEBUG_PRINT("parsing hops\n");
     char *saveptr = NULL, *token = NULL;
     char *inputptr = args.khops_arg;
     for (;;) {
          token = strtok_r (inputptr, " \n", &saveptr);
          DEBUG_PRINT("hop %d  %p %p %p\n", n_khops, token, inputptr, saveptr);
          inputptr = NULL;
          if (token == NULL) break;
          errno = 0;
          long hop = strtol (token, NULL, 10);
          if (hop <= 0) {
               fprintf (stderr, "Invalid hop value: %ld\n", hop);
               exit (1);
          }
          if (errno) {
               fprintf (stderr, "Error parsing hop %d", n_khops+1);
               perror ("");
               exit (1);
          }
          khops[n_khops++] = hop;
     }
     DEBUG_PRINT("done parsing hops\n");

     VERBOSE_PRINT("Starting GrB-mxm-timer\n");

     init_globals (args.scale_arg, args.edgefactor_arg, 255,
                   1, // unused
                   args.A_arg, args.B_arg, args.noisefact_arg);

     GrB_Info info;
     GrB_Matrix A, Bini, B;

     info = GrB_init (GrB_NONBLOCKING);
     if (info != GrB_SUCCESS) {
          fprintf (stderr, "Error initializing GraphBLAS: %ld\n", (long)info);
          exit (1);
     }

     VERBOSE_PRINT("Creating A... ");
     if (!args.no_time_A_flag) {
          hooks_set_attr_i64 ("scale", SCALE);
          hooks_set_attr_i64 ("edgefactor", EF);
          hooks_set_attr_f64 ("A", args.A_arg);
          hooks_set_attr_f64 ("B", args.B_arg);
          hooks_set_attr_f64 ("noisefact", NOISEFACT);
          hooks_region_begin ("Generating matrix A");
     }

     info = make_A (&A, NV, NE, args.NE_chunk_size_arg);

     double A_time = 0.0;
     if (!args.no_time_A_flag) A_time = hooks_region_end ();
     if (info != GrB_SUCCESS) {
          fprintf (stderr, "Error making A: %ld\n", (long)info);
          exit (1);
     }
     VERBOSE_PRINT("%g ms\n", A_time);

     VERBOSE_PRINT("Creating Bini... ");
     if (!args.no_time_B_flag) {
          hooks_set_attr_i64 ("b-ncols", args.b_ncols_arg);
          hooks_set_attr_i64 ("b-used-ncols", args.b_used_ncols_arg);
          hooks_set_attr_i64 ("b-nents-col", args.b_nents_col_arg);
          hooks_region_begin ("Generating Bini");
     }

     info = make_B (&Bini, NV, args.b_ncols_arg, args.b_used_ncols_arg, args.b_nents_col_arg);

     double Bini_time = 0.0;
     if (!args.no_time_B_flag) hooks_region_end ();
     if (info != GrB_SUCCESS) {
          fprintf (stderr, "Error making Bini: %ld\n", (long)info);
          exit (1);
     }
     VERBOSE_PRINT("%g ms\n", Bini_time);

     for (int k = 0; k < n_khops; ++k) {
          info = GrB_Matrix_dup (&B, Bini);
          if (info != GrB_SUCCESS) {
               fprintf (stderr, "Error copying B = Bini on hop value %d\n", k);
               exit (1);
          }

          VERBOSE_PRINT("Running hop #%d for %ld steps... ", k, khops[k]);
          if (!args.no_time_iter_flag) {
               hooks_set_attr_i64 ("khop", khops[k]);
               hooks_region_begin ("Iterating");
          }

          info = timed_loop (B, A, khops[k]);

          double iter_time = 0.0;
          if (!args.no_time_iter_flag) hooks_region_end ();
          VERBOSE_PRINT("%g ms\n", iter_time);

          GrB_free (&B);
     }

     VERBOSE_PRINT("DONE\n");
     GrB_finalize ();
}
