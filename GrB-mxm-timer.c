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

static GrB_Info
make_A (GrB_Matrix *A, const GrB_Index NV, const GrB_Index NE, const GrB_Index NE_chunk_size)
{
     GrB_Info info = GrB_SUCCESS;

     GrB_Matrix tmpA;
     GrB_Index *I = NULL, *J = NULL;
     uint64_t *V = NULL;
     const GrB_Index nchunks = (NE + NE_chunk_size - 1)/NE_chunk_size;

     info = GrB_Matrix_new (&tmpA, GrB_UINT64, NV, NV);
     if (info != GrB_SUCCESS) return info;

     info = GrB_Matrix_new (A, GrB_UINT64, NV, NV);
     if (info != GrB_SUCCESS) goto done;

     I = malloc (NE_chunk_size * sizeof (*I));
     J = malloc (NE_chunk_size * sizeof (*J));
     V = malloc (NE_chunk_size * sizeof (*V));
     if (!I || !J || !V) { info = GrB_OUT_OF_MEMORY; goto done; }

     for (GrB_Index ck = 0; ck < nchunks; ++ck) {
          if (verbose > 1) { printf ("  chunk %ld/%ld  ", (long)ck+1, (long)nchunks); fflush (stdout); }
          const GrB_Index ngen = (ck * NE_chunk_size <= NE? NE_chunk_size : NE - (ck-1) * NE_chunk_size);
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
          info = GrB_mxm (B, GrB_NULL, GrB_NULL, GrB_MIN_FIRST_SEMIRING_FP64, A, B, GrB_DESC_R);
          if (info != GrB_SUCCESS) return info;
     }
     GrB_wait();
     return info;
}

static struct gengetopt_args_info args;

int
main (int argc, char **argv) 
{
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
     char *saveptr = NULL, *token;
     for (;;) {
          token = strtok_r ((saveptr == NULL ? args.khops_arg : NULL), " \n", &saveptr);
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

     if (verbose) { printf ("Creating A... "); fflush (stdout); }
     hooks_set_attr_i64 ("scale", SCALE);
     hooks_set_attr_i64 ("edgefactor", EF);
     hooks_set_attr_f64 ("A", args.A_arg);
     hooks_set_attr_f64 ("B", args.B_arg);
     hooks_set_attr_f64 ("noisefact", NOISEFACT);
     hooks_region_begin ("Generating matrix A");
     {
          info = make_A (&A, NV, NE, args.NE_chunk_size_arg);
     }
     const double A_time = hooks_region_end ();
     if (info != GrB_SUCCESS) {
          fprintf (stderr, "Error making A: %ld\n", (long)info);
          exit (1);
     }
     if (verbose) { printf ("%g ms\n", A_time); fflush (stdout); }

     if (verbose) { printf ("Creating Bini... "); fflush (stdout); }
     hooks_set_attr_i64 ("b-ncols", args.b_ncols_arg);
     hooks_set_attr_i64 ("b-used-ncols", args.b_used_ncols_arg);
     hooks_set_attr_i64 ("b-nents-col", args.b_nents_col_arg);
     hooks_region_begin ("Generating Bini");
     {
          info = make_B (&Bini, NV, args.b_ncols_arg, args.b_used_ncols_arg, args.b_nents_col_arg);
     }
     const double Bini_time = hooks_region_end ();
     if (info != GrB_SUCCESS) {
          fprintf (stderr, "Error making Bini: %ld\n", (long)info);
          exit (1);
     }
     if (verbose) { printf ("%g ms\n", Bini_time); fflush (stdout); }

     for (int k = 0; k < n_khops; ++k) {
          info = GrB_Matrix_dup (&B, Bini);
          if (info != GrB_SUCCESS) {
               fprintf (stderr, "Error copying B = Bini on hop value %d\n", k);
               exit (1);
          }

          if (verbose) { printf ("Running hop #%d for %ld steps... ", k, khops[k]); fflush (stdout); }
          hooks_set_attr_i64 ("khop", khops[k]);
          hooks_region_begin ("Iterating");
          {
               info = timed_loop (B, A, khops[k]);
          }
          const double iter_time = hooks_region_end ();
          if (verbose) { printf ("%g ms\n", iter_time); fflush (stdout); }

          GrB_free (&B);
     }

     if (verbose) { printf ("DONE\n"); fflush (stdout); }
     GrB_finalize ();
}
