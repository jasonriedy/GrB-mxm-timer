#define _POSIX_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>
#include <errno.h>

#include <assert.h>

#include <GraphBLAS.h>
#if !defined(USE_SUITESPARSE)
#include <LucataGraphBLAS.h>
#endif

#include "cmdline.h"

#include "compat.h"
#include "globals.h"

extern struct gengetopt_args_info args;

static const char filetag[] = "mxmtimer";
static const char reverse_filetag[] = "remitmxm";

FILE*
open_filename (const char* filename) {
    FILE *f;
    if (args.dump_flag && args.run_powers_flag)
        DIE("Dumping not compatible with running the powers kernel\n");

    if (args.dump_flag) {
        if (!strcmp(filename, "-"))
            DIE("Cannot write to stdout");
        f = fopen (args.filename_arg, "w+");
    } else {
        if (!strcmp(filename, "-"))
            f = stdin;
        else
            f = fopen (args.filename_arg, "r");
    }
    if (!f)
        DIE_PERROR("Error opening \"%s\": ", args.filename_arg);
    return f;
}

GrB_Info
make_mtx_from_file (GrB_Matrix *A_out, GrB_Index * NV_out, GrB_Index * NE_out, FILE *f)
{
    GrB_Info info = GrB_SUCCESS;

    char name[1024]; // Yeah, bad idea;

    GrB_Index nrows, ncols, nvals;
    GrB_Index *off = NULL;
    GrB_Index *colind = NULL;
    uint64_t *val = NULL;

    fscanf (f, "%s", name);

    // Next input: nr, nc, nvals
    {
        long nrl, ncl, nvl;
        fscanf (f, "%ld %ld %ld", &nrl, &ncl, &nvl);
        nrows = nrl; ncols = ncl; nvals = nvl;
        DEBUG_PRINT("Read %s dims %ld %ld %ld\n", name, nrl, ncl, nvl);
    }

    if (NV_out) *NV_out = nrows;
    if (NE_out) *NE_out = nvals;

    off = malloc((nrows+1) * sizeof(*off));
    colind = malloc(nvals * sizeof(*colind));
    val = malloc(nvals * sizeof(*val));
    if (!off || !colind || !val)
        DIE_PERROR("Memory allocation failed reading matrix %s: ", name);

    // The nrows+1 offsets
    for (size_t k = 0; k <= nrows; ++k) {
        long tmp;
        fscanf (f, "%ld", &tmp);
        off[k] = tmp;
    }
    if (off[nrows] != nvals)
        DIE("off[nrows] != nvals reading %s\n", name);

    // Now the nvals colinds
    for (size_t k = 0; k < nvals; ++k) {
        long tmp;
        fscanf (f, "%ld", &tmp);
        if (tmp >= ncols)
            DIE("Out of range column reading %s\n", name);
        colind[k] = tmp;
    }

    // Finally the nvals values
    for (size_t k = 0; k < nvals; ++k)
        fscanf (f, "%" SCNu64, &val[k]);

    GrB_Matrix A;

#if !defined(USE_SUITESPARSE)
    info = LGB_Matrix_import_CSR_UINT64 (&A, GrB_UINT64, nrows, ncols, off, colind, val, 0);
#else
    info = GxB_Matrix_import_CSR (&A, GrB_UINT64, nrows, ncols, &off, &colind, (void**)&val, (nrows+1)*sizeof(GrB_Index), nvals*sizeof(GrB_Index), nvals*sizeof(uint64_t), 0, 0, GrB_NULL);
#endif

    if (info != GrB_SUCCESS)
        DIE("Importing matrix %s failed: %ld\n", name, (long)info);

    *A_out = A;

#if !defined(USE_SUITESPARSE)
    free (val); free (colind); free (off);
#endif

    return GrB_SUCCESS;
}

static inline uint64_t
ensure_byteorder64 (uint64_t x, bool needs_bs)
{
    if (!needs_bs) return x;
#ifdef USE_BUILTIN_BSWAP64
    return __builtin_bswap64(x);
#else
    x = (x >> 32) | (x << 32);
    x = ((x >> 16) & UINT64_C(0x0000FFFF0000FFFF)) |
        ((x & UINT64_C(0x0000FFFF0000FFFF)) << 16);
    x = ((x >>  8) & UINT64_C(0x00FF00FF00FF00FF)) |
        ((x & UINT64_C(0x00FF00FF00FF00FF)) <<  8);
    return x;
#endif
}

GrB_Info
make_mtx_from_binfile (GrB_Matrix *A_out, GrB_Index * NV_out, GrB_Index * NE_out, FILE *f)
{
    GrB_Info info = GrB_SUCCESS;

    char tag[9];
    memset (tag, 0, sizeof(tag));
    bool needs_bs = false;

    uint64_t namelen;
    char name[1025]; // Yeah, bad idea;
    memset (name, 0, sizeof(name));

    GrB_Index nrows, ncols, nvals;
    GrB_Index *off = NULL;
    GrB_Index *colind = NULL;
    uint64_t *val = NULL;

    fread (tag, 8, 1, f);
    if (!strcmp(tag, reverse_filetag))
        needs_bs = true;
    else if (strcmp(tag, filetag))
        DIE("Unrecognized file tag %s\n", tag);

    fread (&namelen, 8, 1, f);
    if (namelen > sizeof(name))
        DIE("Name too long.\n");

    fread (name, 1, namelen, f);
    DEBUG_PRINT("Name: %s\n", name);

    {
        uint64_t dims[3];
        fread (dims, 8, 3, f);

        nrows = ensure_byteorder64(dims[0], needs_bs);
        ncols = ensure_byteorder64(dims[1], needs_bs);
        nvals = ensure_byteorder64(dims[2], needs_bs);
    }

    if (NV_out) *NV_out = nrows;
    if (NE_out) *NE_out = nvals;

    off = malloc((nrows+1) * sizeof(*off));
    colind = malloc(nvals * sizeof(*colind));
    val = malloc(nvals * sizeof(*val));
    if (!off || !colind || !val)
        DIE_PERROR("Memory allocation failed reading matrix %s: ", name);

    // The nrows+1 offsets
    fread (off, 8, nrows+1, f);
    if (needs_bs) {
        parfor (size_t k = 0; k <= nrows; ++k)
            off[k] = ensure_byteorder64(off[k], true);
    }
    if (off[nrows] != nvals)
        DIE("off[nrows] != nvals reading %s\n", name);

    // Now the nvals colinds
    fread (colind, 8, nvals, f);
    if (needs_bs) {
        parfor (size_t k = 0; k < nvals; ++k) {
            colind[k] = ensure_byteorder64(colind[k], true);
            if (colind[k] >= ncols)
                DIE("Out of range column reading %s\n", name);
        }
    }

    // Finally the nvals values
    fread (val, 8, nvals, f);
    if (needs_bs)
        parfor (size_t k = 0; k < nvals; ++k)
            val[k] = ensure_byteorder64(val[k], true);

    GrB_Matrix A;

#if !defined(USE_SUITESPARSE)
    info = LGB_Matrix_import_CSR_UINT64 (&A, GrB_UINT64, nrows, ncols, off, colind, val, 0);
#else
    info = GxB_Matrix_import_CSR (&A, GrB_UINT64, nrows, ncols, &off, &colind, (void**)&val, (nrows+1)*sizeof(GrB_Index), nvals*sizeof(GrB_Index), nvals*sizeof(uint64_t), 0, 0, GrB_NULL);
#endif

    if (info != GrB_SUCCESS)
        DIE ("Importing matrix %s failed: %ld\n", name, (long)info);

    *A_out = A;

#if !defined(USE_SUITESPARSE)
    free (val); free (colind); free (off);
#endif

    return GrB_SUCCESS;
}

void
make_file_from_mtx (GrB_Matrix A, const char *name, FILE *f)
{
    GrB_Info info = GrB_SUCCESS;

    GrB_Index nrows, ncols;
    GrB_Index *off = NULL;
    GrB_Index *colind = NULL;
    uint64_t *val = NULL;

#if !defined(USE_SUITESPARSE)
    info = LGB_Matrix_export_CSR_UINT64 (A, NULL, &nrows, &ncols, &off, &colind, &val, NULL);

    if (info != GrB_SUCCESS)
        DIE("Export of %s failed: %ld\n", name, (long)info);
#else
    {
        GrB_Matrix dupA;
        info = GrB_Matrix_dup (&dupA, A);
        if (info != GrB_SUCCESS)
            DIE("Duplicating matrix %s: %ld\n", name, (long)info);

        GrB_Type type;
        GrB_Index off_size, colind_size, val_size;
        bool is_uniformed;
        bool jumbled;

        info = GxB_Matrix_export_CSR (&dupA, &type, &nrows, &ncols, &off, &colind, (void**)&val, &off_size, &colind_size, &val_size, &is_uniformed, &jumbled, GrB_NULL);
        if (info != GrB_SUCCESS)
            DIE("Export of %s failed: %ld\n", name, (long)info);

        /* Not only is this not implemented, it's wrong. */
        /* if (is_uniformed) */
        /*     DIE("%s is wearing a funny hat.\n", name); */
        if (jumbled)
            DIE("Assumed %s is sorted and it ain't.", name);

        GrB_free (&dupA);
    }
#endif

    GrB_Index nnz = off[nrows];
    DEBUG_PRINT("Writing name %s  dims %ld %ld %ld\n", name, (long)nrows, (long)ncols, (long)nnz);

    fprintf (f, "%s\n", name);
    fprintf (f, "%ld %ld %ld\n", (long)nrows, (long)ncols, (long)nnz);

    if (nrows > 0) {
        fprintf (f, "%ld", (long)off[0]);
        for (size_t k = 1; k <= nrows; ++k)
            fprintf (f, " %ld", (long)off[k]);
        fprintf (f, "\n");
    }

    if (nnz > 0) {
        for (size_t i = 0; i < nrows; ++i) {
            const size_t begin = off[i];
            const size_t end = off[i+1];
            if (begin != end) {
                fprintf (f, "%ld", (long)colind[begin]);
                for (size_t k = begin+1; k < end; ++k)
                    fprintf (f, " %ld", (long)colind[k]);
                fprintf (f, "\n\n");
            }
        }

        fprintf (f, "\n\n\n\n");

        for (size_t i = 0; i < nrows; ++i) {
            const size_t begin = off[i];
            const size_t end = off[i+1];
            if (begin != end) {
                fprintf (f, "%" PRIu64, val[begin]);
                for (size_t k = begin+1; k < end; ++k)
                    fprintf (f, " %" PRIu64, val[k]);
                fprintf (f, "\n\n");
            }
        }

        fprintf (f, "\n");
    }

    free (val);
    free (colind);
    free (off);
}

void
make_binfile_from_mtx (GrB_Matrix A, const char *name, FILE *f)
{
    GrB_Info info = GrB_SUCCESS;

    GrB_Index nrows, ncols;
    GrB_Index *off = NULL;
    GrB_Index *colind = NULL;
    uint64_t *val = NULL;

#if !defined(USE_SUITESPARSE)
    info = LGB_Matrix_export_CSR_UINT64 (A, NULL, &nrows, &ncols, &off, &colind, &val, NULL);

    if (info != GrB_SUCCESS)
        DIE("Export of %s failed: %ld\n", name, (long)info);
#else
    {
        GrB_Matrix dupA;
        info = GrB_Matrix_dup (&dupA, A);
        if (info != GrB_SUCCESS)
            DIE("Duplicating matrix %s: %ld\n", name, (long)info);

        GrB_Type type;
        GrB_Index off_size, colind_size, val_size;
        bool is_uniformed;
        bool jumbled;

        info = GxB_Matrix_export_CSR (&dupA, &type, &nrows, &ncols, &off, &colind, (void**)&val, &off_size, &colind_size, &val_size, &is_uniformed, &jumbled, GrB_NULL);
        if (info != GrB_SUCCESS)
            DIE("Export of %s failed: %ld\n", name, (long)info);

        /* Not only is this not implemented, it's wrong. */
        /* if (is_uniformed) */
        /*     DIE("%s is wearing a funny hat.\n", name); */
        if (jumbled)
            DIE("Assumed %s is sorted and it ain't.", name);

        GrB_free (&dupA);
    }
#endif

    GrB_Index nnz = off[nrows];
    DEBUG_PRINT("Writing name %s  dims %ld %ld %ld\n", name, (long)nrows, (long)ncols, (long)nnz);

    fwrite (filetag, 8, 1, f);
    uint64_t namelen = strlen(name)+1;
    fwrite (&namelen, 8, 1, f);
    fwrite (name, 1, namelen, f);
    fwrite (&nrows, 1, 8, f);
    fwrite (&ncols, 1, 8, f);
    fwrite (&nnz, 1, 8, f);

    fwrite (off, 8, nrows+1, f);

    if (nnz > 0) {
        fwrite (colind, 8, nnz, f);
        fwrite (val, 8, nnz, f);
    }

    free (val);
    free (colind);
    free (off);
}
