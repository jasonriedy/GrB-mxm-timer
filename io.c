#define _POSIX_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>
#include <errno.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <assert.h>

#include <intrinsics.h>

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

int
open_filename (const char* filename) {
    int f;
    if (args.dump_flag && args.run_powers_flag)
        DIE("Dumping not compatible with running the powers kernel\n");

    errno = 0;
    if (args.dump_flag) {
        if (!strcmp(filename, "-"))
            DIE("Cannot write to stdout");
        f = open (args.filename_arg, O_CREAT|O_WRONLY|O_TRUNC, 0666); //fopen (args.filename_arg, "w+");
    } else {
        if (!strcmp(filename, "-"))
            f = dup(0); // stdin
        else
            f = open (args.filename_arg, O_RDONLY); //fopen (args.filename_arg, "r");
    }
    if (errno)
        DIE_PERROR("Error opening \"%s\": ", args.filename_arg);
    return f;
}

GrB_Info
make_mtx_from_file (GrB_Matrix *A_out, GrB_Index * NV_out, GrB_Index * NE_out, int fd)
{
    GrB_Info info = GrB_SUCCESS;

    char name[1024]; // Yeah, bad idea;

    GrB_Index nrows, ncols, nvals;
    GrB_Index *off = NULL;
    GrB_Index *colind = NULL;
    uint64_t *val = NULL;

    // For text, "reopen" for the f* family.
    FILE *f = fdopen (fd, "r");
    if (!f)
        DIE_PERROR("Cannot move fd %d to FILE*", fd);

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

    fclose (f);

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
make_mtx_from_binfile (GrB_Matrix *A_out, GrB_Index * NV_out, GrB_Index * NE_out, int fd)
{
    GrB_Info info = GrB_SUCCESS;

    fprintf (stderr, "In make_mtx_from_binfile\n");

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

    volatile long clock1 = CLOCK();

    read (fd, tag, 8);
    if (!strcmp(tag, reverse_filetag))
        needs_bs = true;
    else if (strcmp(tag, filetag))
        DIE("Unrecognized file tag %s\n", tag);

    read (fd, &namelen, 8);
    if (namelen > sizeof(name))
        DIE("Name too long.\n");

    read (fd, name, namelen);
    DEBUG_PRINT("Name: %s\n", name);

    {
        uint64_t dims[3];
        read (fd, dims, 8*3);

        nrows = ensure_byteorder64(dims[0], needs_bs);
        ncols = ensure_byteorder64(dims[1], needs_bs);
        nvals = ensure_byteorder64(dims[2], needs_bs);
    }

    if (NV_out) *NV_out = nrows;
    if (NE_out) *NE_out = nvals;

    volatile long clock2 = CLOCK();

    off = malloc((nrows+1) * sizeof(*off));
    colind = malloc(nvals * sizeof(*colind));
    val = malloc(nvals * sizeof(*val));
    if (!off || !colind || !val)
        DIE_PERROR("Memory allocation failed reading matrix %s: ", name);

    volatile long clock3 = CLOCK();

    // The nrows+1 offsets
    read (fd, off, 8 * (nrows+1));
    if (needs_bs) {
        parfor (size_t k = 0; k <= nrows; ++k)
            off[k] = ensure_byteorder64(off[k], true);
    }
    if (off[nrows] != nvals)
        DIE("off[nrows] != nvals reading %s\n", name);

    volatile long clock4 = CLOCK();

    // Now the nvals colinds
    read (fd, colind, 8 * nvals);
    if (needs_bs) {
        parfor (size_t k = 0; k < nvals; ++k) {
            colind[k] = ensure_byteorder64(colind[k], true);
            if (colind[k] >= ncols)
                DIE("Out of range column reading %s\n", name);
        }
    }

    volatile long clock5 = CLOCK();

    // Finally the nvals values
    read (fd, val, 8 * nvals);
    if (needs_bs)
        parfor (size_t k = 0; k < nvals; ++k)
            val[k] = ensure_byteorder64(val[k], true);

    volatile long clock6 = CLOCK();

    GrB_Matrix A;

#if !defined(USE_SUITESPARSE)
    info = LGB_Matrix_import_CSR_UINT64 (&A, GrB_UINT64, nrows, ncols, off, colind, val, 0);
#else
    info = GxB_Matrix_import_CSR (&A, GrB_UINT64, nrows, ncols, &off, &colind, (void**)&val, (nrows+1)*sizeof(GrB_Index), nvals*sizeof(GrB_Index), nvals*sizeof(uint64_t), 0, 0, GrB_NULL);
#endif

    volatile long clock7 = CLOCK();

    if (info != GrB_SUCCESS)
        DIE ("Importing matrix %s failed: %ld\n", name, (long)info);

    *A_out = A;

#if !defined(USE_SUITESPARSE)
    free (val); free (colind); free (off);
#endif

    fprintf (stderr, "clocks %ld %ld %ld %ld %ld %ld %ld\n", clock1, clock2, clock3, clock4, clock5, clock6, clock7);

    return GrB_SUCCESS;
}

void
make_file_from_mtx (GrB_Matrix A, const char *name, int fd)
{
    GrB_Info info = GrB_SUCCESS;

    GrB_Index nrows, ncols;
    GrB_Index *off = NULL;
    GrB_Index *colind = NULL;
    uint64_t *val = NULL;

    FILE *f = fdopen (fd, "w");

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

    fclose (f);
}

void
make_binfile_from_mtx (GrB_Matrix A, const char *name, int fd)
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

    write (fd, filetag, 8);
    uint64_t namelen = strlen(name)+1;
    write (fd, &namelen, 8);
    write (fd, name, namelen);
    write (fd, &nrows, 8);
    write (fd, &ncols, 8);
    write (fd, &nnz, 8);

    write (fd, off, 8 * (nrows+1));

    if (nnz > 0) {
        write (fd, colind, 8 * nnz);
        write (fd, val, 8 * nnz);
    }

    free (val);
    free (colind);
    free (off);
}
