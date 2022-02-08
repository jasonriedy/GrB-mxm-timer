#if !defined(GLOBALS_HEADER_)
#define GLOBALS_HEADER_

#if defined(NDEBUG)
#define DEBUG_PRINT(...)
#else
#define DEBUG_PRINT(...) do { fprintf (stderr, __VA_ARGS__); } while (0)
#endif

#define DIE(...) do { fprintf (stderr, __VA_ARGS__); exit (EXIT_FAILURE); } while (0)
#define DIE_PERROR(...) do { fprintf (stderr, __VA_ARGS__); perror (""); exit (EXIT_FAILURE); } while (0)

#define VERBOSELVL_PRINT(lvl, ...) do { if (verbose >= lvl) { fprintf (stderr, __VA_ARGS__); fflush (stderr); } } while (0)
#define VERBOSE_PRINT(...) VERBOSELVL_PRINT(1, __VA_ARGS__)

#if !defined(SCALE_MAX)
#define SCALE_MAX 40
#elif SCALE_MAX > 53
#error "Not using sufficient randomness in root sampling for SCALE_MAX > 53."
#endif

#if !defined(SCALE_DEFAULT)
#define SCALE_DEFAULT 12
#endif

#if !defined(EF_DEFAULT)
#define EF_DEFAULT 16
#endif

#if !defined(MAXWEIGHT_DEFAULT)
#define MAXWEIGHT_DEFAULT 255
#endif

#if !defined(NROOT_MAX)
#define NROOT_MAX 32
#endif

#if !defined(NROOT_DEFAULT)
#define NROOT_DEFAULT 8
#endif

#if !defined(A_DEFAULT)
#define A_DEFAULT 0.55
#endif

#if !defined(B_DEFAULT)
#define B_DEFAULT 0.1
#endif

#if !defined(NOISEFACT_DEFAULT)
#define NOISEFACT_DEFAULT 0.1
#endif

#if !defined(IN_GLOBALS_C)
extern const int SCALE, EF, NROOT, MAXWEIGHT;
extern int64_t NV, NE;
extern const int64_t Z, Zinv;
extern const uint64_t Z_hi, Z_low, Zinv_hi, Zinv_low;
extern const float A, B, NOISEFACT;
extern const int SCALE_BIG_THRESH;
extern int gen_tree;
#endif

void init_globals (int, int, int, int, float, float, float, int);

#endif /* GLOBALS_HEADER_ */
