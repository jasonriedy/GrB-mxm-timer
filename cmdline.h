/** @file cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.23
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt */

#ifndef CMDLINE_H
#define CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE "GrB-mxm-timer"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "GrB-mxm-timer"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "0"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  int scale_arg;	/**< @brief Scale (log2 # vertices in A) (default='16').  */
  char * scale_orig;	/**< @brief Scale (log2 # vertices in A) original value given at command line.  */
  const char *scale_help; /**< @brief Scale (log2 # vertices in A) help description.  */
  int edgefactor_arg;	/**< @brief Edge factor, so # edges = ef * 2^scale (default='8').  */
  char * edgefactor_orig;	/**< @brief Edge factor, so # edges = ef * 2^scale original value given at command line.  */
  const char *edgefactor_help; /**< @brief Edge factor, so # edges = ef * 2^scale help description.  */
  float A_arg;	/**< @brief R-MAT upper left quadrant probability (default='0.55').  */
  char * A_orig;	/**< @brief R-MAT upper left quadrant probability original value given at command line.  */
  const char *A_help; /**< @brief R-MAT upper left quadrant probability help description.  */
  float B_arg;	/**< @brief R-MAT upper right & lower left quadrant probability (default='0.1').  */
  char * B_orig;	/**< @brief R-MAT upper right & lower left quadrant probability original value given at command line.  */
  const char *B_help; /**< @brief R-MAT upper right & lower left quadrant probability help description.  */
  float noisefact_arg;	/**< @brief Noise factor on each recursion (default='0.1').  */
  char * noisefact_orig;	/**< @brief Noise factor on each recursion original value given at command line.  */
  const char *noisefact_help; /**< @brief Noise factor on each recursion help description.  */
  int run_powers_flag;	/**< @brief Run powers of the generated A matrix rather than applying A to B (default=off).  */
  const char *run_powers_help; /**< @brief Run powers of the generated A matrix rather than applying A to B help description.  */
  char * filename_arg;	/**< @brief Filename to read/write for a CSR format.  */
  char * filename_orig;	/**< @brief Filename to read/write for a CSR format original value given at command line.  */
  const char *filename_help; /**< @brief Filename to read/write for a CSR format help description.  */
  int dump_flag;	/**< @brief Write a file to read (default=off).  */
  const char *dump_help; /**< @brief Write a file to read help description.  */
  int binary_flag;	/**< @brief File is in binary format (default=off).  */
  const char *binary_help; /**< @brief File is in binary format help description.  */
  int b_ncols_arg;	/**< @brief Number of columns in B (default='16').  */
  char * b_ncols_orig;	/**< @brief Number of columns in B original value given at command line.  */
  const char *b_ncols_help; /**< @brief Number of columns in B help description.  */
  int b_used_ncols_arg;	/**< @brief Number of columns actually used in the initial B (default='1').  */
  char * b_used_ncols_orig;	/**< @brief Number of columns actually used in the initial B original value given at command line.  */
  const char *b_used_ncols_help; /**< @brief Number of columns actually used in the initial B help description.  */
  int b_nents_col_arg;	/**< @brief Number of entries per column in the initial B (default='1').  */
  char * b_nents_col_orig;	/**< @brief Number of entries per column in the initial B original value given at command line.  */
  const char *b_nents_col_help; /**< @brief Number of entries per column in the initial B help description.  */
  char * khops_arg;	/**< @brief Number of iterations / hops (can be a space-delim list) (default='2 4 8').  */
  char * khops_orig;	/**< @brief Number of iterations / hops (can be a space-delim list) original value given at command line.  */
  const char *khops_help; /**< @brief Number of iterations / hops (can be a space-delim list) help description.  */
  int NE_chunk_size_arg;	/**< @brief Number of edges to generate in a chunk. (default='1048576').  */
  char * NE_chunk_size_orig;	/**< @brief Number of edges to generate in a chunk. original value given at command line.  */
  const char *NE_chunk_size_help; /**< @brief Number of edges to generate in a chunk. help description.  */
  int verbose_arg;	/**< @brief Provide status updates via stdout. (default='1').  */
  char * verbose_orig;	/**< @brief Provide status updates via stdout. original value given at command line.  */
  const char *verbose_help; /**< @brief Provide status updates via stdout. help description.  */
  int no_time_A_flag;	/**< @brief Do not time A (default=off).  */
  const char *no_time_A_help; /**< @brief Do not time A help description.  */
  int no_time_B_flag;	/**< @brief Do not time B (default=off).  */
  const char *no_time_B_help; /**< @brief Do not time B help description.  */
  int no_time_iter_flag;	/**< @brief Do not time iteration (default=off).  */
  const char *no_time_iter_help; /**< @brief Do not time iteration help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int scale_given ;	/**< @brief Whether scale was given.  */
  unsigned int edgefactor_given ;	/**< @brief Whether edgefactor was given.  */
  unsigned int A_given ;	/**< @brief Whether A was given.  */
  unsigned int B_given ;	/**< @brief Whether B was given.  */
  unsigned int noisefact_given ;	/**< @brief Whether noisefact was given.  */
  unsigned int run_powers_given ;	/**< @brief Whether run-powers was given.  */
  unsigned int filename_given ;	/**< @brief Whether filename was given.  */
  unsigned int dump_given ;	/**< @brief Whether dump was given.  */
  unsigned int binary_given ;	/**< @brief Whether binary was given.  */
  unsigned int b_ncols_given ;	/**< @brief Whether b-ncols was given.  */
  unsigned int b_used_ncols_given ;	/**< @brief Whether b-used-ncols was given.  */
  unsigned int b_nents_col_given ;	/**< @brief Whether b-nents-col was given.  */
  unsigned int khops_given ;	/**< @brief Whether khops was given.  */
  unsigned int NE_chunk_size_given ;	/**< @brief Whether NE-chunk-size was given.  */
  unsigned int verbose_given ;	/**< @brief Whether verbose was given.  */
  unsigned int no_time_A_given ;	/**< @brief Whether no-time-A was given.  */
  unsigned int no_time_B_given ;	/**< @brief Whether no-time-B was given.  */
  unsigned int no_time_iter_given ;	/**< @brief Whether no-time-iter was given.  */

} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief the description string of the program */
extern const char *gengetopt_args_info_description;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
