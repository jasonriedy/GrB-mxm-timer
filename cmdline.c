/*
  File autogenerated by gengetopt version 2.23
  generated with the following command:
  gengetopt 

  The developers of gengetopt consider the fixed text that goes in all
  gengetopt output files to be in the public domain:
  we make no copyright claims on it.
*/

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FIX_UNUSED
#define FIX_UNUSED(X) (void) (X) /* avoid warnings for unused params */
#endif

#include <getopt.h>

#include "cmdline.h"

const char *gengetopt_args_info_purpose = "";

const char *gengetopt_args_info_usage = "Usage: GrB-mxm-timer [OPTION]...";

const char *gengetopt_args_info_versiontext = "Copyright 2021, Lucata Corporation";

const char *gengetopt_args_info_description = "Times the iteration of B = A*B";

const char *gengetopt_args_info_help[] = {
  "  -h, --help               Print help and exit",
  "  -V, --version            Print version and exit",
  "  -s, --scale=INT          Scale (log2 # vertices in A)  (default=`16')",
  "  -e, --edgefactor=INT     Edge factor, so # edges = ef * 2^scale\n                             (default=`8')",
  "  -A, --A=FLOAT            R-MAT upper left quadrant probability\n                             (default=`0.55')",
  "  -B, --B=FLOAT            R-MAT upper right & lower left quadrant probability\n                             (default=`0.1')",
  "  -N, --noisefact=FLOAT    Noise factor on each recursion  (default=`0.1')",
  "      --run-powers         Run powers of the generated A matrix rather than\n                             applying A to B  (default=off)",
  "",
  "  -f, --filename=STRING    Filename to read/write for a CSR format",
  "      --dump               Write a file to read  (default=off)",
  "      --binary             File is in binary format  (default=off)",
  "",
  "  -c, --b-ncols=INT        Number of columns in B  (default=`16')",
  "  -C, --b-used-ncols=INT   Number of columns actually used in the initial B\n                             (default=`1')",
  "  -E, --b-nents-col=INT    Number of entries per column in the initial B\n                             (default=`1')",
  "",
  "  -k, --khops=STRING       Number of iterations / hops (can be a space-delim\n                             list)  (default=`2 4 8')",
  "",
  "      --NE-chunk-size=INT  Number of edges to generate in a chunk.\n                             (default=`1048576')",
  "      --verbose[=INT]      Provide status updates via stdout.  (default=`1')",
  "      --no-time-A          Do not time A  (default=off)",
  "      --no-time-B          Do not time B  (default=off)",
  "      --no-time-iter       Do not time iteration  (default=off)",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_STRING
  , ARG_INT
  , ARG_FLOAT
} cmdline_parser_arg_type;

static
void clear_given (struct gengetopt_args_info *args_info);
static
void clear_args (struct gengetopt_args_info *args_info);

static int
cmdline_parser_internal (int argc, char **argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error);


static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->scale_given = 0 ;
  args_info->edgefactor_given = 0 ;
  args_info->A_given = 0 ;
  args_info->B_given = 0 ;
  args_info->noisefact_given = 0 ;
  args_info->run_powers_given = 0 ;
  args_info->filename_given = 0 ;
  args_info->dump_given = 0 ;
  args_info->binary_given = 0 ;
  args_info->b_ncols_given = 0 ;
  args_info->b_used_ncols_given = 0 ;
  args_info->b_nents_col_given = 0 ;
  args_info->khops_given = 0 ;
  args_info->NE_chunk_size_given = 0 ;
  args_info->verbose_given = 0 ;
  args_info->no_time_A_given = 0 ;
  args_info->no_time_B_given = 0 ;
  args_info->no_time_iter_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  FIX_UNUSED (args_info);
  args_info->scale_arg = 16;
  args_info->scale_orig = NULL;
  args_info->edgefactor_arg = 8;
  args_info->edgefactor_orig = NULL;
  args_info->A_arg = 0.55;
  args_info->A_orig = NULL;
  args_info->B_arg = 0.1;
  args_info->B_orig = NULL;
  args_info->noisefact_arg = 0.1;
  args_info->noisefact_orig = NULL;
  args_info->run_powers_flag = 0;
  args_info->filename_arg = NULL;
  args_info->filename_orig = NULL;
  args_info->dump_flag = 0;
  args_info->binary_flag = 0;
  args_info->b_ncols_arg = 16;
  args_info->b_ncols_orig = NULL;
  args_info->b_used_ncols_arg = 1;
  args_info->b_used_ncols_orig = NULL;
  args_info->b_nents_col_arg = 1;
  args_info->b_nents_col_orig = NULL;
  args_info->khops_arg = gengetopt_strdup ("2 4 8");
  args_info->khops_orig = NULL;
  args_info->NE_chunk_size_arg = 1048576;
  args_info->NE_chunk_size_orig = NULL;
  args_info->verbose_arg = 1;
  args_info->verbose_orig = NULL;
  args_info->no_time_A_flag = 0;
  args_info->no_time_B_flag = 0;
  args_info->no_time_iter_flag = 0;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->scale_help = gengetopt_args_info_help[2] ;
  args_info->edgefactor_help = gengetopt_args_info_help[3] ;
  args_info->A_help = gengetopt_args_info_help[4] ;
  args_info->B_help = gengetopt_args_info_help[5] ;
  args_info->noisefact_help = gengetopt_args_info_help[6] ;
  args_info->run_powers_help = gengetopt_args_info_help[7] ;
  args_info->filename_help = gengetopt_args_info_help[9] ;
  args_info->dump_help = gengetopt_args_info_help[10] ;
  args_info->binary_help = gengetopt_args_info_help[11] ;
  args_info->b_ncols_help = gengetopt_args_info_help[13] ;
  args_info->b_used_ncols_help = gengetopt_args_info_help[14] ;
  args_info->b_nents_col_help = gengetopt_args_info_help[15] ;
  args_info->khops_help = gengetopt_args_info_help[17] ;
  args_info->NE_chunk_size_help = gengetopt_args_info_help[19] ;
  args_info->verbose_help = gengetopt_args_info_help[20] ;
  args_info->no_time_A_help = gengetopt_args_info_help[21] ;
  args_info->no_time_B_help = gengetopt_args_info_help[22] ;
  args_info->no_time_iter_help = gengetopt_args_info_help[23] ;
  
}

void
cmdline_parser_print_version (void)
{
  printf ("%s %s\n",
     (strlen(CMDLINE_PARSER_PACKAGE_NAME) ? CMDLINE_PARSER_PACKAGE_NAME : CMDLINE_PARSER_PACKAGE),
     CMDLINE_PARSER_VERSION);

  if (strlen(gengetopt_args_info_versiontext) > 0)
    printf("\n%s\n", gengetopt_args_info_versiontext);
}

static void print_help_common(void)
{
	size_t len_purpose = strlen(gengetopt_args_info_purpose);
	size_t len_usage = strlen(gengetopt_args_info_usage);

	if (len_usage > 0) {
		printf("%s\n", gengetopt_args_info_usage);
	}
	if (len_purpose > 0) {
		printf("%s\n", gengetopt_args_info_purpose);
	}

	if (len_usage || len_purpose) {
		printf("\n");
	}

	if (strlen(gengetopt_args_info_description) > 0) {
		printf("%s\n\n", gengetopt_args_info_description);
	}
}

void
cmdline_parser_print_help (void)
{
  int i = 0;
  print_help_common();
  while (gengetopt_args_info_help[i])
    printf("%s\n", gengetopt_args_info_help[i++]);
}

void
cmdline_parser_init (struct gengetopt_args_info *args_info)
{
  clear_given (args_info);
  clear_args (args_info);
  init_args_info (args_info);
}

void
cmdline_parser_params_init(struct cmdline_parser_params *params)
{
  if (params)
    { 
      params->override = 0;
      params->initialize = 1;
      params->check_required = 1;
      params->check_ambiguity = 0;
      params->print_errors = 1;
    }
}

struct cmdline_parser_params *
cmdline_parser_params_create(void)
{
  struct cmdline_parser_params *params = 
    (struct cmdline_parser_params *)malloc(sizeof(struct cmdline_parser_params));
  cmdline_parser_params_init(params);  
  return params;
}

static void
free_string_field (char **s)
{
  if (*s)
    {
      free (*s);
      *s = 0;
    }
}


static void
cmdline_parser_release (struct gengetopt_args_info *args_info)
{

  free_string_field (&(args_info->scale_orig));
  free_string_field (&(args_info->edgefactor_orig));
  free_string_field (&(args_info->A_orig));
  free_string_field (&(args_info->B_orig));
  free_string_field (&(args_info->noisefact_orig));
  free_string_field (&(args_info->filename_arg));
  free_string_field (&(args_info->filename_orig));
  free_string_field (&(args_info->b_ncols_orig));
  free_string_field (&(args_info->b_used_ncols_orig));
  free_string_field (&(args_info->b_nents_col_orig));
  free_string_field (&(args_info->khops_arg));
  free_string_field (&(args_info->khops_orig));
  free_string_field (&(args_info->NE_chunk_size_orig));
  free_string_field (&(args_info->verbose_orig));
  
  

  clear_given (args_info);
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, const char *values[])
{
  FIX_UNUSED (values);
  if (arg) {
    fprintf(outfile, "%s=\"%s\"\n", opt, arg);
  } else {
    fprintf(outfile, "%s\n", opt);
  }
}


int
cmdline_parser_dump(FILE *outfile, struct gengetopt_args_info *args_info)
{
  int i = 0;

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot dump options to stream\n", CMDLINE_PARSER_PACKAGE);
      return EXIT_FAILURE;
    }

  if (args_info->help_given)
    write_into_file(outfile, "help", 0, 0 );
  if (args_info->version_given)
    write_into_file(outfile, "version", 0, 0 );
  if (args_info->scale_given)
    write_into_file(outfile, "scale", args_info->scale_orig, 0);
  if (args_info->edgefactor_given)
    write_into_file(outfile, "edgefactor", args_info->edgefactor_orig, 0);
  if (args_info->A_given)
    write_into_file(outfile, "A", args_info->A_orig, 0);
  if (args_info->B_given)
    write_into_file(outfile, "B", args_info->B_orig, 0);
  if (args_info->noisefact_given)
    write_into_file(outfile, "noisefact", args_info->noisefact_orig, 0);
  if (args_info->run_powers_given)
    write_into_file(outfile, "run-powers", 0, 0 );
  if (args_info->filename_given)
    write_into_file(outfile, "filename", args_info->filename_orig, 0);
  if (args_info->dump_given)
    write_into_file(outfile, "dump", 0, 0 );
  if (args_info->binary_given)
    write_into_file(outfile, "binary", 0, 0 );
  if (args_info->b_ncols_given)
    write_into_file(outfile, "b-ncols", args_info->b_ncols_orig, 0);
  if (args_info->b_used_ncols_given)
    write_into_file(outfile, "b-used-ncols", args_info->b_used_ncols_orig, 0);
  if (args_info->b_nents_col_given)
    write_into_file(outfile, "b-nents-col", args_info->b_nents_col_orig, 0);
  if (args_info->khops_given)
    write_into_file(outfile, "khops", args_info->khops_orig, 0);
  if (args_info->NE_chunk_size_given)
    write_into_file(outfile, "NE-chunk-size", args_info->NE_chunk_size_orig, 0);
  if (args_info->verbose_given)
    write_into_file(outfile, "verbose", args_info->verbose_orig, 0);
  if (args_info->no_time_A_given)
    write_into_file(outfile, "no-time-A", 0, 0 );
  if (args_info->no_time_B_given)
    write_into_file(outfile, "no-time-B", 0, 0 );
  if (args_info->no_time_iter_given)
    write_into_file(outfile, "no-time-iter", 0, 0 );
  

  i = EXIT_SUCCESS;
  return i;
}

int
cmdline_parser_file_save(const char *filename, struct gengetopt_args_info *args_info)
{
  FILE *outfile;
  int i = 0;

  outfile = fopen(filename, "w");

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot open file for writing: %s\n", CMDLINE_PARSER_PACKAGE, filename);
      return EXIT_FAILURE;
    }

  i = cmdline_parser_dump(outfile, args_info);
  fclose (outfile);

  return i;
}

void
cmdline_parser_free (struct gengetopt_args_info *args_info)
{
  cmdline_parser_release (args_info);
}

/** @brief replacement of strdup, which is not standard */
char *
gengetopt_strdup (const char *s)
{
  char *result = 0;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}

int
cmdline_parser (int argc, char **argv, struct gengetopt_args_info *args_info)
{
  return cmdline_parser2 (argc, argv, args_info, 0, 1, 1);
}

int
cmdline_parser_ext (int argc, char **argv, struct gengetopt_args_info *args_info,
                   struct cmdline_parser_params *params)
{
  int result;
  result = cmdline_parser_internal (argc, argv, args_info, params, 0);

  if (result == EXIT_FAILURE)
    {
      cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
cmdline_parser2 (int argc, char **argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required)
{
  int result;
  struct cmdline_parser_params params;
  
  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;

  result = cmdline_parser_internal (argc, argv, args_info, &params, 0);

  if (result == EXIT_FAILURE)
    {
      cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name)
{
  FIX_UNUSED (args_info);
  FIX_UNUSED (prog_name);
  return EXIT_SUCCESS;
}


static char *package_name = 0;

/**
 * @brief updates an option
 * @param field the generic pointer to the field to update
 * @param orig_field the pointer to the orig field
 * @param field_given the pointer to the number of occurrence of this option
 * @param prev_given the pointer to the number of occurrence already seen
 * @param value the argument for this option (if null no arg was specified)
 * @param possible_values the possible values for this option (if specified)
 * @param default_value the default value (in case the option only accepts fixed values)
 * @param arg_type the type of this option
 * @param check_ambiguity @see cmdline_parser_params.check_ambiguity
 * @param override @see cmdline_parser_params.override
 * @param no_free whether to free a possible previous value
 * @param multiple_option whether this is a multiple option
 * @param long_opt the corresponding long option
 * @param short_opt the corresponding short option (or '-' if none)
 * @param additional_error possible further error specification
 */
static
int update_arg(void *field, char **orig_field,
               unsigned int *field_given, unsigned int *prev_given, 
               char *value, const char *possible_values[],
               const char *default_value,
               cmdline_parser_arg_type arg_type,
               int check_ambiguity, int override,
               int no_free, int multiple_option,
               const char *long_opt, char short_opt,
               const char *additional_error)
{
  char *stop_char = 0;
  const char *val = value;
  int found;
  char **string_field;
  FIX_UNUSED (field);

  stop_char = 0;
  found = 0;

  if (!multiple_option && prev_given && (*prev_given || (check_ambiguity && *field_given)))
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: `--%s' (`-%c') option given more than once%s\n", 
               package_name, long_opt, short_opt,
               (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: `--%s' option given more than once%s\n", 
               package_name, long_opt,
               (additional_error ? additional_error : ""));
      return 1; /* failure */
    }

  FIX_UNUSED (default_value);
    
  if (field_given && *field_given && ! override)
    return 0;
  if (prev_given)
    (*prev_given)++;
  if (field_given)
    (*field_given)++;
  if (possible_values)
    val = possible_values[found];

  switch(arg_type) {
  case ARG_FLAG:
    *((int *)field) = !*((int *)field);
    break;
  case ARG_INT:
    if (val) *((int *)field) = strtol (val, &stop_char, 0);
    break;
  case ARG_FLOAT:
    if (val) *((float *)field) = (float)strtod (val, &stop_char);
    break;
  case ARG_STRING:
    if (val) {
      string_field = (char **)field;
      if (!no_free && *string_field)
        free (*string_field); /* free previous string */
      *string_field = gengetopt_strdup (val);
    }
    break;
  default:
    break;
  };

  /* check numeric conversion */
  switch(arg_type) {
  case ARG_INT:
  case ARG_FLOAT:
    if (val && !(stop_char && *stop_char == '\0')) {
      fprintf(stderr, "%s: invalid numeric value: %s\n", package_name, val);
      return 1; /* failure */
    }
    break;
  default:
    ;
  };

  /* store the original value */
  switch(arg_type) {
  case ARG_NO:
  case ARG_FLAG:
    break;
  default:
    if (value && orig_field) {
      if (no_free) {
        *orig_field = value;
      } else {
        if (*orig_field)
          free (*orig_field); /* free previous string */
        *orig_field = gengetopt_strdup (value);
      }
    }
  };

  return 0; /* OK */
}


int
cmdline_parser_internal (
  int argc, char **argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error)
{
  int c;	/* Character of the parsed option.  */

  int error_occurred = 0;
  struct gengetopt_args_info local_args_info;
  
  int override;
  int initialize;
  int check_required;
  int check_ambiguity;
  
  package_name = argv[0];
  
  /* TODO: Why is this here? It is not used anywhere. */
  override = params->override;
  FIX_UNUSED(override);

  initialize = params->initialize;
  check_required = params->check_required;

  /* TODO: Why is this here? It is not used anywhere. */
  check_ambiguity = params->check_ambiguity;
  FIX_UNUSED(check_ambiguity);

  if (initialize)
    cmdline_parser_init (args_info);

  cmdline_parser_init (&local_args_info);

  optarg = 0;
  optind = 0;
  opterr = params->print_errors;
  optopt = '?';

  while (1)
    {
      int option_index = 0;

      static struct option long_options[] = {
        { "help",	0, NULL, 'h' },
        { "version",	0, NULL, 'V' },
        { "scale",	1, NULL, 's' },
        { "edgefactor",	1, NULL, 'e' },
        { "A",	1, NULL, 'A' },
        { "B",	1, NULL, 'B' },
        { "noisefact",	1, NULL, 'N' },
        { "run-powers",	0, NULL, 0 },
        { "filename",	1, NULL, 'f' },
        { "dump",	0, NULL, 0 },
        { "binary",	0, NULL, 0 },
        { "b-ncols",	1, NULL, 'c' },
        { "b-used-ncols",	1, NULL, 'C' },
        { "b-nents-col",	1, NULL, 'E' },
        { "khops",	1, NULL, 'k' },
        { "NE-chunk-size",	1, NULL, 0 },
        { "verbose",	2, NULL, 0 },
        { "no-time-A",	0, NULL, 0 },
        { "no-time-B",	0, NULL, 0 },
        { "no-time-iter",	0, NULL, 0 },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVs:e:A:B:N:f:c:C:E:k:", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          cmdline_parser_print_help ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'V':	/* Print version and exit.  */
          cmdline_parser_print_version ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 's':	/* Scale (log2 # vertices in A).  */
        
        
          if (update_arg( (void *)&(args_info->scale_arg), 
               &(args_info->scale_orig), &(args_info->scale_given),
              &(local_args_info.scale_given), optarg, 0, "16", ARG_INT,
              check_ambiguity, override, 0, 0,
              "scale", 's',
              additional_error))
            goto failure;
        
          break;
        case 'e':	/* Edge factor, so # edges = ef * 2^scale.  */
        
        
          if (update_arg( (void *)&(args_info->edgefactor_arg), 
               &(args_info->edgefactor_orig), &(args_info->edgefactor_given),
              &(local_args_info.edgefactor_given), optarg, 0, "8", ARG_INT,
              check_ambiguity, override, 0, 0,
              "edgefactor", 'e',
              additional_error))
            goto failure;
        
          break;
        case 'A':	/* R-MAT upper left quadrant probability.  */
        
        
          if (update_arg( (void *)&(args_info->A_arg), 
               &(args_info->A_orig), &(args_info->A_given),
              &(local_args_info.A_given), optarg, 0, "0.55", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "A", 'A',
              additional_error))
            goto failure;
        
          break;
        case 'B':	/* R-MAT upper right & lower left quadrant probability.  */
        
        
          if (update_arg( (void *)&(args_info->B_arg), 
               &(args_info->B_orig), &(args_info->B_given),
              &(local_args_info.B_given), optarg, 0, "0.1", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "B", 'B',
              additional_error))
            goto failure;
        
          break;
        case 'N':	/* Noise factor on each recursion.  */
        
        
          if (update_arg( (void *)&(args_info->noisefact_arg), 
               &(args_info->noisefact_orig), &(args_info->noisefact_given),
              &(local_args_info.noisefact_given), optarg, 0, "0.1", ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "noisefact", 'N',
              additional_error))
            goto failure;
        
          break;
        case 'f':	/* Filename to read/write for a CSR format.  */
        
        
          if (update_arg( (void *)&(args_info->filename_arg), 
               &(args_info->filename_orig), &(args_info->filename_given),
              &(local_args_info.filename_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "filename", 'f',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* Number of columns in B.  */
        
        
          if (update_arg( (void *)&(args_info->b_ncols_arg), 
               &(args_info->b_ncols_orig), &(args_info->b_ncols_given),
              &(local_args_info.b_ncols_given), optarg, 0, "16", ARG_INT,
              check_ambiguity, override, 0, 0,
              "b-ncols", 'c',
              additional_error))
            goto failure;
        
          break;
        case 'C':	/* Number of columns actually used in the initial B.  */
        
        
          if (update_arg( (void *)&(args_info->b_used_ncols_arg), 
               &(args_info->b_used_ncols_orig), &(args_info->b_used_ncols_given),
              &(local_args_info.b_used_ncols_given), optarg, 0, "1", ARG_INT,
              check_ambiguity, override, 0, 0,
              "b-used-ncols", 'C',
              additional_error))
            goto failure;
        
          break;
        case 'E':	/* Number of entries per column in the initial B.  */
        
        
          if (update_arg( (void *)&(args_info->b_nents_col_arg), 
               &(args_info->b_nents_col_orig), &(args_info->b_nents_col_given),
              &(local_args_info.b_nents_col_given), optarg, 0, "1", ARG_INT,
              check_ambiguity, override, 0, 0,
              "b-nents-col", 'E',
              additional_error))
            goto failure;
        
          break;
        case 'k':	/* Number of iterations / hops (can be a space-delim list).  */
        
        
          if (update_arg( (void *)&(args_info->khops_arg), 
               &(args_info->khops_orig), &(args_info->khops_given),
              &(local_args_info.khops_given), optarg, 0, "2 4 8", ARG_STRING,
              check_ambiguity, override, 0, 0,
              "khops", 'k',
              additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
          /* Run powers of the generated A matrix rather than applying A to B.  */
          if (strcmp (long_options[option_index].name, "run-powers") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->run_powers_flag), 0, &(args_info->run_powers_given),
                &(local_args_info.run_powers_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "run-powers", '-',
                additional_error))
              goto failure;
          
          }
          /* Write a file to read.  */
          else if (strcmp (long_options[option_index].name, "dump") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->dump_flag), 0, &(args_info->dump_given),
                &(local_args_info.dump_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "dump", '-',
                additional_error))
              goto failure;
          
          }
          /* File is in binary format.  */
          else if (strcmp (long_options[option_index].name, "binary") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->binary_flag), 0, &(args_info->binary_given),
                &(local_args_info.binary_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "binary", '-',
                additional_error))
              goto failure;
          
          }
          /* Number of edges to generate in a chunk..  */
          else if (strcmp (long_options[option_index].name, "NE-chunk-size") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->NE_chunk_size_arg), 
                 &(args_info->NE_chunk_size_orig), &(args_info->NE_chunk_size_given),
                &(local_args_info.NE_chunk_size_given), optarg, 0, "1048576", ARG_INT,
                check_ambiguity, override, 0, 0,
                "NE-chunk-size", '-',
                additional_error))
              goto failure;
          
          }
          /* Provide status updates via stdout..  */
          else if (strcmp (long_options[option_index].name, "verbose") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->verbose_arg), 
                 &(args_info->verbose_orig), &(args_info->verbose_given),
                &(local_args_info.verbose_given), optarg, 0, "1", ARG_INT,
                check_ambiguity, override, 0, 0,
                "verbose", '-',
                additional_error))
              goto failure;
          
          }
          /* Do not time A.  */
          else if (strcmp (long_options[option_index].name, "no-time-A") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->no_time_A_flag), 0, &(args_info->no_time_A_given),
                &(local_args_info.no_time_A_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "no-time-A", '-',
                additional_error))
              goto failure;
          
          }
          /* Do not time B.  */
          else if (strcmp (long_options[option_index].name, "no-time-B") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->no_time_B_flag), 0, &(args_info->no_time_B_given),
                &(local_args_info.no_time_B_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "no-time-B", '-',
                additional_error))
              goto failure;
          
          }
          /* Do not time iteration.  */
          else if (strcmp (long_options[option_index].name, "no-time-iter") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->no_time_iter_flag), 0, &(args_info->no_time_iter_given),
                &(local_args_info.no_time_iter_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "no-time-iter", '-',
                additional_error))
              goto failure;
          
          }
          
          break;
        case '?':	/* Invalid option.  */
          /* `getopt_long' already printed an error message.  */
          goto failure;

        default:	/* bug: option not considered.  */
          fprintf (stderr, "%s: option unknown: %c%s\n", CMDLINE_PARSER_PACKAGE, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */



	FIX_UNUSED(check_required);

  cmdline_parser_release (&local_args_info);

  if ( error_occurred )
    return (EXIT_FAILURE);

  return 0;

failure:
  
  cmdline_parser_release (&local_args_info);
  return (EXIT_FAILURE);
}
/* vim: set ft=c noet ts=8 sts=8 sw=8 tw=80 nojs spell : */
