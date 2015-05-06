/*
 * \file     dts.c
 * \brief    Linear solvers of Symmetric Toeplitz matrices
 * \author   Pablo Martínez Naredo and Pedro Alonso Jordá
 * \date     15/01/11
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include <cutils.h>

/* BLAS */
#if USE_MKL
#include <mkl_blas.h>
#else
#include <blas.h>
#endif

/* csp_dts module */
#include "csp_dts.h"

/* Command line parser */
#include "cmdline_dts.h"
static struct gengetopt_args_info_dts args_info;


/**
 * Parse command line commands into \a args_info structure (uses gengetopt).
 *
 * \param argc Number of params in command line.
 * \param argv Params of command line.
 * \param args_info Structure to be filled with parsed input params.
*/
void parse_cmdline(int argc, char *argv[], struct gengetopt_args_info_dts *args_info)
{
  if (cmdline_parser(argc, argv, args_info) != 0) {
    fprintf(stderr, "Run %s --help ('-h') to see the list of options.\n", argv[0]);
    exit(1);
  }
  
  if (args_info->help_given) {
    cmdline_parser_print_help();
    exit(0);
  }

  if (args_info->version_given) {
    cmdline_parser_print_version();
    exit(0);
  }
}

/**
 * Read input data to a vector (including memory allocation).
 *
 * \param size Vector size.
 * \param data Pointer to vector data.
 * \param name Name of the vector being read (used in error reporting only).
 * \param file_given Command line flag indicating if input should be stdin or a file.
 * \param file_arg Command line parameter indicating input data file.
 * \param random_given Command line flag indicating if data should be randomly-generated between -1 and 1.
 * \param save_random_given Command line flag indicating if randomly-generated data should be saved.
 * \param save_random_arg Command line parameter indicating where to save generated data.
 * \return \c 0 if successful, non-zero otherwise.
*/
int read_vector_data(unsigned int size, double **data, 
		     const char *name, unsigned int file_given, const char *file_arg,
		     unsigned int random_given, 
		     unsigned int save_random_given, const char *save_random_arg)
{
  unsigned int i = 0;

  /* Select input source */
  if (file_given) {
    if (!strcmp(file_arg, "-")) {
      /* Read from standard input */
      if (!read_input((int) size, data, stdin)) {
        fprintf(stderr, "Error reading %s data from standard input.\n", name);
        return 1;
      }
    }
    else {
      /* Read from a specified file */
      FILE *input = fopen(file_arg, "r");
      if (!input) {
        fprintf(stderr, "Error opening %s input file \"%s\"\n", name, file_arg);
        return 1;
      }

      if (!read_input((int) size, data, input)) {
        fprintf(stderr, "Error reading %s input file \"%s\"\n", name, file_arg);
        fclose(input);
        return 1;
      }

      fclose(input);
    }
  }

  /* Generate input data randomly */
  else {
    if (random_given) {
      /* Allocate memory */
      unsigned int size_align = align_mem(size);
      *data = (double *) malloc(size_align * sizeof (double));

      /* Open output file if vectors need to be saved */
      FILE *out = NULL;
      if (save_random_given) {
	if (!(out = fopen(save_random_arg, "w"))) {
	  fprintf(stderr, "Error opening %s output file \"%s\"\n", name, save_random_arg);
	  free(*data);
	  return 1;
	}
      }
      
      /* Fill input data randomly */
      double *pos = *data;
      for (i = 0; i < size; ++i) {
	/* pos[i] = rand() / (double) RAND_MAX; */
	pos[i] = 2.0 * ( rand() / (double) RAND_MAX ) - 1.0;
	/* pos[i] *= (rand() & 1) ? (double) 1 : (double) -1; */
	/* printf("%lf\n",(rand() & 1) ? (double) 1 : (double) -1); */
        /* printf("t[%d] = %lf\n",i,pos[i]); */
      }
      
      /* Save randomly generated data */
      if (save_random_given) {
	fprintf(out, "%u ", size);
	for(i=0; i<size; ++i) fprintf(out, "%.8f ", pos[i]);
	fprintf(out, "\n");

	/* Close output file */
	fclose(out);
      }
    }
  }
  return 0;
}

/**
 * Read Toeplitz and rhs vector input data.
 *
 * \param n  Pointer to problem size.
 * \param t  Pointer to the first column.
 * \param b  Pointer to rhs data.
 * \param nb Pointer to the block size.
 * \return \c 0 if successful, non-zero otherwise.
*/
int read_input_data(int *n, double **t, double **b, int *nb)
{
  if (args_info.size_arg <= 0) {
    fprintf(stderr, "Invalid problem size.\n");
    return 1;
  }
  *n = args_info.size_arg;

  /* Initialize random number generator */
  if (args_info.random_seed_given)
    srand(args_info.random_seed_arg);
  else
    srand(time(NULL));

  /* Toeplitz first column matrix */
  args_info.toeplitz_random_given = !args_info.first_column_file_given;
  if (read_vector_data(*n, t, "toeplitz", args_info.first_column_file_given, args_info.first_column_file_arg, 
		       args_info.toeplitz_random_given, args_info.first_column_save_random_given, 
		       args_info.first_column_save_random_arg)) {
    return 1;
  }

  /* rhs vector */
  args_info.rhs_random_given = !args_info.rhs_file_given;
  if (read_vector_data(*n, b, "rhs", args_info.rhs_file_given, args_info.rhs_file_arg, 
		       args_info.rhs_random_given, args_info.rhs_save_random_given, 
		       args_info.rhs_save_random_arg)) {
    return 1;
  }

  /* Block size */
  if (args_info.block_size_given) {
    *(nb) = args_info.block_size_arg;
  }
  else {
    /* Block size is 64 by default */
    *(nb) = 64;
  }

  return 0;
}


/* Entry point */
int main(int argc, char *argv[])
{
  int n=0, nb;
  double elapsed, ucpu, scpu, gtime;
  double *t=NULL, *b=NULL;

  /* Parse command line */
  parse_cmdline(argc, argv, &args_info);

  if (read_input_data(&n, &t, &b, &nb)) {
    return 1;
  }

  /* Check arguments */

  unsigned int n_align = align_mem(n);
  double *x = (double *) malloc(n_align * sizeof (double));

  if (args_info.raw_headers_flag) {
    printf("#      n     Time (sec.)   Forward error    Backward error  \n");
    printf("#===========================================================\n");
  }

  memcpy(x, b, n_align * sizeof (double));
  gtime=0.0;
  ctimer_(&elapsed, &ucpu, &scpu, &gtime);
  csp_dts_sv(n, t, x, b, nb, args_info.pivoting_flag, args_info.core_affinity_flag );
  ctimer_(&elapsed, &ucpu, &scpu, &gtime);

  /* Show results */
  if (args_info.show_results_flag)
    print_vector("x", x, n);

  if (args_info.raw_results_flag) {
    /* Backward error using 1-norm */
    double *r=NULL, backwd_error=0.0, forwd_error=0.0, norm=0.0, n1x=0.0, n1b=0.0;
    int i;
    
    r = (double *) malloc(n_align * sizeof (double));
    memcpy(r, b, n_align * sizeof (double));
    csp_dts_gemv(n, -1.0, t, x, r);

    norm = csp_dts_nrm1(n, t);
    for (i = 0; i<n; i++) {
      backwd_error += fabs(r[i]);
      n1x          += fabs(x[i]);
      n1b          += fabs(b[i]);
    }
    backwd_error /= (norm * n1x) + n1b;

    /* Forward error using 1-norm */
    for (i = 0; i<n; i++)
      x[i] = 1.0;
    
    csp_dts_gemv  (n, 1.0, t, x, r);
    csp_dts_sv(n, t, x, r, nb, args_info.pivoting_flag, args_info.core_affinity_flag );

    for (i = 0; i<n; i++)
      x[i] -= 1.0;

    for (i = 0; i<n; i++)
      forwd_error += fabs(x[i]);

    forwd_error /= n;
    printf("%8u%14.2E%15.2E%18.2E\n", n, elapsed, forwd_error, backwd_error);
    free(r);
  } 

  if (args_info.time_flag & !args_info.raw_results_flag ) {
    printf("Time = %10.4lf sec.\n", elapsed);
  }

  /* Free command line parser data */
  cmdline_parser_free(&args_info);

  return 0;
}
