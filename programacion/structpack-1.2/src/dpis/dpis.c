/*
 *  \file     dpis.c
 *  \brief    Linear solvers of Symmetric Tridiagonal Toeplitz matrices
 *  \author   Pedro Alonso Jordá, Antonio M. Vidal Maciá, Pablo Martinez Naredo
 *  \date     30/05/10
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include <cutils.h>

/* BLAS */
#if USE_MKL
#include <mkl_blas.h>
#else
#include <blas.h>
#endif

/* csp_dpis module header */
#include "csp_dpis.h"

/* Command line parser */
#include "cmdline_dpis.h"
static struct gengetopt_args_info_dpis args_info;

/*
 * Parse command line commands into \a args_info structure (uses gengetopt).
 *
 * \param argc Number of params in command line.
 * \param argv Params of command line.
 * \param args_info Structure to be filled with parsed input params.
*/
void parse_cmdline(int argc, char *argv[], struct gengetopt_args_info_dpis *args_info)
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

/*
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

      /* Generate multiple vectors */
      double *pos = *data;
      
      /* Fill input data randomly */
      for (i = 0; i < size; ++i) {
	pos[i] = rand() / (double) RAND_MAX;
	pos[i] *= (rand() & 1) ? (double) 1 : (double) -1;
      }
      
      /* Save randomly generated data */
      if (save_random_given) {
	fprintf(out, "%u ", size);
	for (i=0; i<size; ++i)
	  fprintf(out, "%.8f ", pos[i]);
	fprintf(out, "\n");

        /* Close output file */
        fclose(out);
      }
    }
  }

  return 0;
}

/*
 * Read Toeplitz and rhs vector input data.
 *
 * \param n  Pointer to problem size.
 * \param t0 Pointer to diagonal entry.
 * \param t1 Pointer to off-diagonal entries.
 * \param b  Pointer to rhs data.
 * \return \c 0 if successful, non-zero otherwise.
*/
int read_input_data(int *n, double *t0, double *t1, double **b)
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

  /* Toeplitz matrix */
  if (args_info.diagonal_given) {
    *(t0) = args_info.diagonal_arg;
    *(t1) = args_info.off_diagonal_arg;
  }
  else {
    /* Generate input data for Toeplitz matrix */
    *(t0) = rand() / (double) RAND_MAX;
    *(t1) = rand() / (double) RAND_MAX;
  }

  /* rhs vector */
  args_info.rhs_random_given = !args_info.rhs_file_given;
  if (read_vector_data(*n, b, "rhs", args_info.rhs_file_given, args_info.rhs_file_arg, 
		       args_info.rhs_random_given, args_info.rhs_save_random_given, 
		       args_info.rhs_save_random_arg)) {
    return 1;
  }

  return 0;
}


int solver_rojo(unsigned int n, double t0, double t1, double* b)
{
  csp_dpis_rojosv(n, t0, t1, b);

  return 0;
}

int solver_dst_method(unsigned int n, double t0, double t1, double* b)
{
  csp_dpis_dstsv(n, t0, t1, b);

  return 0;
}

int solver_ldlt_method(unsigned int n, double t0, double t1, double* b)
{
  return csp_dpis_ldltsv(n, t0, t1, b);
}


/* Entry point */
int main(int argc, char *argv[])
{
  int n=0, inc1=1, rojo_set=0, dst_method_set=0, ldlt_method_set=0;
  double *b=NULL, t0, t1, MINUSONE=-1.0;
  double gtime, elapsed, ucpu, scpu;    

  /* Parse command line input */
  parse_cmdline(argc, argv, &args_info);

  if (read_input_data(&n, &t0, &t1, &b)) {
    return 1;
  }

  unsigned int n_align = align_mem(n);
  double *x = (double *) malloc(n_align * sizeof (double));

  /* Method */
  rojo_set        = args_info.rojo_given;
  dst_method_set  = args_info.dst_given;
  ldlt_method_set = args_info.ldlt_given;
  rojo_set        = rojo_set || (!dst_method_set && !ldlt_method_set);

  if (args_info.raw_headers_flag) {
     printf("#      n     Time (sec.)   Forward error    Backward error  Max. svd    Min. svd        Method\n");
     printf("#===================================================================================================\n");
  }

  gtime=0.0;

  /* Rojo method */
  if (rojo_set) {
      memcpy(x, b, n_align * sizeof (double));
      ctimer_(&elapsed, &ucpu, &scpu, &gtime);
        if (solver_rojo(n, t0, t1, x)) {
          fprintf(stderr, "Problem with the solver.\n");
          return 1;
        }
      ctimer_(&elapsed, &ucpu, &scpu, &gtime);

      /* Show results */
      if (args_info.show_results_flag)
	print_vector("x", x, n);

      if (args_info.raw_results_flag) {
        /* Backward error */
        double *r = (double *) malloc(n_align * sizeof (double));
        memcpy(r, b, n_align * sizeof (double));
        csp_dpis_gemv(n, -1.0, t0, t1, x, r);
        
        double backwd_error = dnrm2_(&n, r, &inc1);
        double maxsvd       = csp_dpis_maxsvd(n, t0, t1);
        backwd_error       /= (maxsvd * dnrm2(&n, x, &inc1))+ dnrm2(&n, b, &inc1);

        /* Forward error */
        int i;
        for (i = 0; i < n; i++) x[i] = 0.0;
        for (i = 0; i < n; i++) r[i] = 0.0;
        x[0] = 1.0;
        csp_dpis_gemv(n, 1.0, t0, t1, x, r);
        solver_rojo      (n, t0, t1, r);
        daxpy_           (&n, &MINUSONE, x, &inc1, r, &inc1);
        
        double forwd_error = dnrm2_(&n, r, &inc1);
        forwd_error       /= dnrm2_(&n, x, &inc1);
        double minsvd      = csp_dpis_minsvd(n, t0, t1);
        printf("%8u%14.2E%14.2E%17.2E%15.2E%12.2E\t(rojo)\n", n, elapsed,
	       forwd_error, backwd_error, maxsvd, minsvd);
        free(r);
      } 
  }

  /* DST method */
  if (dst_method_set) {
      memcpy( x, b, n_align * sizeof(double) );
      ctimer_(&elapsed, &ucpu, &scpu, &gtime);
        if(solver_dst_method( n, t0, t1, x )) {
          fprintf(stderr, "Problem with the solver.\n");
          return 1;
        }
      ctimer_(&elapsed, &ucpu, &scpu, &gtime);

      /* Show results */
      if (args_info.show_results_flag)
        print_vector("x", x, n);

      if( args_info.raw_results_flag ) {
        /* Backward error */
        double* r = (double *) malloc(n_align * sizeof(double));
        memcpy(r, b, n_align * sizeof(double));
        csp_dpis_gemv(n, -1.0, t0, t1, x, r);
        
        double backwd_error = dnrm2_(&n, r, &inc1);
        double maxsvd       = csp_dpis_maxsvd(n,t0,t1);
        backwd_error       /= (maxsvd * dnrm2(&n, x, &inc1)) + dnrm2(&n, b, &inc1);

        /* Forward error */
        int i;
        for( i=0; i<n; i++ ) x[i] = 0.0;
        for( i=0; i<n; i++ ) r[i] = 0.0;
        x[0] = 1.0;
        csp_dpis_gemv(n, 1.0, t0, t1, x, r);
        solver_dst_method(n, t0, t1, r);
        daxpy_           (&n, &MINUSONE, x, &inc1, r, &inc1);

        double forwd_error = dnrm2_(&n, r, &inc1);
        forwd_error       /= dnrm2_(&n, x, &inc1);
        double minsvd      = csp_dpis_minsvd(n,t0,t1);
        printf("%8u%14.2E%14.2E%17.2E%15.2E%12.2E\t(dst_method)\n",n,elapsed,forwd_error,backwd_error,maxsvd,minsvd);
        free(r);
      } 
  }


  /* LDLT method */
  if (ldlt_method_set) {
      memcpy(x, b, n_align * sizeof (double));
      ctimer_(&elapsed, &ucpu, &scpu, &gtime);
        if (solver_ldlt_method(n, t0, t1, x)) {
          fprintf(stderr,"      --> Impossible to solve the problem with this method.\t\t\t\t(ldlt_method)\n");
        } 
      ctimer_(&elapsed, &ucpu, &scpu, &gtime);

      /* Show results */
      if (args_info.show_results_flag)
        print_vector("x", x, n);

      if (args_info.raw_results_flag) {
        /* Backward error */
        double *r = (double *)malloc(n_align * sizeof(double));
        memcpy(r, b, n_align * sizeof(double));
        csp_dpis_gemv(n, -1.0, t0, t1, x, r);
        
        double backwd_error = dnrm2_(&n, r, &inc1);
        double maxsvd       = csp_dpis_maxsvd(n,t0,t1);
        backwd_error       /= (maxsvd * dnrm2_(&n, x, &inc1)) + dnrm2_(&n, b, &inc1);

        /* Forward error */
        int i;
        for( i=0; i<n; i++ ) x[i] = 1.0;
        for( i=0; i<n; i++ ) r[i] = 0.0;
        csp_dpis_gemv (n, 1.0, t0, t1, x, r);
        solver_ldlt_method(n, t0, t1, r);
        daxpy_            (&n, &MINUSONE, x, &inc1, r, &inc1);
        
        double forwd_error  = dnrm2_(&n, r, &inc1);
        forwd_error        /= dnrm2_(&n, x, &inc1);
        double minsvd       = csp_dpis_minsvd(n,t0,t1);
        printf("%8u%14.2E%14.2E%17.2E%15.2E%12.2E\t(ldlt_method)\n",n,elapsed,forwd_error,backwd_error,maxsvd,minsvd);
        free(r);
      } 
  }

  if (args_info.time_flag && !args_info.raw_results_flag) {
    printf("Time=%E\n", elapsed);
  }

  /* Free command line parser data */
  cmdline_parser_free(&args_info);

  return 0;
}
