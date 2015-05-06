/***************************************************************************
 *   Copyright (C) 2010 by Daniel Arguelles Martino                        *
 *    daniel.arguelles.martino@gmail.com                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************
 *
 *
 *  \file     dtspg.c
 *  \brief    Linear solvers of Toeplitz matrices
 *  \author   Daniel Arguelles Martino
 *  \date     12/05/11
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

#include "csp_dtspg.h"
/* Command line parser */

#include "cmdline_dtspg.h"
static struct gengetopt_args_info_dtspg args_info;

double dlange_( char *, int *, int *, double *, int *, int * );
/**
 * Parse command line commands into \a args_info structure (uses gengetopt).
 *
 * \param argc Number of params in command line.
 * \param argv Params of command line.
 * \param args_info Structure to be filled with parsed input params.
 */
void parse_cmdline(int argc, char *argv[], struct gengetopt_args_info_dtspg *args_info)
{
	if (cmdline_parser(argc, argv, args_info) != 0) {
		fprintf(stderr, "Run %s --help ('-h') to see the list of options.\n", argv[0]);
		exit(1);
	}

	if (args_info->help_given) {
		cmdline_parser_print_help();
		exit(0);
	}

	if(args_info->version_given) {
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
int read_vector_data(unsigned const int size, double **data, 
		const char *name, unsigned int file_given, const char *file_arg,
		unsigned int random_given, 
		unsigned int save_random_given, const char *save_random_arg, unsigned const int bSize )
{

	unsigned int i = 0;
	unsigned int j = 0;
	/* Select input source */
	if (file_given) {
		/* Read from standard input */
		if (!strcmp(file_arg, "-")) {
			if (!read_input((int) size, data, stdin)) {
				fprintf(stderr, "Error reading %s data from standard input.\n", name);
				return 1;
			}

			/* Read from a specified file */
		}
		else {
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
			double * tmp;
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
			if( bSize == 0 ){
				for (i = 0; i < size; ++i) {
					pos[i] = rand() / (double) RAND_MAX;
					pos[i] *= (rand() & 1) ? (double) 1 : (double) -1;
				}
			}else{

				tmp = (double *) malloc( (size/bSize)*(size/bSize) * sizeof( double ) );
				csp_dtspg_rg( size/(bSize*bSize), bSize, args_info.random_seed_arg, tmp );
				/*memcpy(pos, tmp, size_align);*/

				for( i= 0; i<bSize ; i++ ){
					for( j= 0; j<size/(bSize) ; j++ ){
						/*pos[i] = tmp[i];  */
						pos[j*bSize+i] = tmp[i*size/(bSize)+j];
					}
				}
				free(tmp);

			}

			/* Save randomly generated data */
			if(save_random_given) {
				fprintf(out, "%u ", size);
				for(i=0; i<size; ++i) fprintf(out, "%.8f ", pos[i]);
				fprintf(out, "\n");
			}

			/* Close output file */
			if(save_random_given) fclose(out);
		}
	}

	return 0;
}

int read_input_data(int *n, double **T, double **b, int *bSize, int *nBlocks)
{
    int myN, myBs, myNb;
    
    if( args_info.size_given ){
        myN = args_info.size_arg;
        if( args_info.block_size_given ){
            myBs = args_info.block_size_arg;
            myNb = myN/myBs;
        }else if( args_info.number_blocks_given ){
            myNb = args_info.number_blocks_arg;
            myBs = myN/myNb;
        }else{
            fprintf(stderr, "Unknown problem size.\n");
            return 1;
        }
    }else if( args_info.block_size_given ){
        myBs = args_info.block_size_arg;
        if( args_info.number_blocks_given ){
            myNb = args_info.number_blocks_arg;
            myN = myNb*myBs;
        }else{
            fprintf(stderr, "Unknown problem size.\n");
            return 1;
        }
    }
    
    

	if (myN <= 0) {
		fprintf(stderr, "Invalid problem size.\n");
		return 1;
	}

	/* Initialize random number generator */
	if (args_info.random_seed_given)
		srand(args_info.random_seed_arg);
	else
		srand(time(NULL));

	/* rhs vector */
	args_info.rhs_random_given = !args_info.rhs_file_given;
	if (read_vector_data(myN*3, b, "rhs", args_info.rhs_file_given, args_info.rhs_file_arg, 
				args_info.rhs_random_given, args_info.rhs_save_random_given, 
				args_info.rhs_save_random_arg, 0)) {
		return 1;
	}

	args_info.toeplitz_random_given = !args_info.toeplitz_file_given;
	if (read_vector_data(myN*myBs, T, "toeplitz ", args_info.toeplitz_file_given, args_info.toeplitz_file_arg, 
				args_info.toeplitz_random_given, args_info.toeplitz_save_random_given, 
				args_info.toeplitz_save_random_arg, myBs)) {
		return 1;
	}
	*n       = myN;
	*bSize   = myBs;
    *nBlocks = myNb;

	return 0;
}




int main(int argc, char *argv[]){

	int n,
	    bSize,
        nBlocks;

	int i,j;
	double *T = NULL, *b = NULL;

	double elapsed, ucpu, scpu, gtime;

	/*csp_dtspg_sv( 'T', nBlocks *bSize, bSize, T, b, x);*/


	parse_cmdline(argc, argv, &args_info);

	if (read_input_data(&n, &T, &b, &bSize, &nBlocks)) {
		return 1;
	}

	unsigned int n_align = align_mem(n);
	double *x = (double *) malloc(n_align* 3 * sizeof (double));
	if (args_info.raw_headers_flag) {
		printf("#      n    nBlock  bSize   Time (sec.)   Forward error    Backward error  \n");
		printf("#=========================================================================\n");
	}

	/* Main method */
	gtime = 0;

	ctimer_(&elapsed, &ucpu, &scpu, &gtime);
	csp_dtspg_sv( n, bSize, T, b, x,3);
	ctimer_(&elapsed, &ucpu, &scpu, &gtime);


	if (args_info.show_results_flag)
		print_vector("x", x, n);


	if (args_info.raw_results_flag) {
		double forwd_error = 0.0, backwd_error = 0.0, xnorm1 = 0.0, bnorm1 = 0.0, dnorm1 = 0.0;
		/* Backward error. Using 1-norma */
		double *aux = malloc( n*sizeof(double) );
		int i;

		memcpy( aux, b, n*sizeof(double)); 
		csp_dtspg_gemv( nBlocks, bSize, -1.0, T, x, 1.0, aux);

		for( i = 0; i<n; i++ ){
			dnorm1 += fabs( aux[i] );
			xnorm1 += fabs(x[i]);
			bnorm1 += fabs(b[i]);
		}

		backwd_error = dnorm1 / ( csp_dtspg_nrm1( nBlocks, bSize, T ) * xnorm1 + bnorm1 );


		/* Forward error. Using 1-norma */
		for(i = 0; i<n; i++) x[i] = 1.0;

		csp_dtspg_gemv( nBlocks, bSize, 1.0, T, x, 1.0, aux);
		csp_dtspg_sv( n, bSize, T, aux, x,1);

		for( i = 0; i<n; i++ ) x[i] -= 1.0;

		for(i = 0; i<n; i++) forwd_error += fabs(x[i]);
		forwd_error /= n;
		printf("%8u%8u%8u%12.2e%17.2e%17.2e\n", n, nBlocks, bSize, elapsed, forwd_error, backwd_error);

		free(aux);
	} 

	if (args_info.time_flag & !args_info.raw_results_flag ) {
		printf("Time = %E\n", elapsed);
	}

	/* Free command line parser data */
	cmdline_parser_free(&args_info);

	return 0;








}
