/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff.c                                                  **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauß-Seidel and    **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>

#include "partdiff.h"

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time; /* time when program started                      */
struct timeval comp_time;  /* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
	uint64_t i;

	for (i = 0; i < arguments->num_matrices; i++)
	{
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
		exit(1);
	}

	return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

		for (j = 0; j <= N; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{
			for (i = 0; i <= N; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * i);
				Matrix[g][i][N] = h * i;
				Matrix[g][0][i] = 1.0 - (h * i);
				Matrix[g][N][i] = h * i;
			}

			Matrix[g][N][0] = 0.0;
			Matrix[g][0][N] = 0.0;
		}
	}
}

/* ************************************************************************ */
/*                                                                          */
/* @brief Solves the equation using the specified method and options.       */
/*                                                                          */
/* This function performs the calculation to solve                          */
/* the equation based on the provided arguments,                            */
/* results, and options. It supports different methods                      */
/* (e.g., Jacobi) and termination criteria                                  */
/* (e.g., precision, iteration count).                                      */
/*                                                                          */
/* @param arguments Pointer to the structure                                */
/*      containing calculation arguments.                                   */
/* @param results Pointer to the structure                                  */
/*      to store calculation results.                                       */
/* @param options Pointer to the structure                                  */
/*      containing options for the calculation.                             */
/* ************************************************************************ */
static
void
calculate_seq (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	int i, j;           /* local variables for loops */
	int m1, m2;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

    /* Calculate the pih and fpisin values if the function is FPISIN */
	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

    /* Perform the calculation based on the termination method:*/
	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxResiduum = 0;

		/* over all rows */
		for (i = 1; i < N; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}

		results->stat_iteration++;
		results->stat_precision = maxResiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxResiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	results->m = m2;
}

/* ************************************************************************ */
/* @brief solves the pde for the complete matrix useing openmp in a row vise*/
/*  fashion.                                                                */
/*                                                                          */
/* This function performs the calculation on the matrix based               */
/* on the provided arguments, results, and options. It supports parallel    */
/* execution using OpenMP if the Jacobi method is selected.                 */
/*                                                                          */
/* @param arguments Pointer to the structure containing                     */
/*  calculation arguments.                                                  */
/* @param results Pointer to the structure containing                       */
/*  calculation results.                                                    */
/* @param options Pointer to the structure containing                       */
/*  options for the calculation.                                            */
/*                                                                          */
/* The function uses the following parameters from the structures:          */
/* - arguments->N: The size of the matrix.                                  */
/* - arguments->h: The step size.                                           */
/* - arguments->Matrix: The matrices to be calculated on.                   */
/* - options->method: The method to be used (Jacobi or Gauss-Seidel).       */
/* - options->inf_func: The function to be used for initialization.         */
/* - options->term_iteration: The number of iterations for the calculation. */
/* - options->termination: The termination condition (precision/iteration). */
/* - options->term_precision: The precision for termination.                */
/* - options->number: The number of threads for parallel execution.         */
/* - results->stat_iteration: The number of iterations performed.           */
/* - results->stat_precision: The precision achieved.                       */
/* - results->m: The index of the matrix used for the final result.         */
/*                                                                          */
/* The function performs the following steps:                               */
/* 1. Initialize local variables and parameters based on the options.       */
/* 2. Set up parallel execution if the Jacobi method's selected.            */
/* 3. Perform the calculation for each row and column of the matrix.
/*    Each row gets assinged a thread and each thread works on all colloms. */
/* 4. Update the maximum residuum value for convergence checking.           */
/* 5. Exchange the old and new matrices for the next iteration.             */
/* 6. Check for stopping conditions based on the termination method.        */
/* 7. Update the results structure with iteration count and precision.      */
/* ************************************************************************ */
static
void
calculate_row (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	int i, j;           /* local variables for loops */
	int m1, m2;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */
	double maxLocalResiduum; /* maximum residuum value of a slave in iteration in a thread */
										//maxLocalResiduum updated for each calculation, doesnt require sync. between Threads
        								//code is faster
	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	double** Matrix_Out;
	double** Matrix_In;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}					//wenn if true wird parallelisiert            //threads über ersten parameter gesetzt
#pragma omp parallel if(options->method == METH_JACOBI) default(none) num_threads(options->number) 		/* private oder shared mit default deklaration wird erzwungen */ \
    private(i, j, star, residuum, maxLocalResiduum) 		/*nur relevant für eine Rechung im jeweiligen Thread, werden verändert */ \
    shared(arguments, m1, m2, N, options, pih, fpisin, term_iteration, results, Matrix_Out, Matrix_In, maxResiduum)
        //werden nicht verändert oder nur syncronisiert verändert
{
	while (term_iteration > 0)
	{
        #pragma omp master				//wird nur im Master thread ausgeführt
        {
		 	Matrix_Out  = arguments->Matrix[m1];
			Matrix_In   = arguments->Matrix[m2];
		    maxResiduum = 0;
		}
		maxLocalResiduum = 0;
        #pragma omp barrier			//sync punkt, alle threads müssen diesen erreichen bevor weiter gerechnet wird
									//außerhalb von schleifen, zwecks performance und Threadunabhänigkeit
		/* over all rows */			//Reihen hier paralellisiert
        #pragma omp for				//for schleife über alle threads parallel , ++++++++ Aufgabe 2b scheduler(x,y)++++++++++++
		for (i = 1; i < N; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
                    if (residuum > maxLocalResiduum) {maxLocalResiduum = residuum;}			//maxLocal verwendet
				}

				Matrix_Out[i][j] = star;
			}
		}
		#pragma omp atomic compare												//compiler stellt korrektheit vom nächsten If statement sicher
		if (maxLocalResiduum > maxResiduum) {maxResiduum = maxLocalResiduum;}
        #pragma omp barrier //sync
		#pragma omp master
		{

			results->stat_iteration++;
			results->stat_precision = maxResiduum;

			/* exchange m1 and m2 */
			i = m1;
			m1 = m2;
			m2 = i;

			/* check for stopping calculation depending on termination method */
			if (options->termination == TERM_PREC)
			{
				if (maxResiduum < options->term_precision)
				{
					term_iteration = 0;
				}
			}
			else if (options->termination == TERM_ITER)
			{
				term_iteration--;
			}
        }
        #pragma omp barrier //sync
	}
} /* End #pragma parallel */

	results->m = m2;
}

/* ************************************************************************ */
/* @brief Calculate the pde of the matrix in a collumn vise fashion         */
/*                                                                          */
/* This function performs the calculation on the matrix based               */
/* on the provided arguments, results, and options. It supports parallel    */
/* execution using OpenMP if the Jacobi method is selected.                 */
/*                                                                          */
/* @param arguments Pointer to the structure containing                     */
/*  calculation arguments.                                                  */
/* @param results Pointer to the structure containing                       */
/*  calculation results.                                                    */
/* @param options Pointer to the structure containing                       */
/*  options for the calculation.                                            */
/*                                                                          */
/* The function uses the following parameters from the structures:          */
/* - arguments->N: The size of the matrix.                                  */
/* - arguments->h: The step size.                                           */
/* - arguments->Matrix: The matrix to be calculated.                        */
/* - options->method: The method to be used (Jacobi or Gauss-Seidel).       */
/* - options->inf_func: The function to be used for initialization.         */
/* - options->term_iteration: The number of iterations for the calculation. */
/* - options->termination: The termination condition (precision/iteration). */
/* - options->term_precision: The precision for termination.                */
/* - options->number: The number of threads for parallel execution.         */
/* - results->stat_iteration: The number of iterations performed.           */
/* - results->stat_precision: The precision achieved.                       */
/* - results->m: The index of the matrix used for the final result.         */
/*                                                                          */
/* The function performs the following steps:                               */
/* 1. Initialize local variables and parameters based on the options.       */
/* 2. Set up parallel execution if the Jacobi method's selected.            */
/* 3. Perform the calculation for each row and column of the matrix.        */
/* 4. Update the maximum residuum value for convergence checking.           */
/* 5. Exchange the old and new matrices for the next iteration.             */
/* 6. Check for stopping conditions based on the termination method.        */
/* 7. Update the results structure with iteration count and precision.      */
/* ************************************************************************ */
static
void
calculate_column (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	int i, j;           /* local variables for loops */
	int m1, m2;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */
	double maxLocalResiduum; /* maximum residuum value of a slave in iteration in a thread */
    
	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;
	double fpisin_i = 0.0;

	double** Matrix_Out;
	double** Matrix_In;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}
#pragma omp parallel if(options->method == METH_JACOBI) default(none) num_threads(options->number) \
    private(i, j, star, residuum, maxLocalResiduum) 	/*fpisin_i vorher in private, jetzt in shared. Gleich über alle threads, weil reihe für reihe berechnet wird */ \
    shared(arguments, m1, m2, N, options, pih, fpisin, term_iteration, results, Matrix_Out, Matrix_In, maxResiduum, fpisin_i)
{
	while (term_iteration > 0)
	{
        /* Master thread initializes the matrices for the current iteration and resets maxResiduum */
        #pragma omp master
        {
		 	Matrix_Out  = arguments->Matrix[m1];
			Matrix_In   = arguments->Matrix[m2];
		    maxResiduum = 0;
		}
		maxLocalResiduum = 0;
        #pragma omp barrier

		/* over all rows */
		for (i = 1; i < N; i++)
		{
            #pragma omp single										//block wird nur nch einmal ausgeführt
            {
				fpisin_i = 0.0;

				if (options->inf_func == FUNC_FPISIN)
				{
					fpisin_i = fpisin * sin(pih * (double)i);
				}
            }

			/* over all columns */					//spalten hier parallelisiert
            #pragma omp for
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
                    if (residuum > maxLocalResiduum) {maxLocalResiduum = residuum;}
				}

				Matrix_Out[i][j] = star;
			}
		}
        /*thread calculates the maxLocalresiduum value */
		#pragma omp atomic compare
		if (maxLocalResiduum > maxResiduum) {maxResiduum = maxLocalResiduum;}
        #pragma omp barrier
		#pragma omp master
		{

			results->stat_iteration++;
			results->stat_precision = maxResiduum;

			/* exchange m1 and m2 */
			i = m1;
			m1 = m2;
			m2 = i;

			/* check for stopping calculation depending on termination method */
			if (options->termination == TERM_PREC)
			{
                /* If the maximum residuum is smaller than the precision, the calculation stops */
				if (maxResiduum < options->term_precision)
				{
					term_iteration = 0;
				}
			}
			else if (options->termination == TERM_ITER)
			{
				term_iteration--;
			}
        }
        #pragma omp barrier
	}
} /* End #pragma parallel */

	results->m = m2;
}

/* ************************************************************************ */
/* @brief Calculate an element of the matrix using                          */
/*  the specified method and options.                                       */
/*                                                                          */
/* This function performs the calculation of an element in the matrix based */
/* on the provided arguments, results, and options. It supports parallel    */
/* execution using OpenMP if the Jacobi method is selected.                 */
/*                                                                          */
/* @param arguments Pointer to the structure containing                     */
/*  calculation arguments.                                                  */
/* @param results Pointer to the structure containing                       */
/*  calculation results.                                                    */
/* @param options Pointer to the structure containing                       */
/*  options for the calculation.                                            */
/*                                                                          */
/* The function uses the following parameters from the structures:          */
/* - arguments->N: The size of the matrix.                                  */
/* - arguments->h: The step size.                                           */
/* - arguments->Matrix: The matrix to be calculated.                        */
/* - options->method: The method to be used (Jacobi or Gauss-Seidel).       */
/* - options->inf_func: The function to be used for initialization.         */
/* - options->term_iteration: The number of iterations for the calculation. */
/* - options->termination: The termination condition (precision/iteration). */
/* - options->term_precision: The precision for termination.                */
/* - options->number: The number of threads for parallel execution.         */
/* - results->stat_iteration: The number of iterations performed.           */
/* - results->stat_precision: The precision achieved.                       */
/* - results->m: The index of the matrix used for the final result.         */
/*                                                                          */
/* The function performs the following steps:                               */
/* 1. Initialize local variables and parameters based on the options.       */
/* 2. Set up parallel execution if the Jacobi method's selected.            */
/* 3. Perform the calculation for each element of the matrix.               */
/* 4. Update the maximum residuum value for convergence checking.           */
/* 5. Exchange the old and new matrices for the next iteration.             */
/* 6. Check for stopping conditions based on the termination method.        */
/* 7. Update the results structure with iteration count and precision.      */
/* ************************************************************************ */
static
void
calculate_element (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	int i, j;           /* local variables for loops */
	int m1, m2;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */
	double maxLocalResiduum; /* maximum residuum value of a slave in iteration in a thread */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	double** Matrix_Out;
	double** Matrix_In;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

    /* Calculate the pih and fpisin values if the function is FPISIN */
	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}
#pragma omp parallel if(options->method == METH_JACOBI) default(none) num_threads(options->number) \
    private(i, j, star, residuum, maxLocalResiduum) 	/*fpisin_i wird nicht mehr verwendet, keine unterteilung reihen und spalten, fpisin in innerer for s. */ \
    shared(arguments, m1, m2, N, options, pih, fpisin, term_iteration, results, Matrix_Out, Matrix_In, maxResiduum)
{
	while (term_iteration > 0)
	{
        /* Master thread initializes the matrices for the current iteration and resets maxResiduum */
        #pragma omp master
        {
		 	Matrix_Out  = arguments->Matrix[m1];
			Matrix_In   = arguments->Matrix[m2];
		    maxResiduum = 0;
		}
        /* Each thread initializes its local maxLocalResiduum */
		maxLocalResiduum = 0;
        #pragma omp barrier

		/* over all rows */				//hier alle elemente parallelisert
        #pragma omp for collapse(2) //collapse(2) führt die nächsten 2 schleifen gleichzeitig parallel aus
		for (i = 1; i < N; i++)
		{
			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin * sin(pih * (double)i) * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
                    if (residuum > maxLocalResiduum) {maxLocalResiduum = residuum;}
				}

				Matrix_Out[i][j] = star;
			}
		}
        /* thread calculates the maximum residuum value */
		#pragma omp atomic compare
		if (maxLocalResiduum > maxResiduum) {maxResiduum = maxLocalResiduum;}
        #pragma omp barrier
		#pragma omp master
		{

			results->stat_iteration++;
			results->stat_precision = maxResiduum;

			/* exchange m1 and m2 */
			i = m1;
			m1 = m2;
			m2 = i;

			/* check for stopping calculation depending on termination method */
			if (options->termination == TERM_PREC)
			{
				if (maxResiduum < options->term_precision)
				{
					term_iteration = 0;
				}
			}
			else if (options->termination == TERM_ITER)
			{
				term_iteration--;
			}
        }
        #pragma omp barrier
	}
} /* End #pragma parallel */

	results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauß-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %" PRIu64 "\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0)
	{
		printf("f(x,y) = 0");
	}
	else if (options->inf_func == FUNC_FPISIN)
	{
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:       ");

	if (options->termination == TERM_PREC)
	{
		printf("Hinreichende Genaugkeit");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %.11e\n", results->stat_precision);
	printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion displayMatrix:                               **/
/**                                                                        **/
/** Die Funktion displayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static
void
displayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		for (x = 0; x < 9; x++)
		{
			printf ("%11.8f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
		}

		printf ("\n");
	}

	fflush (stdout);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	askParams(&options, argc, argv);

	initVariables(&arguments, &results, &options);

	allocateMatrices(&arguments);
	initMatrices(&arguments, &options);

	gettimeofday(&start_time, NULL);
    /* Perform the calculation based on the defined method:
    * - SPALTE: Calculate by column
    * - ZEILE: Calculate by row
    * - ELEMENT: Calculate by element
    * - Default: Sequential calculation
    */
#ifdef SPALTE
	calculate_column(&arguments, &results, &options);
#elif defined(ZEILE)
	calculate_row(&arguments, &results, &options);
#elif defined(ELEMENT)
	calculate_element(&arguments, &results, &options);
#else
	calculate_seq(&arguments, &results, &options);
#endif
	gettimeofday(&comp_time, NULL);

	displayStatistics(&arguments, &results, &options);
	displayMatrix(&arguments, &results, &options);

	freeMatrices(&arguments);

	return 0;
}
