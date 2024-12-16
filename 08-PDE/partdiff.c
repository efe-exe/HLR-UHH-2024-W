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
#include <mpi.h> // NEW

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

// NEW
struct mpi_parameters
{
	int world_size; // Anzahl der Prozesse
	int world_rank; // Rang des Prozesses
	uint64_t start_r; // Startzeile des Prozesses
	uint64_t end_r; // Endzeile des Prozesses
	uint64_t N_r; // Anzahl der Zeilen des Prozesses
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
allocateMatrices (struct calculation_arguments* arguments, struct mpi_parameters* parameters) /* NEW */
{
    uint64_t i, j;

    uint64_t const N = arguments->N;
    uint64_t const N_r = parameters->N_r; /* NEW */

	arguments->M = allocateMemory(arguments->num_matrices * (N_r + 2) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++) {
		arguments->Matrix[i] = allocateMemory((N_r + 2) * sizeof(double*)); /* NEW */
		for (j = 0; j <= N_r + 1; j++) {
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N_r + 2)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options, struct mpi_parameters* parameters) /* NEW */
{
	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	uint64_t const N_r = parameters->N_r; /* NEW */
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N_r + 1; i++)
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
		uint64_t start_r = parameters->start_r; /* NEW */
		for(i = 0; i < N_r + 2; i++) /* NEW */
		{
			for(j = 0; j < arguments->num_matrices; j++)
			{
				Matrix[j][i][0] = 1.0 + (1.0 - (h * (i + start_r - 1))); // linker Rand /* NEW */
				Matrix[j][i][N] = 1.0 - h * (i + start_r - 1); // rechter Rand /* NEW */
			}
		}

		if (parameters->world_rank == 0) /* NEW */
		{
			for(i = 0; i < N; i++)
			{
				for (j = 0; j < arguments->num_matrices; j++)
				{
					Matrix[j][0][N - i] = 1 + h * i; // Obere Kante
				}
			}
		}

		if (parameters->world_rank == parameters->world_size - 1) /* NEW */
		{
			for(i = 0; i < N; i++)
			{
				for (j = 0; j < arguments->num_matrices; j++)
				{
					Matrix[j][N_r + 1][i] = 1 - (h * i); // Untere Kante
				}
			}
		}
    }
}


/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
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

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

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
/* calculateJacobiMPI: solves the equation for Jacobi with MPI              */
/* ************************************************************************ */
static
void
calculateJacobiMPI(struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, struct mpi_parameters* parameters)
{
    int i, j;           /* local variables for loops */
    int m1, m2;         /* used as indices for old and new matrices */
    double star;        /* four times center value minus 4 neigh.b values */
    double residuum;    /* residuum of current iteration */
    double maxResiduum; /* maximum residuum value of a slave in iteration */

    int const N_r = parameters->N_r; /* NEW */
    int const N = arguments->N;
    double const h = arguments->h;

    double pih = 0.0;
    double fpisin = 0.0;

    int rank = parameters->world_rank; /* NEW */
    int size = parameters->world_size; /* NEW */

    int term_iteration = options->term_iteration;

    /* initialize m1 and m2 depending on Jacobi NEW */
    m1 = 0;
    m2 = 1;

    if (options->inf_func == FUNC_FPISIN)
    {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    while (term_iteration > 0)
    {
        double** Matrix_Out = arguments->Matrix[m1];
        double** Matrix_In  = arguments->Matrix[m2];

        maxResiduum = 0;

        /* over all rows */
        for (i = 1; i < N_r + 1; i++)
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

    	int prev, next;           /* rank ids of the previous/next rank or MPI_PROC_NULL if it doesn't exists NEW */
        prev = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
        next = (rank < size - 1) ? rank + 1 : MPI_PROC_NULL;

		/* Use MPI_Sendrecv to exchange boundary rows NEW */
        MPI_Sendrecv(&Matrix_Out[1][0], N, MPI_DOUBLE, prev, 0,
                     &Matrix_Out[N_r + 1][0], N, MPI_DOUBLE, next, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&Matrix_Out[N_r][0], N, MPI_DOUBLE, next, 1,
                     &Matrix_Out[0][0], N, MPI_DOUBLE, prev, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (options->termination == TERM_PREC || term_iteration == 1)
        {
            double all_maxResiduum = 0;
            MPI_Allreduce(&maxResiduum, &all_maxResiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            maxResiduum = all_maxResiduum;
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
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments* arguments, struct calculation_results const* results, struct options const* options, double total_memory)
{
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", total_memory);
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
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
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
displayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, struct mpi_parameters* parameters)
{
	int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	uint64_t start_r = parameters->start_r;
	uint64_t end_r = parameters->end_r;

	if (parameters->world_rank == 0)
	{
		printf("Matrix:\n");
		// erste Zeile
		for (x = 0; x < 9; x++)
			{
				printf ("%7.4f", Matrix[0][x * (interlines + 1)]);
			}
			printf ("\n");
	}

	for (y = 1; y < 8; y++)
	{
		// Rückwärtszählen, da die erste Zeile bereits ausgegeben wurde
		uint64_t line = y * (interlines + 1);
		if (line < end_r && line >= start_r)
		{
			for (x = 0; x < 9; x++)
			{
				printf ("%7.4f", Matrix[line - start_r + 1][x * (interlines + 1)]);
			}
			printf ("\n");
		}
	}

	if (parameters->world_rank == parameters->world_size - 1)
	{
		// letzte Zeile
		for (x = 0; x < 9; x++)
			{
				printf ("%7.4f", Matrix[parameters->N_r + 1][x * (interlines + 1)]);
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
	MPI_Init(NULL, NULL);

	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;
	struct mpi_parameters parameters;

    MPI_Comm_size(MPI_COMM_WORLD, &parameters.world_size); // size of process NEW
    MPI_Comm_rank(MPI_COMM_WORLD, &parameters.world_rank); // rank of process NEW

	// Signal for displaying matrix NEW
	int signal = 0;

	askParams(&options, argc, argv);
	initVariables(&arguments, &results, &options);

	parameters.start_r = (parameters.world_rank*(arguments.N - 1))/parameters.world_size + 1;
	parameters.end_r = ((parameters.world_rank + 1)*(arguments.N - 1))/parameters.world_size + 1;
	parameters.N_r = parameters.end_r - parameters.start_r;

	if (options.method == METH_GAUSS_SEIDEL || parameters.world_size == 1) // sequenziell
	{
		if (parameters.world_rank == 0)
		{
			parameters.start_r = 1;
			parameters.end_r = arguments.N;
			parameters.N_r = parameters.end_r - parameters.start_r;

			allocateMatrices(&arguments, &parameters);
			initMatrices(&arguments, &options, &parameters);
			gettimeofday(&start_time, NULL);
			calculate(&arguments, &results, &options);
			gettimeofday(&comp_time, NULL);

			
			double total_memory = (arguments.N + 1) * (parameters.N_r + 2) * sizeof(double) * arguments.num_matrices / 1024.0 / 1024.0;
			displayStatistics(&arguments, &results, &options, total_memory);
			displayMatrix(&arguments, &results, &options, &parameters);
			MPI_Barrier(MPI_COMM_WORLD);
			freeMatrices(&arguments);
		}
		else
		{
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	else
	{
		// parallel
		allocateMatrices(&arguments, &parameters);
		initMatrices(&arguments, &options, &parameters);

		gettimeofday(&start_time, NULL);
		// barrier for time calc
		MPI_Barrier(MPI_COMM_WORLD);

		calculateJacobiMPI(&arguments, &results, &options, &parameters);

		MPI_Barrier(MPI_COMM_WORLD);
		gettimeofday(&comp_time, NULL);

		// total_memory ist Summe der einzelnen memories
		double memory_r = (arguments.N + 1) * (parameters.N_r + 2) * sizeof(double) * arguments.num_matrices / 1024.0 / 1024.0;
		double total_memory;
		MPI_Allreduce(&memory_r, &total_memory, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		if (parameters.world_rank == 0)
		{
			displayStatistics(&arguments, &results, &options, total_memory);
			displayMatrix(&arguments, &results, &options, &parameters);
			MPI_Ssend(&signal, 1, MPI_INT, parameters.world_rank + 1, 0, MPI_COMM_WORLD);
		}
		else if (parameters.world_rank < parameters.world_size && parameters.world_rank > 0)
		{
			MPI_Recv(&signal, 1, MPI_INT, parameters.world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			displayMatrix(&arguments, &results, &options, &parameters);
			if (parameters.world_rank < parameters.world_size - 1)
			{
				MPI_Ssend(&signal, 1, MPI_INT, parameters.world_rank + 1, 0, MPI_COMM_WORLD);
			}
		}
			MPI_Barrier(MPI_COMM_WORLD);
			freeMatrices(&arguments);
	}
	MPI_Finalize();
	return 0;
}