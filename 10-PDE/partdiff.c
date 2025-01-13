/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institute of Computer Science              **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff.c                                                  **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauss-Seidel and   **/
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
#include <mpi.h>
#include <omp.h>

#include "partdiff.h"

struct calculation_arguments
{
    uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
    uint64_t  num_matrices;   /* number of matrices                             */
    double    h;              /* length of a space between two lines            */
    double    ***Matrix;      /* index matrix used for addressing M             */
    double    *M;             /* two matrices with real values                  */
    uint64_t  comm_size;      /* number of processes the program handles        */
    uint64_t  rang;           /* rank of the process                            */
    uint64_t  rows;           /* number of rows each process handles            */
    int       offset;
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

static
void
calc_chunk (struct calculation_arguments* arguments)
{
    uint64_t const N = arguments->N;

    int world_size = arguments->comm_size;

    /* Calculation of the chunks each process handles. Corresponds to the number of matrix rows held by each process, NOT the number of rows to be calculated.*/
    uint64_t rows = (N + 1) / world_size;

    /* Since the matrix cannot be evenly divided among all processes except in rare cases, some processes need to handle one more row than others. In this case, the first x processes take one additional row.*/
    if((N + 1) % world_size != 0)
    {
        if(arguments->rang < ((N + 1) % world_size))
        {
            rows++;
        }
    }

    /* Since each process needs the bottom row of the previous process and the top row of the next process for calculations,
    2 additional rows are allocated to the chunk. The first and last processes only need one additional row.*/
    if(arguments->rang == 0 || (arguments->rang + 1) == arguments->comm_size)
    {
        rows++;
    }
    else
    {
        rows += 2;
    }

    /* offset describes the position of the respective row in the total matrix. Needed for the sine calculation in calculate_parallel as well as in init_matrices_parallel and display_matrices_parallel.
    i + offset should point to the global rank of the local row with index i */
    int offset = 0;
    MPI_Status status;

    if(arguments->rang == 0)
    {
        arguments->offset = offset;

        offset += rows - 2; // Subtract 2 rows, one for this process's boundary row, one for the next process's boundary row

        MPI_Send((void*) &offset, 1, MPI_INT, arguments->rang + 1, 0, MPI_COMM_WORLD);
    }
    else if(arguments->rang + 1 == arguments->comm_size)
    {
        MPI_Recv((void*) &offset, 1, MPI_INT, arguments->rang - 1, 0, MPI_COMM_WORLD, &status);

        arguments->offset = offset;
    }
    else
    {
        MPI_Recv((void*) &offset, 1, MPI_INT, arguments->rang - 1, 0, MPI_COMM_WORLD, &status);

        arguments->offset = offset;

        offset += rows - 2; // Subtract 2 rows, one for this process's boundary row, one for the next process's boundary row

        MPI_Send((void*) &offset, 1, MPI_INT, arguments->rang + 1, 0, MPI_COMM_WORLD);
    }

    //printf("rank: %" PRIu64 "\nrows: %" PRIu64 "\noffset: %d\n\n", arguments->rang, rows, arguments->offset);

    arguments->rows = rows;
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
        printf("Memory issues! (%" PRIu64 " Bytes requested)\n", size);
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
/* allocateMatrices_parallel: allocates memory for a chunk of the matrices  */
/* ************************************************************************ */
static
void
allocateMatrices_parallel (struct calculation_arguments* arguments)
{
    uint64_t i, j;

    uint64_t const N = arguments->N;
    uint64_t const rows = arguments->rows;

    arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * rows * sizeof(double));
    arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

    for (i = 0; i < arguments->num_matrices; i++)
    {
        /* Memory allocation for the row pointers */
        arguments->Matrix[i] = allocateMemory(rows * sizeof(double*));

        for (j = 0; j < rows; j++)
        {
            /* Memory allocation for the "rows" rows */
            arguments->Matrix[i][j] = arguments->M + (i * rows * (N + 1)) + (j * (N + 1));
        }
    }
}

/* ********************************************************************************* */
/* initMatrices: Initialize chunk of matrix/matrices and some global variables       */
/* ********************************************************************************* */
static
void
initMatrices_parallel (struct calculation_arguments* arguments, struct options const* options)
{
    uint64_t g, i, j; /* local variables for loops */

    uint64_t const N = arguments->N;
    uint64_t rows = arguments->rows;
    double const h = arguments->h;
    double*** Matrix = arguments->Matrix;
    int offset = arguments->offset;

    /* initialize matrix/matrices with zeros */
    for (g = 0; g < arguments->num_matrices; g++)
    {
        for (i = 0; i < rows; i++)
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
        /* Only the first and last process need to initialize their first (or last) row, but all processes must initialize the left and right boundaries */
        if(arguments->rang == 0)
        {
            for(g = 0; g < arguments->num_matrices; g++)
            {
                for(i = 0; i < rows; i++)
                {
                    /* left and right boundaries */
                    Matrix[g][i][0] = 1.0 - (h * i); // offset not needed as rank = 0
                    Matrix[g][i][N] = h * i;
                }
                Matrix[g][0][N] = 0.0;

                for(i = 0; i <= N; i++)
                {
                    /* top row */
                    Matrix[g][0][i] = 1.0 - (h * i);
                }

            }
        }
        else if(arguments->rang == arguments->comm_size - 1)
        {
            for(g = 0; g < arguments->num_matrices; g++)
            {
                for(i = 0; i < rows; i++)
                {
                    /* left and right boundaries */
                    Matrix[g][i][0] = 1.0 - (h * (i + offset));
                    Matrix[g][i][N] = h * (i + offset);
                }

                for(i = 0; i <= N; i++)
                {
                    /* bottom row */
                    Matrix[g][rows - 1][i] = h * i;
                }
            }
        }
        else
        {
            for(g = 0; g < arguments->num_matrices; g++)
            {
                for(i = 0; i < rows; i++)
                {
                    /* left and right boundaries */
                    Matrix[g][i][0] = 1.0 - (h * (i + offset));
                    Matrix[g][i][N] = h * (i + offset);
                }
            }
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
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
    int i, j;           /* local variables for loops */
    int m1, m2;         /* used as indices for old and new matrices */
    double star;        /* four times center value minus 4 neighbor values */
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
/* calculate_parallel_gauss: Solves the equation in parallel                */
/* ************************************************************************ */
static
void
calculate_parallel_gauss(
    struct calculation_arguments const* arguments,
    struct calculation_results* results,
    struct options const* options)
{
    int i, j;           /* local variables for loops */
    double star;        /* four times center value minus 4 neighbor values */
    double residuum;    /* residuum of current iteration */
    double maxResiduum; /* maximum residuum value of a slave in iteration */
    
    /* initialize m1 and m2 for Gauss */
    // int m1 = 0;
    // int m2 = 0;

    int const N = arguments->N;
    double const h = arguments->h;

    int const offset = arguments->offset;

    /* MPI stuff */
    uint64_t const rank = arguments->rang;
    uint64_t const comm_size = arguments->comm_size;
    MPI_Status status;

    /* Decide send directions */
    int const send_up           = (rank < (comm_size - 1));
    int const send_down         = (rank > 0);
    int const receive_from_up   = send_up;
    int const receive_from_down = send_down;

    double pih = 0.0;
    double fpisin = 0.0;

    int term_iteration = options->term_iteration;

    if (options->inf_func == FUNC_FPISIN)
    {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    while (term_iteration > 0)
    {
        double** Matrix = arguments->Matrix[0];

        MPI_Request receive_up;
        MPI_Request receive_down;

        if (receive_from_up && (results->stat_iteration > 0))
        {
            MPI_Irecv(
                (void*) &Matrix[arguments->rows - 1][0],
                N + 1,
                MPI_DOUBLE,
                rank + 1,
                0,
                MPI_COMM_WORLD,
                &receive_up
            );
        }

        if (receive_from_down)
        {
            MPI_Irecv(
                (void*) &Matrix[0][0],
                N + 1,
                MPI_DOUBLE,
                rank - 1,
                1,
                MPI_COMM_WORLD,
                &receive_down
            );
        }

        if (receive_from_down)
        {
            MPI_Wait(&receive_down, &status);
        }

        if (receive_from_up && (results->stat_iteration > 0))
        {
            MPI_Wait(&receive_up, &status);
        }

        // Border line request
        MPI_Request request_up;
        MPI_Request request_down;
        // MPI_Request request_residuum;

        maxResiduum = 0;

        /* over all rows */
        for (i = 1; (uint64_t)i < arguments->rows - 1; i++)
        {
            double fpisin_i = 0.0;

            if (options->inf_func == FUNC_FPISIN)
            {
                fpisin_i = fpisin * sin(pih * (double)(i + offset));
            }

            /* over all columns */
            for (j = 1; j < N; j++)
            {
                star = (Matrix[i-1][j] + Matrix[i][j-1] +
                        Matrix[i][j+1] + Matrix[i+1][j]) * 0.25;

                if (options->inf_func == FUNC_FPISIN)
                {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options->termination == TERM_PREC || term_iteration == 1)
                {
                    residuum = Matrix[i][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
                }

                Matrix[i][j] = star;
            }

            if ((i == 1) && send_up)
            {
                MPI_Isend(
                    (void*) &Matrix[arguments->rows - 2][0],
                    N + 1,
                    MPI_DOUBLE,
                    rank + 1,
                    1,
                    MPI_COMM_WORLD,
                    &request_up
                );
            }
        }

        results->stat_iteration++;
        results->stat_precision = maxResiduum;

        if (send_down)
        {
            MPI_Isend(
                (void*) &Matrix[1][0],
                N + 1,
                MPI_DOUBLE,
                rank - 1,
                0,
                MPI_COMM_WORLD,
                &request_down
            );
        }

        if (send_up)
        {
            MPI_Wait(&request_up, &status);
        }

        if (send_down)
        {
            MPI_Wait(&request_down, &status);
        }

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

    results->m = 0;
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

    printf("Computation time:    %f s \n", time);
    printf("Memory usage:        %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
    printf("Computation method:  ");

    if (options->method == METH_GAUSS_SEIDEL)
    {
        printf("Gauss-Seidel");
    }
    else if (options->method == METH_JACOBI)
    {
        printf("Jacobi");
    }

    printf("\n");
    printf("Interlines:          %" PRIu64 "\n",options->interlines);
    printf("Source function:     ");

    if (options->inf_func == FUNC_F0)
    {
        printf("f(x,y) = 0");
    }
    else if (options->inf_func == FUNC_FPISIN)
    {
        printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
    }

    printf("\n");
    printf("Termination:         ");

    if (options->termination == TERM_PREC)
    {
        printf("Sufficient Precision");
    }
    else if (options->termination == TERM_ITER)
    {
        printf("Number of Iterations");
    }

    printf("\n");
    printf("Number of Iterations: %" PRIu64 "\n", results->stat_iteration);
    printf("Error Norm:           %e\n", results->stat_precision);
    printf("\n");
}

/* ********************************************************************************* */
/*  displayStatistics_parallel: displays some statistics about the calculation       */
/* ********************************************************************************* */
static
void
displayStatistics_parallel (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
    int N = arguments->N;
    double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

    double memory = (N + 1) * arguments->rows * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0;
    double total_memory;

    /* All processes synchronously send their memory usage to process 0. It sums the values. */
    MPI_Reduce((void*) &memory, &total_memory, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double error = 0.0;

    /* All processes synchronously send their error to process 0. It stores the maximum. */
    MPI_Reduce((void*) &results->stat_precision, &error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(arguments->rang == 0)
    {
        printf("Computation time:    %f s \n", time);
        printf("Memory usage:        %f MiB\n", total_memory);
        printf("Computation method:  ");

        if (options->method == METH_GAUSS_SEIDEL)
        {
            printf("Gauss-Seidel");
        }
        else if (options->method == METH_JACOBI)
        {
            printf("Jacobi");
        }

        printf("\n");
        printf("Interlines:          %" PRIu64 "\n",options->interlines);
        printf("Source function:     ");

        if (options->inf_func == FUNC_F0)
        {
            printf("f(x,y) = 0");
        }
        else if (options->inf_func == FUNC_FPISIN)
        {
            printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
        }

        printf("\n");
        printf("Termination:         ");

        if (options->termination == TERM_PREC)
        {
            printf("Sufficient Precision");
        }
        else if (options->termination == TERM_ITER)
        {
            printf("Number of Iterations");
        }

        printf("\n");
        printf("Number of Iterations: %" PRIu64 "\n", results->stat_iteration);
        printf("Error Norm:           %e\n", error);
        printf("\n");
    }
}

/****************************************************************************/
/** Description of the function displayMatrix:                             **/
/**                                                                        **/
/** The function displayMatrix prints a matrix                             **/
/** in a "clear" manner to the standard output.                            **/
/**                                                                        **/
/** Clarity is achieved by only printing a part of the matrix,              **/
/** specifically the boundary rows/columns and seven intermediate rows.   **/
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
            printf ("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
        }

        printf ("\n");
    }

    fflush (stdout);
}

/****************************************************************************/
/** Description of the function displayMatrix_parallel                       **/
/**                                                                        **/
/** The function displayMatrix_parallel prints a matrix                    **/
/** in a "clear" manner to the standard output.                            **/
/**                                                                        **/
/** Clarity is achieved by only printing a part of the matrix,              **/
/** specifically the boundary rows/columns and seven intermediate rows.   **/
/****************************************************************************/
static
void
displayMatrix_parallel (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
    int x, y;

    double** Matrix = arguments->Matrix[results->m];

    int const interlines = options->interlines;

    MPI_Status status;

    if(arguments->rang == 0)
    {
        for (y = 0; y < 9 && (uint64_t)(y * (interlines + 1)) < arguments->rows; y++)
        {
            for (x = 0; x < 9; x++)
            {
                printf ("%7.4f", Matrix[(y * (interlines + 1))][x * (interlines + 1)]);
            }

            printf ("\n");
        }

        fflush (stdout);

        MPI_Send((void*) &y, 1, MPI_INT, arguments->rang + 1, 6, MPI_COMM_WORLD);
    }
    else if(arguments->rang == arguments->comm_size - 1)
    {
        MPI_Recv((void*) &y, 1, MPI_INT, arguments->rang - 1, 6, MPI_COMM_WORLD, &status);

        for (; y < 9 && (uint64_t)(y * (interlines + 1)) < arguments->offset + arguments->rows && (y * (interlines + 1)) > arguments->offset; y++)
        {
            for (x = 0; x < 9; x++)
            {
                printf ("%7.4f", Matrix[(y * (interlines + 1)) - arguments->offset][x * (interlines + 1)]);
            }

            printf ("\n");
        }

        fflush (stdout);
    }
    else
    {
        MPI_Recv((void*) &y, 1, MPI_INT, arguments->rang - 1, 6, MPI_COMM_WORLD, &status);

        for (; y < 9 && (uint64_t)(y * (interlines + 1)) < arguments->offset + arguments->rows && (y * (interlines + 1)) > arguments->offset; y++)
        {
            for (x = 0; x < 9; x++)
            {
                printf ("%7.4f", Matrix[(y * (interlines + 1)) - arguments->offset][x * (interlines + 1)]);
            }

            printf ("\n");
        }

        fflush (stdout);

        MPI_Send((void*) &y, 1, MPI_INT, arguments->rang + 1, 6, MPI_COMM_WORLD);
    }
}

static
void
sequential_block(struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
    allocateMatrices(arguments);
    initMatrices(arguments, options);

    calculate(arguments, results, options);

    gettimeofday(&comp_time, NULL);

    displayStatistics(arguments, results, options);
    displayMatrix(arguments, results, options);

    freeMatrices(arguments);
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

    MPI_Init(NULL, NULL);

    /* Query and store the process rank */
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    arguments.rang = world_rank;

    /* Query and store the total number of processes */
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    arguments.comm_size = world_size;

    if (world_size > 1)
    {
        /* Process 0 takes the parameters */
        if (world_rank == 0)
        {
            askParams(&options, argc, argv);
        }

        /* Process 0 broadcasts the parameters to all other processes */
        MPI_Bcast((void*) &options.number,         1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Bcast((void*) &options.method,         1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Bcast((void*) &options.interlines,     1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Bcast((void*) &options.inf_func,       1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Bcast((void*) &options.termination,    1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Bcast((void*) &options.term_iteration, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Bcast((void*) &options.term_precision, 1, MPI_DOUBLE,   0, MPI_COMM_WORLD);

        initVariables(&arguments, &results, &options);

        /* Calculate the number of rows each process handles */
        calc_chunk(&arguments);

        allocateMatrices_parallel(&arguments);

        initMatrices_parallel(&arguments, &options);

        gettimeofday(&start_time, NULL);

        if (options.method == METH_JACOBI)
        {
            calculate(&arguments, &results, &options);
        }
        else
        {
            calculate_parallel_gauss(&arguments, &results, &options);
        }

        /* Barrier to keep measurement consistent */
        MPI_Barrier(MPI_COMM_WORLD);

        gettimeofday(&comp_time, NULL);

        displayStatistics_parallel(&arguments, &results, &options);
        //printf("%d going to display_matrix\n", world_rank);
        displayMatrix_parallel(&arguments, &results, &options);

        freeMatrices(&arguments);
    }
    else
    {
        askParams(&options, argc, argv);
        initVariables(&arguments, &results, &options);
        sequential_block(&arguments, &results, &options);
    }

    MPI_Finalize();

    return 0;
}