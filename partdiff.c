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
    arguments->N = (options->interlines * 8) + 9 - 1; // Berechnet die Anzahl der Zeilen basierend auf den Interlines-Optionen
    arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1; // Setzt die Anzahl der Matrizen basierend auf der Methode
    arguments->h = 1.0 / arguments->N; // Berechnet die Länge eines Raums zwischen zwei Zeilen

    results->m = 0; // Initialisiert die Matrix-Nummer
    results->stat_iteration = 0; // Initialisiert die Iterationsanzahl
    results->stat_precision = 0; // Initialisiert die Präzision
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
    uint64_t i; // Deklariert eine Schleifenvariable

    for (i = 0; i < arguments->num_matrices; i++) // Schleife über die Anzahl der Matrizen
    {
        free(arguments->Matrix[i]); // Gibt den Speicher für jede Matrix frei
    }

    free(arguments->Matrix); // Gibt den Speicher für das Matrix-Array frei
    free(arguments->M); // Gibt den Speicher für die Matrix-Daten frei
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
    void *p; // Deklariert einen Zeiger für den Speicher

    if ((p = malloc(size)) == NULL) // Versucht, Speicher zuzuweisen und prüft auf Fehler
    {
        printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size); // Gibt eine Fehlermeldung aus
        exit(1); // Beendet das Programm bei Speicherfehler
    }

    return p; // Gibt den zugewiesenen Speicher zurück
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
    uint64_t i, j; // Deklariert Schleifenvariablen

    uint64_t const N = arguments->N; // Holt die Anzahl der Zeilen

    arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double)); // Weist Speicher für die Matrix-Daten zu
    arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**)); // Weist Speicher für das Matrix-Array zu

    for (i = 0; i < arguments->num_matrices; i++) // Schleife über die Anzahl der Matrizen
    {
        arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*)); // Weist Speicher für jede Matrix zu

        for (j = 0; j <= N; j++) // Schleife über die Zeilen jeder Matrix
        {
            arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1)); // Setzt die Zeiger für jede Zeile
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
    uint64_t g, i, j; // Deklariert Schleifenvariablen

    uint64_t const N = arguments->N; // Holt die Anzahl der Zeilen
    double const h = arguments->h; // Holt die Länge eines Raums zwischen zwei Zeilen
    double*** Matrix = arguments->Matrix; // Holt das Matrix-Array

    /* initialize matrix/matrices with zeros */
    for (g = 0; g < arguments->num_matrices; g++) // Schleife über die Anzahl der Matrizen
    {
        for (i = 0; i <= N; i++) // Schleife über die Zeilen jeder Matrix
        {
            for (j = 0; j <= N; j++) // Schleife über die Spalten jeder Matrix
            {
                Matrix[g][i][j] = 0.0; // Setzt jedes Element auf 0
            }
        }
    }

    /* initialize borders, depending on function (function 2: nothing to do) */
    if (options->inf_func == FUNC_F0) // Prüft, ob die Funktion FUNC_F0 ist
    {
        for (g = 0; g < arguments->num_matrices; g++) // Schleife über die Anzahl der Matrizen
        {
            for (i = 0; i <= N; i++) // Schleife über die Zeilen jeder Matrix
            {
                Matrix[g][i][0] = 1.0 - (h * i); // Setzt den linken Rand
                Matrix[g][i][N] = h * i; // Setzt den rechten Rand
                Matrix[g][0][i] = 1.0 - (h * i); // Setzt den oberen Rand
                Matrix[g][N][i] = h * i; // Setzt den unteren Rand
            }

            Matrix[g][N][0] = 0.0; // Setzt die untere linke Ecke
            Matrix[g][0][N] = 0.0; // Setzt die obere rechte Ecke
        }
    }
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate_seq (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
    int i, j;           // Lokale Variablen für Schleifen
    int m1, m2;         // Indizes für alte und neue Matrizen
    double star;        // Viermal der Mittelwert minus 4 Nachbarwerte
    double residuum;    // Residuum der aktuellen Iteration
    double maxResiduum; // Maximales Residuum eines Slaves in der Iteration

    int const N = arguments->N; // Anzahl der Zeilen
    double const h = arguments->h; // Länge eines Raums zwischen zwei Zeilen

    double pih = 0.0; // Variable für pi*h
    double fpisin = 0.0; // Variable für fpisin

    int term_iteration = options->term_iteration; // Anzahl der Iterationen

    /* initialize m1 and m2 depending on algorithm */
    if (options->method == METH_JACOBI) // Prüft, ob die Methode Jacobi ist
    {
        m1 = 0; // Setzt m1 auf 0
        m2 = 1; // Setzt m2 auf 1
    }
    else // Wenn die Methode nicht Jacobi ist
    {
        m1 = 0; // Setzt m1 auf 0
        m2 = 0; // Setzt m2 auf 0
    }

    if (options->inf_func == FUNC_FPISIN) // Prüft, ob die Funktion FUNC_FPISIN ist
    {
        pih = PI * h; // Berechnet pi*h
        fpisin = 0.25 * TWO_PI_SQUARE * h * h; // Berechnet fpisin
    }

    while (term_iteration > 0) // Schleife über die Iterationen
    {
        double** Matrix_Out = arguments->Matrix[m1]; // Holt die Ausgangsmatrix
        double** Matrix_In  = arguments->Matrix[m2]; // Holt die Eingabematrix

        maxResiduum = 0; // Setzt das maximale Residuum auf 0

        /* over all rows */
        for (i = 1; i < N; i++) // Schleife über die Zeilen
        {
            double fpisin_i = 0.0; // Variable für fpisin*i

            if (options->inf_func == FUNC_FPISIN) // Prüft, ob die Funktion FUNC_FPISIN ist
            {
                fpisin_i = fpisin * sin(pih * (double)i); // Berechnet fpisin*i
            }

            /* over all columns */
            for (j = 1; j < N; j++) // Schleife über die Spalten
            {
                star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]); // Berechnet den neuen Wert

                if (options->inf_func == FUNC_FPISIN) // Prüft, ob die Funktion FUNC_FPISIN ist
                {
                    star += fpisin_i * sin(pih * (double)j); // Fügt den fpisin-Wert hinzu
                }

                if (options->termination == TERM_PREC || term_iteration == 1) // Prüft, ob die Termination PREC ist oder die letzte Iteration
                {
                    residuum = Matrix_In[i][j] - star; // Berechnet das Residuum
                    residuum = (residuum < 0) ? -residuum : residuum; // Nimmt den absoluten Wert des Residuum
                    maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum; // Aktualisiert das maximale Residuum
                }

                Matrix_Out[i][j] = star; // Setzt den neuen Wert in die Ausgangsmatrix
            }
        }

        results->stat_iteration++; // Erhöht die Iterationsanzahl
        results->stat_precision = maxResiduum; // Setzt die Präzision

        /* exchange m1 and m2 */
        i = m1; // Tauscht m1 und m2
        m1 = m2; // Tauscht m1 und m2
        m2 = i; // Tauscht m1 und m2

        /* check for stopping calculation depending on termination method */
        if (options->termination == TERM_PREC) // Prüft, ob die Termination PREC ist
        {
            if (maxResiduum < options->term_precision) // Prüft, ob das maximale Residuum kleiner als die Präzision ist
            {
                term_iteration = 0; // Setzt die Iterationen auf 0
            }
        }
        else if (options->termination == TERM_ITER) // Prüft, ob die Termination ITER ist
        {
            term_iteration--; // Verringert die Iterationen
        }
    }

    results->m = m2; // Setzt die Matrix-Nummer
}
static
void
calculate_row (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
    int i, j;           // Lokale Variablen für Schleifen
    int m1, m2;         // Indizes für alte und neue Matrizen
    double star;        // Viermal der Mittelwert minus 4 Nachbarwerte
    double residuum;    // Residuum der aktuellen Iteration
    double maxResiduum; // Maximales Residuum eines Slaves in der Iteration
    double maxLocalResiduum; // Maximales Residuum eines Slaves in der Iteration in einem Thread

    int const N = arguments->N; // Anzahl der Zeilen
    double const h = arguments->h; // Länge eines Raums zwischen zwei Zeilen

    double pih = 0.0; // Variable für pi*h
    double fpisin = 0.0; // Variable für fpisin

    double** Matrix_Out; // Zeiger auf die Ausgangsmatrix
    double** Matrix_In; // Zeiger auf die Eingabematrix

    int term_iteration = options->term_iteration; // Anzahl der Iterationen

    /* initialize m1 and m2 depending on algorithm */
    if (options->method == METH_JACOBI) // Prüft, ob die Methode Jacobi ist
    {
        m1 = 0; // Setzt m1 auf 0
        m2 = 1; // Setzt m2 auf 1
    }
    else // Wenn die Methode nicht Jacobi ist
    {
        m1 = 0; // Setzt m1 auf 0
        m2 = 0; // Setzt m2 auf 0
    }

    if (options->inf_func == FUNC_FPISIN) // Prüft, ob die Funktion FUNC_FPISIN ist
    {
        pih = PI * h; // Berechnet pi*h
        fpisin = 0.25 * TWO_PI_SQUARE * h * h; // Berechnet fpisin
    }
#pragma omp parallel if(options->method == METH_JACOBI) default(none) num_threads(options->number) \
    private(i, j, star, residuum, maxLocalResiduum) \
    shared(arguments, m1, m2, N, options, pih, fpisin, term_iteration, results, Matrix_Out, Matrix_In, maxResiduum)
{
    while (term_iteration > 0) // Schleife über die Iterationen
    {
        #pragma omp master
        {
         	Matrix_Out  = arguments->Matrix[m1]; // Holt die Ausgangsmatrix
            Matrix_In   = arguments->Matrix[m2]; // Holt die Eingabematrix
            maxResiduum = 0; // Setzt das maximale Residuum auf 0
        }
        maxLocalResiduum = 0; // Setzt das lokale maximale Residuum auf 0
        #pragma omp barrier

        /* over all rows */
        #pragma omp for
        for (i = 1; i < N; i++) // Schleife über die Zeilen
        {
            double fpisin_i = 0.0; // Variable für fpisin*i

            if (options->inf_func == FUNC_FPISIN) // Prüft, ob die Funktion FUNC_FPISIN ist
            {
                fpisin_i = fpisin * sin(pih * (double)i); // Berechnet fpisin*i
            }

            /* over all columns */
            for (j = 1; j < N; j++) // Schleife über die Spalten
            {
                star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]); // Berechnet den neuen Wert

                if (options->inf_func == FUNC_FPISIN) // Prüft, ob die Funktion FUNC_FPISIN ist
                {
                    star += fpisin_i * sin(pih * (double)j); // Fügt den fpisin-Wert hinzu
                }

                if (options->termination == TERM_PREC || term_iteration == 1) // Prüft, ob die Termination PREC ist oder die letzte Iteration
                {
                    residuum = Matrix_In[i][j] - star; // Berechnet das Residuum
                    residuum = (residuum < 0) ? -residuum : residuum; // Nimmt den absoluten Wert des Residuum
                    if (residuum > maxLocalResiduum) {maxLocalResiduum = residuum;} // Aktualisiert das lokale maximale Residuum
                }

                Matrix_Out[i][j] = star; // Setzt den neuen Wert in die Ausgangsmatrix
            }
        }
        #pragma omp atomic compare
        if (maxLocalResiduum > maxResiduum) {maxResiduum = maxLocalResiduum;} // Aktualisiert das maximale Residuum
        #pragma omp barrier
        #pragma omp master
        {

            results->stat_iteration++; // Erhöht die Iterationsanzahl
            results->stat_precision = maxResiduum; // Setzt die Präzision

            /* exchange m1 and m2 */
            i = m1; // Tauscht m1 und m2
            m1 = m2; // Tauscht m1 und m2
            m2 = i; // Tauscht m1 und m2

            /* check for stopping calculation depending on termination method */
            if (options->termination == TERM_PREC) // Prüft, ob die Termination PREC ist
            {
                if (maxResiduum < options->term_precision) // Prüft, ob das maximale Residuum kleiner als die Präzision ist
                {
                    term_iteration = 0; // Setzt die Iterationen auf 0
                }
            }
            else if (options->termination == TERM_ITER) // Prüft, ob die Termination ITER ist
            {
                term_iteration--; // Verringert die Iterationen
            }
        }
        #pragma omp barrier
    }
} /* End #pragma parallel */

    results->m = m2; // Setzt die Matrix-Nummer
}
static
void
calculate_column (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
    int i, j;           // Lokale Variablen für Schleifen
    int m1, m2;         // Indizes für alte und neue Matrizen
    double star;        // Vier/****************************************************************************/
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
    arguments->N = (options->interlines * 8) + 9 - 1; // Berechnet die Anzahl der Zeilen basierend auf den Interlines-Optionen
    arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1; // Setzt die Anzahl der Matrizen basierend auf der Methode
    arguments->h = 1.0 / arguments->N; // Berechnet die Länge eines Raums zwischen zwei Zeilen

    results->m = 0; // Initialisiert die Matrix-Nummer
    results->stat_iteration = 0; // Initialisiert die Iterationsanzahl
    results->stat_precision = 0; // Initialisiert die Präzision
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
    uint64_t i; // Deklariert eine Schleifenvariable

    for (i = 0; i < arguments->num_matrices; i++) // Schleife über die Anzahl der Matrizen
    {
        free(arguments->Matrix[i]); // Gibt den Speicher für jede Matrix frei
    }

    free(arguments->Matrix); // Gibt den Speicher für das Matrix-Array frei
    free(arguments->M); // Gibt den Speicher für die Matrix-Daten frei
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
    void *p; // Deklariert einen Zeiger für den Speicher

    if ((p = malloc(size)) == NULL) // Versucht, Speicher zuzuweisen und prüft auf Fehler
    {
        printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size); // Gibt eine Fehlermeldung aus
        exit(1); // Beendet das Programm bei Speicherfehler
    }

    return p; // Gibt den zugewiesenen Speicher zurück
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
    uint64_t i, j; // Deklariert Schleifenvariablen

    uint64_t const N = arguments->N; // Holt die Anzahl der Zeilen

    arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double)); // Weist Speicher für die Matrix-Daten zu
    arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**)); // Weist Speicher für das Matrix-Array zu

    for (i = 0; i < arguments->num_matrices; i++) // Schleife über die Anzahl der Matrizen
    {
        arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*)); // Weist Speicher für jede Matrix zu

        for (j = 0; j <= N; j++) // Schleife über die Zeilen jeder Matrix
        {
            arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1)); // Setzt die Zeiger für jede Zeile
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
    uint64_t g, i, j; // Deklariert Schleifenvariablen

    uint64_t const N = arguments->N; // Holt die Anzahl der Zeilen
    double const h = arguments->h; // Holt die Länge eines Raums zwischen zwei Zeilen
    double*** Matrix = arguments->Matrix; // Holt das Matrix-Array

    /* initialize matrix/matrices with zeros */
    for (g = 0; g < arguments->num_matrices; g++) // Schleife über die Anzahl der Matrizen
    {
        for (i = 0; i <= N; i++) // Schleife über die Zeilen jeder Matrix
        {
            for (j = 0; j <= N; j++) // Schleife über die Spalten jeder Matrix
            {
                Matrix[g][i][j] = 0.0; // Setzt jedes Element auf 0
            }
        }
    }

    /* initialize borders, depending on function (function 2: nothing to do) */
    if (options->inf_func == FUNC_F0) // Prüft, ob die Funktion FUNC_F0 ist
    {
        for (g = 0; g < arguments->num_matrices; g++) // Schleife über die Anzahl der Matrizen
        {
            for (i = 0; i <= N; i++) // Schleife über die Zeilen jeder Matrix
            {
                Matrix[g][i][0] = 1.0 - (h * i); // Setzt den linken Rand
                Matrix[g][i][N] = h * i; // Setzt den rechten Rand
                Matrix[g][0][i] = 1.0 - (h * i); // Setzt den oberen Rand
                Matrix[g][N][i] = h * i; // Setzt den unteren Rand
            }

            Matrix[g][N][0] = 0.0; // Setzt die untere linke Ecke
            Matrix[g][0][N] = 0.0; // Setzt die obere rechte Ecke
        }
    }
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate_seq (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
    int i, j;           // Lokale Variablen für Schleifen
    int m1, m2;         // Indizes für alte und neue Matrizen
    double star;        // Viermal der Mittelwert minus 4 Nachbarwerte
    double residuum;    // Residuum der aktuellen Iteration
    double maxResiduum; // Maximales Residuum eines Slaves in der Iteration

    int const N = arguments->N; // Anzahl der Zeilen
    double const h = arguments->h; // Länge eines Raums zwischen zwei Zeilen

    double pih = 0.0; // Variable für pi*h
    double fpisin = 0.0; // Variable für fpisin

    int term_iteration = options->term_iteration; // Anzahl der Iterationen

    /* initialize m1 and m2 depending on algorithm */
    if (options->method == METH_JACOBI) // Prüft, ob die Methode Jacobi ist
    {
        m1 = 0; // Setzt m1 auf 0
        m2 = 1; // Setzt m2 auf 1
    }
    else // Wenn die Methode nicht Jacobi ist
    {
        m1 = 0; // Setzt m1 auf 0
        m2 = 0; // Setzt m2 auf 0
    }

    if (options->inf_func == FUNC_FPISIN) // Prüft, ob die Funktion FUNC_FPISIN ist
    {
        pih = PI * h; // Berechnet pi*h
        fpisin = 0.25 * TWO_PI_SQUARE * h * h; // Berechnet fpisin
    }

    while (term_iteration > 0) // Schleife über die Iterationen
    {
        double** Matrix_Out = arguments->Matrix[m1]; // Holt die Ausgangsmatrix
        double** Matrix_In  = arguments->Matrix[m2]; // Holt die Eingabematrix

        maxResiduum = 0; // Setzt das maximale Residuum auf 0

        /* over all rows */
        for (i = 1; i < N; i++) // Schleife über die Zeilen
        {
            double fpisin_i = 0.0; // Variable für fpisin*i

            if (options->inf_func == FUNC_FPISIN) // Prüft, ob die Funktion FUNC_FPISIN ist
            {
                fpisin_i = fpisin * sin(pih * (double)i); // Berechnet fpisin*i
            }

            /* over all columns */
            for (j = 1; j < N; j++) // Schleife über die Spalten
            {
                star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]); // Berechnet den neuen Wert

                if (options->inf_func == FUNC_FPISIN) // Prüft, ob die Funktion FUNC_FPISIN ist
                {
                    star += fpisin_i * sin(pih * (double)j); // Fügt den fpisin-Wert hinzu
                }

                if (options->termination == TERM_PREC || term_iteration == 1) // Prüft, ob die Termination PREC ist oder die letzte Iteration
                {
                    residuum = Matrix_In[i][j] - star; // Berechnet das Residuum
                    residuum = (residuum < 0) ? -residuum : residuum; // Nimmt den absoluten Wert des Residuum
                    maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum; // Aktualisiert das maximale Residuum
                }

                Matrix_Out[i][j] = star; // Setzt den neuen Wert in die Ausgangsmatrix
            }
        }

        results->stat_iteration++; // Erhöht die Iterationsanzahl
        results->stat_precision = maxResiduum; // Setzt die Präzision

        /* exchange m1 and m2 */
        i = m1; // Tauscht m1 und m2
        m1 = m2; // Tauscht m1 und m2
        m2 = i; // Tauscht m1 und m2

        /* check for stopping calculation depending on termination method */
        if (options->termination == TERM_PREC) // Prüft, ob die Termination PREC ist
        {
            if (maxResiduum < options->term_precision) // Prüft, ob das maximale Residuum kleiner als die Präzision ist
            {
                term_iteration = 0; // Setzt die Iterationen auf 0
            }
        }
        else if (options->termination == TERM_ITER) // Prüft, ob die Termination ITER ist
        {
            term_iteration--; // Verringert die Iterationen
        }
    }

    results->m = m2; // Setzt die Matrix-Nummer
}
static
void
calculate_row (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
    int i, j;           // Lokale Variablen für Schleifen
    int m1, m2;         // Indizes für alte und neue Matrizen
    double star;        // Viermal der Mittelwert minus 4 Nachbarwerte
    double residuum;    // Residuum der aktuellen Iteration
    double maxResiduum; // Maximales Residuum eines Slaves in der Iteration
    double maxLocalResiduum; // Maximales Residuum eines Slaves in der Iteration in einem Thread

    int const N = arguments->N; // Anzahl der Zeilen
    double const h = arguments->h; // Länge eines Raums zwischen zwei Zeilen

    double pih = 0.0; // Variable für pi*h
    double fpisin = 0.0; // Variable für fpisin

    double** Matrix_Out; // Zeiger auf die Ausgangsmatrix
    double** Matrix_In; // Zeiger auf die Eingabematrix

    int term_iteration = options->term_iteration; // Anzahl der Iterationen

    /* initialize m1 and m2 depending on algorithm */
    if (options->method == METH_JACOBI) // Prüft, ob die Methode Jacobi ist
    {
        m1 = 0; // Setzt m1 auf 0
        m2 = 1; // Setzt m2 auf 1
    }
    else // Wenn die Methode nicht Jacobi ist
    {
        m1 = 0; // Setzt m1 auf 0
        m2 = 0; // Setzt m2 auf 0
    }

    if (options->inf_func == FUNC_FPISIN) // Prüft, ob die Funktion FUNC_FPISIN ist
    {
        pih = PI * h; // Berechnet pi*h
        fpisin = 0.25 * TWO_PI_SQUARE * h * h; // Berechnet fpisin
    }
#pragma omp parallel if(options->method == METH_JACOBI) default(none) num_threads(options->number) \
    private(i, j, star, residuum, maxLocalResiduum) \
    shared(arguments, m1, m2, N, options, pih, fpisin, term_iteration, results, Matrix_Out, Matrix_In, maxResiduum)
{
    while (term_iteration > 0) // Schleife über die Iterationen
    {
        #pragma omp master
        {
         	Matrix_Out  = arguments->Matrix[m1]; // Holt die Ausgangsmatrix
            Matrix_In   = arguments->Matrix[m2]; // Holt die Eingabematrix
            maxResiduum = 0; // Setzt das maximale Residuum auf 0
        }
        maxLocalResiduum = 0; // Setzt das lokale maximale Residuum auf 0
        #pragma omp barrier

        /* over all rows */
        #pragma omp for
        for (i = 1; i < N; i++) // Schleife über die Zeilen
        {
            double fpisin_i = 0.0; // Variable für fpisin*i

            if (options->inf_func == FUNC_FPISIN) // Prüft, ob die Funktion FUNC_FPISIN ist
            {
                fpisin_i = fpisin * sin(pih * (double)i); // Berechnet fpisin*i
            }

            /* over all columns */
            for (j = 1; j < N; j++) // Schleife über die Spalten
            {
                star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]); // Berechnet den neuen Wert

                if (options->inf_func == FUNC_FPISIN) // Prüft, ob die Funktion FUNC_FPISIN ist
                {
                    star += fpisin_i * sin(pih * (double)j); // Fügt den fpisin-Wert hinzu
                }

                if (options->termination == TERM_PREC || term_iteration == 1) // Prüft, ob die Termination PREC ist oder die letzte Iteration
                {
                    residuum = Matrix_In[i][j] - star; // Berechnet das Residuum
                    residuum = (residuum < 0) ? -residuum : residuum; // Nimmt den absoluten Wert des Residuum
                    if (residuum > maxLocalResiduum) {maxLocalResiduum = residuum;} // Aktualisiert das lokale maximale Residuum
                }

                Matrix_Out[i][j] = star; // Setzt den neuen Wert in die Ausgangsmatrix
            }
        }
        #pragma omp atomic compare
        if (maxLocalResiduum > maxResiduum) {maxResiduum = maxLocalResiduum;} // Aktualisiert das maximale Residuum
        #pragma omp barrier
        #pragma omp master
        {

            results->stat_iteration++; // Erhöht die Iterationsanzahl
            results->stat_precision = maxResiduum; // Setzt die Präzision

            /* exchange m1 and m2 */
            i = m1; // Tauscht m1 und m2
            m1 = m2; // Tauscht m1 und m2
            m2 = i; // Tauscht m1 und m2

            /* check for stopping calculation depending on termination method */
            if (options->termination == TERM_PREC) // Prüft, ob die Termination PREC ist
            {
                if (maxResiduum < options->term_precision) // Prüft, ob das maximale Residuum kleiner als die Präzision ist
                {
                    term_iteration = 0; // Setzt die Iterationen auf 0
                }
            }
            else if (options->termination == TERM_ITER) // Prüft, ob die Termination ITER ist
            {
                term_iteration--; // Verringert die Iterationen
            }
        }
        #pragma omp barrier
    }
} /* End #pragma parallel */

    results->m = m2; // Setzt die Matrix-Nummer
}
static
void
calculate_column (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
    int i, j;           // Lokale Variablen für Schleifen
    int m1, m2;         // Indizes für alte und neue Matrizen
    double star;        // Vier
