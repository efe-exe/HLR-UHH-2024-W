#define _DEFAULT_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    struct timeval tv;
    time_t time;
    int micro_sec;
    char time_string[30];
    char output[80];
    char hostname[30];
    int micro_secs[size - 1];

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank != size - 1) {
        gettimeofday(&tv, NULL);
        gethostname(hostname, 30);

        time = tv.tv_sec;
        micro_sec = tv.tv_usec;

        strftime(time_string, 30, "%Y-%m-%d %T", localtime(&time));
        snprintf(output, 80, "[%d] %s @ %s.%d", rank, hostname, time_string, (int)micro_sec);

        MPI_Send(output, 80, MPI_CHAR, size - 1, 0, MPI_COMM_WORLD);
        MPI_Send(&micro_sec, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == size - 1) {
        for (int i = 0; i < size - 1; i++) {
            MPI_Recv(output, 80, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&micro_secs[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("%s\n", output);
        }

        int min_micro_sec = micro_secs[0];
        int max_micro_sec = micro_secs[0];
        for (int i = 1; i < size - 1; i++) {
            if (micro_secs[i] < min_micro_sec) {
                min_micro_sec = micro_secs[i];
            }
            if (micro_secs[i] > max_micro_sec) {
                max_micro_sec = micro_secs[i];
            }
        }
        int max_diff = max_micro_sec - min_micro_sec;

        printf("[%d] Kleinster MS-Anteil: %d\n", rank, min_micro_sec);
        printf("[%d] Größte Differenz: %d\n", rank, max_diff);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    printf("[%d] beendet jetzt!\n", rank);

    MPI_Finalize();
    return 0;
}