#define _DEFAULT_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include </opt/spack/20220821/opt/spack/linux-ubuntu20.04-x86_64/gcc-12.1.0/mpich-4.0.2-ogro7lxzozs6o2djgn6325cqeyzdhil3/include/mpi.h>

int main(int argc, char *argv[]) {

    struct timeval tv;
    time_t time;
    int micro_sec = 0, max_micro, min_micro;
    char time_string[30];
    char output[80];
    char hostname[30];
    int rank, size, length, i;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (rank != size - 1) {
        gettimeofday(&tv, NULL);
        gethostname(hostname, 30);

        time = tv.tv_sec;
        micro_sec = tv.tv_usec;

        strftime(time_string, 30, "%Y-%m-%d %T", localtime(&time));
        length = snprintf(output, 80, "[%d] %s @ %s.%d", rank, hostname, time_string, (int)micro_sec);

        MPI_Send(output, length + 1, MPI_CHAR, size - 1, 1, MPI_COMM_WORLD);
        MPI_Reduce(&micro_sec, NULL, 1, MPI_INT, MPI_MIN, size - 1, MPI_COMM_WORLD);
        MPI_Reduce(&micro_sec, NULL, 1, MPI_INT, MPI_MAX, size - 1, MPI_COMM_WORLD);
    } else {
        for (i = 0; i < size - 1; i++) {
            MPI_Recv(output, sizeof(output), MPI_CHAR, i, 1, MPI_COMM_WORLD, &status);
            printf("%s\n", output);
        }
        min_micro = INT_MAX;
        MPI_Reduce(MPI_IN_PLACE, &min_micro, 1, MPI_INT, MPI_MIN, size - 1, MPI_COMM_WORLD);
        max_micro = INT_MIN;
        MPI_Reduce(MPI_IN_PLACE, &max_micro, 1, MPI_INT, MPI_MAX, size - 1, MPI_COMM_WORLD);
        printf("[%d] Kleinster MS-Anteil: %d\n", rank, min_micro);
        printf("[%d] Größte Differenz: %d\n", rank, max_micro - min_micro);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("[%d] beendet jetzt\n", rank);
    MPI_Finalize();
    return 0;
}