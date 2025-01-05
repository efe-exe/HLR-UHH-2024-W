#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include </opt/spack/20220821/opt/spack/linux-ubuntu20.04-x86_64/gcc-12.1.0/mpich-4.0.2-ogro7lxzozs6o2djgn6325cqeyzdhil3/include/mpi.h>

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

int per_process (int N, int nprocs) {
    return (N + nprocs - 1) / nprocs; // ceil(N / nprocs)
}

int occupied_count (int N, int nprocs, int rank) {
    int pp = per_process(N, nprocs);
    return (N % nprocs == 0) ? pp : pp - (rank >= (N % nprocs)); // The last few process have one slot not filled in
}

int rank_bound(int N, int nprocs, int rank) {
    int needed = N / nprocs;
    int rem = N % nprocs;
    return needed * rank + min(rank, rem);
}

int* init (int N, int nprocs, int rank)
{
	// TODO
    int per_process = (N + nprocs - 1) / nprocs; // ceil(N / nprocs)
	int* buf = (int*)malloc(sizeof(int) * per_process);

	srand(time(NULL) + rank);

	for (int i = 0; i < occupied_count(N, nprocs, rank); i++)
	{
		// Do not modify "% 12"
		buf[i] = rand() % 12;
	}

	return buf;
}

int circle (int* buf, int pp, int nprocs, int rank, int *oc)
{
    if (nprocs == 1) {
        return 1; // If there is only one process, we would always finish after exactly one iteration.
                  // Since we can't even setup the termination condition without special cases, we don't even try
    }
    int term_value, running, iterations, prev, next;
    prev = (rank - 1 + nprocs) % nprocs;
    next = (rank + 1) % nprocs;
    if (rank == 0) {
        term_value = buf[0];
        MPI_Ssend(&buf[0], 1, MPI_INT, nprocs - 1, 200, MPI_COMM_WORLD);
    } else if (rank == nprocs - 1) {
        MPI_Recv(&term_value, 1, MPI_INT, 0, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    running = 1;
    iterations = 0;
    while (running) {
        int prev_oc;
        MPI_Sendrecv(oc, 1, MPI_INT, next, 301, &prev_oc, 1, MPI_INT, prev, 301, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv_replace(buf, pp, MPI_INT, next, 302, prev, 302, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        *oc = prev_oc;
        // Check for oc here before accessing buf[0] so that we don't read invalid data
        // for the case that nprocs > N
        if (rank == nprocs - 1 && *oc >= 1 && buf[0] == term_value) {
            running = 0;
        }
        MPI_Bcast(&running, 1, MPI_INT, nprocs - 1, MPI_COMM_WORLD);
        iterations++;
    }
	return iterations * 12 + term_value;
}

int main (int argc, char** argv)
{
	int N, rank, nprocs, pp, oc;
	int* buf;
    int ret;

    MPI_Init(&argc, &argv);

	if (argc < 2)
	{
		printf("Arguments error!\nPlease specify a buffer size.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		return EXIT_FAILURE;
	}
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// Array length
	N = atoi(argv[1]);
    if (N < 1)
    {
        printf("Arguments error!\nPlease specify a valid buffer size.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        return EXIT_FAILURE;
    }
    pp = per_process(N, nprocs);
    oc = occupied_count(N, nprocs, rank);
	buf = init(N, nprocs, rank);


    // To gurantee a readable output, the active code gurantees console output in the correct order
    // by doing all IO on rank 0. This doesn't quite fit the task.
    // A version that does, but sometimes fails to be in the correct order is here:
    /*         //Just add a / at the start here to switch to the other version
    if (rank == 0) {
        printf("\nBEFORE\n");
    }
    for (int j = 0; j < nprocs; j++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == j) {
            for (int k = 0; k < oc; k++) {
        	    printf("[Before,rank=%d]: %d\n", j, buf[k]);
            }
        }
    }

/*/
	if (rank == 0) {
	    printf("\nBEFORE\n");
        int i = 0;
	    for (int k = 0; k < oc; k++) {
	        printf("[Before,rank=%d,i=%d]: %d\n", rank, i, buf[k]);
            i++;
	    }

	    for (int j = 1; j < nprocs; j++)
	    {
            int this_oc;
            MPI_Recv(&this_oc, 1, MPI_INT, j, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int k = 0; k < this_oc; k++) {
                int value;
                MPI_Recv(&value, 1, MPI_INT, j, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	            printf("[Before,rank=%d,i=%d]: %d\n", j, i, value);
                i++;
            }
        }
	} else {
        MPI_Ssend(&oc, 1, MPI_INT, 0, 101, MPI_COMM_WORLD);
        for (int k = 0; k < oc; k++) {
            MPI_Ssend(&buf[k], 1, MPI_INT, 0, 102, MPI_COMM_WORLD);
        }
    }
//*/
    MPI_Barrier(MPI_COMM_WORLD); // To help find the phases in the vampir graph
	ret = circle(buf, pp, nprocs, rank, &oc);
    MPI_Barrier(MPI_COMM_WORLD); // To help find the phases in the vampir graph

    /*         //See above
    if (rank == 0) {
        printf("\nAFTER %d iterations, termination value=%d\n", ret / 12, ret % 12);
    }
    for (int j = 0; j < nprocs; j++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == j) {
            for (int k = 0; k < oc; k++) {
                printf("[After,rank=%d]: %d\n", j, buf[k]);
            }
        }
    }
    if (rank == 0) {
        printf("\nAFTER %d iterations, termination value=%d\n", ret / 12, ret % 12);
    }

/*/
    if (rank == 0) {
        printf("\nAFTER %d iterations, termination value=%d\n", ret / 12, ret % 12);
        int i = 0;
        for (int k = 0; k < oc; k++) {
            printf("[After,rank=%d,i=%d]: %d\n", rank, i, buf[k]);
            i++;
        }

        for (int j = 1; j < nprocs; j++)
        {
            int this_oc;
            MPI_Recv(&this_oc, 1, MPI_INT, j, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int k = 0; k < this_oc; k++) {
                int value;
                MPI_Recv(&value, 1, MPI_INT, j, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                printf("[After,rank=%d,i=%d]: %d\n", j, i, value);
                i++;
            }
        }
        printf("\nAFTER %d iterations, termination value=%d\n", ret / 12, ret % 12);
    } else {
        MPI_Ssend(&oc, 1, MPI_INT, 0, 101, MPI_COMM_WORLD);
        for (int k = 0; k < oc; k++) {
            MPI_Ssend(&buf[k], 1, MPI_INT, 0, 102, MPI_COMM_WORLD);
        }
    }
//*/
    MPI_Finalize();

	return EXIT_SUCCESS;
}
