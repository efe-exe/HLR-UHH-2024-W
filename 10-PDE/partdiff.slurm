#!/bin/bash
#SBATCH -N 5             # Anzahl der Knoten
#SBATCH -n 25             # Anzahl der totalen Tasks (5 Knoten * 5 Prozesse)
#SBATCH --output=partdiff.out
#SBATCH --error=partdiff.err
#SBATCH --partition=west  # Partition auf "west" festlegen

mpiexec -env I_MPI_PLATFORM 0 ./partdiff 1 2 256 2 2 1000
