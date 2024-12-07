#!/bin/bash
#SBATCH -N 5             # Anzahl der Knoten
#SBATCH -n 25             # Anzahl der totalen Tasks (5 Knoten * 5 Prozesse)
#SBATCH --output=timempi.out
#SBATCH --error=timempi.err
#SBATCH --partition=west  # Partition auf "west" festlegen

# Starte timescript.sh auf jedem Knoten mit srun
mpiexec -env I_MPI_PLATFORM 0 ./timempi
