#!/bin/bash

# Time limit is one minute. See "man sbatch" for other time formats.
#SBATCH --time=00:00:10
#SBATCH -N 10
# Use "west" partition.
#SBATCH --partition=west
# Output goes to "job.out", error messages to "job.err".
#SBATCH --output=timescript.out
#SBATCH --error=job.err

srun ./timescript.sh