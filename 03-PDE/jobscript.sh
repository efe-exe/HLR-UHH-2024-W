#!/bin/bash
#SBATCH -N 5             # Anzahl der Knoten
#SBATCH -n 25            # Anzahl der totalen Tasks (5 Knoten * 5 Prozesse)
#SBATCH -o timescript.out  # Ausgabe in timescript.out speichern
#SBATCH --partition=west  # Partition auf "west" festlegen

# Starte timescript.sh auf jedem Knoten mit srun
srun ./timescript.sh

# Speichert "Fertig" in jobscript.out
echo "Fertig" > jobscript.out

