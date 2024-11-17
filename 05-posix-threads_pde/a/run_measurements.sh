#!/bin/bash

# Array mit den Werten für den ersten Parameter
values=(1 2 3 4 5 6 7 8 9 10 11 12)

# Anzahl der Messungen
num_measurements=3

# Name der Ausgabedatei
output_file="output_all.txt"

# Leere die Ausgabedatei, falls sie bereits existiert
> "$output_file"

# Schleife über die Anzahl der Messungen
for ((m=1; m<=num_measurements; m++)); do
    # Schleife über die Werte
    for value in "${values[@]}"; do
        echo "Measurement $m for value $value" >> "$output_file"
        ./partdiff-posix "$value" 2 512 2 2 900 >> "$output_file"
    done
done
