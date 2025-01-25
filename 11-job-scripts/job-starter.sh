while read p; do
 sbatch $p
done < $1
