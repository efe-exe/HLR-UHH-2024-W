#!/bin/bash
# SCRIPT TO GATHER RESULTS FROM VARIOUS partdiff-par RUNS
# Cedrick Ansorge 2010/12/04
# cedrick@gmx.net
# Michael Kuhn 2013/01/04
#
# DEFINE PRECISION OF OUTPUT
FPREC=4
#HELP STRING:

usage="
SYNOPSIS:

./gather_results <case_file>

Gather time information from partdiff-par - output.
Must be called from the path, where results are located. For
each case, that statistics are written to a file name <case_
name>.dat where <case_name> is the name given to the file in
<case_file> and used for the name of the output file of
partdiff-par.

If for some reason a value cannot be found in a file, or a
file is does not exist, it is set to '-1' (useful for gnuplot:
set datafile missing "-1")

Output is written to files named <case_name>.dat, where <case_name>
is the name given to each case in file <case_file>.

ARGUMENTS:
      <case_file>   is the file, that was used to generate job-
                    scripts with SLURM_generator.sh
"

# Print help and quit if requested
if [[ $1 == "-h" || $1 == "--help" ]]; then
    echo "${usage}"
    exit
fi

# Quit if too few arguments
if [[ ${BASH_ARGC} -lt 1 ]]; then
    echo "ERROR: Argument missing"
    echo "${usage}"
    exit
fi

# Quit if case_file is given, but does not exist
if [[ ! -f $1 ]];then
    echo "ERROR: $1 is not a file"
    echo "${usage}"
    exit
fi

#########################################
# SETUP
#########################################

name_col=1  # column number with case filename
case_col=8  # column number with configurations for each case
sep=";"     # column separator in case_file

case_file=$1
line_count=`wc -l $case_file | cut -d" " -f1`
line_count=`echo "$line_count + 1" | bc`

########################################
# LOOP OVER CASES
########################################

i=1   # Ignore first line in <case_file>
      # (this is the DUMMY_NAME line for documentation purposes)
while [ $i -le $line_count ]; do
    line=`head -$i ${case_file} | tail -n 1 `
    i=`echo "${i} + 1" | bc`

    case=`echo $line | cut -d"${sep}" -f${name_col} | sed -e's/ //g'`
    conf=`echo $line | cut -d"${sep}" -f${case_col} `

    # WRITE Header to ${case}.dat
    header=`echo '# NPROCS NNODES ILINES TIME'`
    cat <<EOF>${case}.dat
${header}
EOF

    ################################################
    # Average time for each configuration
    ################################################
    for c in $conf; do
	echo $c
	nodes=`echo $c | cut -d"-" -f1`
	procs=`echo $c | cut -d"-" -f2`
	ilins=`echo $c | cut -d"-" -f3`
	fname="${case}_${nodes}_${procs}_${ilins}.out"
	times=`grep Berechnungszeit ${fname} | cut -d":" -f2 | sed -e's/Sekunden/ /g' | sed -e's/S/ /g' | sed -e's/s/ /g'`

	sum=0.;	count=0; meant=-1
	for t in $times; do
	    sum=`echo "scale=${FPREC}; ${sum} + ${t}" | bc`
	    count=`echo "${count} + 1" | bc`
	done
	if [[ $sum > 0. ]]; then
	    meant=`echo "scale=${FPREC}; ${sum} / ${count}" | bc`
	fi

	data="${procs} ${nodes} ${ilins} ${meant}"
	cat <<EOF>>${case}.dat
${data}
EOF
	echo ${data}
    done
done
