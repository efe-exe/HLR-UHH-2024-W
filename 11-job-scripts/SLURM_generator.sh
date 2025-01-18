#!/bin/bash
# SCRIPT TO GENERATE SLURM JOBSCRIPTS
# Cedrick Ansorge 2010/12/01
# cedrick@gmx.net
# Michael Kuhn 2013/01/04
#
# SYNOPSIS:
# ./SLURM_generator.sh
#
# OPTIONS:
# -c <case_file_name>     default ./cases
# -o <output_directory>   default ./
# -x <path_to_executable> default ./partdiff-par
#
#
# STRUCTURE of case_file
#   - one line per case
#   - 7 columns per lines, seperated by ':'
#     1 - string:case_name
#     2 - either 1 or 2 (ALGO_MODE:1 - Gauss-Seidel, 2- Jacobi)
#     3 - either 1 or 2 (TERM_MODE:1 - Precision, 2 - Iteration)
#     4 - if TERM_MODE==1: float (precision to quit)
#         if TERM_MODE==2: integer (iteration to quit)
#     5 - either 1 or 2 (BDY_FUNC: 1 - 0/linear at bdys, 2 - sinusodial)
#     6 - number of samples for each configuration
#     7 - wall clock time
#     8 - List of configurations of the form K-P-N
#         K: integer number of nodes to use
#         P: integer number of processors to use
#         N: integer number of interlines to use
usage="Usage: ./SLURM_generator.sh <case_file_name> <executable_file> <output_directory>"

if [[ $1 -eq "--help" || $1 -eq "-help" ]]; then
    echo ${usage}
    exit
fi

if [[ ${BASH_ARGC} -ne 3 ]]; then
    echo "wrong number of arguments (${BASH_ARGC})"
    echo ${usage}
    exit
fi

##############################
#           SETUP
##############################
MPI_CALL="mpiexec"
SEP=";"
NAME_COL=1
ALGO_COL=2
TERM_COL=3
TERV_COL=4
BDYF_COL=5
SMPL_COL=6
TIME_COL=7
CONF_COL=8

##############################
#    PARSE INPUT
##############################

CASE_FILE="$1"
EXE_FILE="$2"
OUT_DIR="$3"

CASE_COUNT=`wc -l ${CASE_FILE} | sed -e"s/${CASE_FILE}//g"| sed -e"s/ //g"`

# IGNORE FIRST LINE - IS USED FOR DOCUMENTATION PURPOSES (DUMMY_NAME)
LINE_COUNT="1"

while [ ${LINE_COUNT} -le ${CASE_COUNT} ]; do
    LINE_COUNT=`echo "${LINE_COUNT} + 1" | bc`
    LINE=`head -${LINE_COUNT} ${CASE_FILE} | tail -1`
    NAME=`echo ${LINE} | cut -d${SEP} -f${NAME_COL} | sed -e"s/ //g"`
    ALGO=`echo ${LINE} | cut -d${SEP} -f${ALGO_COL} | sed -e"s/ //g"`
    TERM=`echo ${LINE} | cut -d${SEP} -f${TERM_COL} | sed -e"s/ //g"`
    TERV=`echo ${LINE} | cut -d${SEP} -f${TERV_COL} | sed -e"s/ //g"`
    BDYF=`echo ${LINE} | cut -d${SEP} -f${BDYF_COL} | sed -e"s/ //g"`
    SMPL=`echo ${LINE} | cut -d${SEP} -f${SMPL_COL} | sed -e"s/ //g"`
    TIME=`echo ${LINE} | cut -d${SEP} -f${TIME_COL} | sed -e"s/ //g"`
    #do not strip off blanks here as they are used for separation
    CONF=`echo ${LINE} | cut -d${SEP} -f${CONF_COL}`

    for c in ${CONF}; do
	NODE_COUNT=`echo ${c}| cut -d'-' -f1`
	PROC_COUNT=`echo ${c}| cut -d'-' -f2`
	ILIN_COUNT=`echo ${c}| cut -d'-' -f3`
	PPN_COUNT=`echo "${PROC_COUNT} / ${NODE_COUNT}" | bc`
	CHK_COUNT=`echo "${PPN_COUNT}*${NODE_COUNT}" | bc`
	JOB_SCRIPT_NAME="${NAME}_${NODE_COUNT}_${PROC_COUNT}_${ILIN_COUNT}"

	if [ ${CHK_COUNT} -lt ${PROC_COUNT} ]; then
	    PPN_COUNT=`echo "${PPN_COUNT}+1" | bc`
	fi

	ARGS="1 ${ALGO} ${ILIN_COUNT} ${BDYF} ${TERM} ${TERV}"

	echo "${JOB_SCRIPT_NAME}-Args: ${ARGS}"
#	echo "Nodes: ${NODE_COUNT}"
#	echo "Procs: ${PROC_COUNT} (${CHK_COUNT})"
#	echo "Iline: ${ILIN_COUNT}"
#	echo "PPern: ${PPN_COUNT}"


	cat <<EOF>${OUT_DIR}/${JOB_SCRIPT_NAME}.job
#!/bin/bash

#SBATCH --time=${TIME}
#SBATCH --partition=west
#SBATCH --nodes=${NODE_COUNT} --tasks-per-node=${PPN_COUNT}
#SBATCH --error=${JOB_SCRIPT_NAME}.err --output=${JOB_SCRIPT_NAME}.out

EOF
	SMPL_COUNT=0
	while [ ${SMPL_COUNT} -lt ${SMPL} ]; do
	    cat <<EOF>>${OUT_DIR}/${JOB_SCRIPT_NAME}.job
${MPI_CALL} -n ${PROC_COUNT} ${EXE_FILE} ${ARGS}
EOF
	    SMPL_COUNT=`echo "${SMPL_COUNT} + 1" | bc`
	done
    done

done

exit
