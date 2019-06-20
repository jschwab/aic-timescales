#!/bin/bash -x
#SBATCH --job-name=AIC
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=40G
#SBATCH --export=ALL
#SBATCH --time=08:00:00
#SBATCH --array 1-6
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jwschwab@ucsc.edu
#SBATCH -A csc116

if [ -n "${SLURM_SUBMIT_DIR}" ]; then
   cd ${SLURM_SUBMIT_DIR}
   module load mesasdk/20190503
   export MESA_DIR=${PROJECT}/mesa-r10398
   export MESA_CACHES_DIR=/oasis/scratch/comet/jschwab/temp_project/mesa-r10398/cache
   export OMP_NUM_THREADS=8
else
    # assume you're already there
    SLURM_SUBMIT_DIR=$(pwd)
    SLURM_ARRAY_TASK_ID=1
fi

# rebuild MESA
# ./clean
#./mk

# make sure there's a place for the movies and output
mkdir -p movies
mkdir -p output

# define a function to run a named model
do_one() {

    # get the relevant line from the batch file
    read ID INLIST_IO INLIST_VARIABLE TRANSITIONS_FILE <<< $(sed "${1}q;d" < $2)

    # use the main inlist in the submit directory
    export MESA_INLIST=${SLURM_SUBMIT_DIR}/inlist

    # make a temporary directory
    TMPDIR=$(mktemp -d)
    cd ${TMPDIR}

    # setup inlists via soft links
    ln -sf ${SLURM_SUBMIT_DIR}/inlist_fixed .
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/${INLIST_VARIABLE} inlist_variable
    ln -sf ${SLURM_SUBMIT_DIR}/inlists/${INLIST_IO} inlist_io
    ln -sf ${SLURM_SUBMIT_DIR}/inlist_pgstar .
    
    # softlink in some more stuff
    ln -sf ${SLURM_SUBMIT_DIR}/history_columns.list .
    ln -sf ${SLURM_SUBMIT_DIR}/profile_columns.list .
    ln -sf ${SLURM_SUBMIT_DIR}/models .
    mkdir -p ${SLURM_SUBMIT_DIR}/output/${ID}
    ln -sf ${SLURM_SUBMIT_DIR}/output/${ID} LOGS

    # set up states and transitions
    ln -sf ${SLURM_SUBMIT_DIR}/urca.net .
    ln -sf ${SLURM_SUBMIT_DIR}/weak.net .
    ln -sf ${SLURM_SUBMIT_DIR}/weak.states .
    ln -sf ${SLURM_SUBMIT_DIR}/transitions/${TRANSITIONS_FILE} weak.transitions
    
    # run MESA
    ${SLURM_SUBMIT_DIR}/star

    # make movie
    DATE=$(date +%F)
    images_to_movie 'png/grid1*.png' ${SLURM_SUBMIT_DIR}/movies/${ID}-${DATE}.mp4

    cd -
    rm -rf ${TMPDIR}

    
}

do_one ${SLURM_ARRAY_TASK_ID} AIC.batch
