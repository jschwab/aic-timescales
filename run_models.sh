#!/bin/bash -x
#PBS -N AIC
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -V
#PBS -m abe
#PBS -M jwschwab@ucsc.edu

if [ -n "${PBS_O_WORKDIR}" ]; then
   cd ${PBS_O_WORKDIR}
   module load mesasdk/20180127
   export MESA_DIR=/pfs/jschwab/mesa-r10398
   export OMP_NUM_THREADS=16
else
    # assume you're already there
    PBS_O_WORKDIR=$(pwd)
    PBS_ARRAYID=1
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
    export MESA_INLIST=${PBS_O_WORKDIR}/inlist

    # make a temporary directory
    TMPDIR=$(mktemp -d)
    cd ${TMPDIR}

    # cache locally
    mkdir -p caches
    export MESA_CACHES_DIR=$(pwd)/caches

    # setup inlists via soft links
    ln -sf ${PBS_O_WORKDIR}/inlist_fixed .
    ln -sf ${PBS_O_WORKDIR}/inlists/${INLIST_VARIABLE} inlist_variable
    ln -sf ${PBS_O_WORKDIR}/inlists/${INLIST_IO} inlist_io
    ln -sf ${PBS_O_WORKDIR}/inlist_pgstar .
    
    # softlink in some more stuff
    ln -sf ${PBS_O_WORKDIR}/history_columns.list .
    ln -sf ${PBS_O_WORKDIR}/profile_columns.list .
    ln -sf ${PBS_O_WORKDIR}/models .
    mkdir -p ${PBS_O_WORKDIR}/output/${ID}
    ln -sf ${PBS_O_WORKDIR}/output/${ID} LOGS

    # set up states and transitions
    ln -sf ${PBS_O_WORKDIR}/urca.net .
    ln -sf ${PBS_O_WORKDIR}/weak.net .
    ln -sf ${PBS_O_WORKDIR}/weak.states .
    ln -sf ${PBS_O_WORKDIR}/transitions/${TRANSITIONS_FILE} weak.transitions
    
    # run MESA
    ${PBS_O_WORKDIR}/star

    # make movie
    DATE=$(date +%F)
    images_to_movie.sh 'png/grid1*.png' ${PBS_O_WORKDIR}/movies/${ID}-${DATE}.mp4

}

do_one ${PBS_ARRAYID} AIC.batch
