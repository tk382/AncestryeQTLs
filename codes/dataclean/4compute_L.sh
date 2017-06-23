#!/bin/bash

########################
#                      #
# Scheduler Directives #
#                      #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N compute_L

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=8:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate four cores on a single node.
#PBS -l procs=10

### Inform the scheduler of the amount of memory you expect to use.
### Use units of 'b', 'kb', 'mb', or 'gb'
#PBS -l mem=20gb

### Set the destination for your program's output.
#PBS -o $HOME/compute_L.out
#PBS -e $HOME/compute_L.err

module load intel/2017
module load gcc/6.2.0
module load R/3.3.2
module load parallel

cd /group/im-lab/nas40t2/tae/differentialeQTLs/codes/dataclean
function computeL {
	tissue="$1"
	var1="/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$tissue/runlamp/intermediateData/AAsubjects"
        var2="/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$tissue/runlamp/intermediateData/EAsubjects"
	naa=$(wc -l < "$var1")
        nea=$(wc -l < "$var2")
        Rscript compute_L.R $tissue $naa $nea	
}

export -f computeL

parallel computeL :::: /group/im-lab/nas40t2/tae/differentialeQTLs/data/tissues
