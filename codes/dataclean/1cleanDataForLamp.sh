#!/bin/bash

########################
#                      #
# Scheduler Directives #
#                      #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N cleanDataForLamp

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=100:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate four cores on a single node.
#PBS -l nodes=1:ppn=12

### Inform the scheduler of the amount of memory you expect to use.
### Use units of 'b', 'kb', 'mb', or 'gb'
#PBS -l mem=100gb

### Set the destination for your program's output.
#PBS -o $HOME/cleanDataForLamp.out
#PBS -e $HOME/cleanDataForLamp.err

module restore
cd /group/im-lab/nas40t2/tae/differentialeQTLs/codes/dataclean
Rscript cleanDataForLamp.R Vagina &
Rscript cleanDataForLamp.R Whole_Blood




