#!/bin/bash

########################
#                      #
# Scheduler Directives #
#                      #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N genotype_perm_7-12+1

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=100:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate four cores on a single node.
#PBS -l nodes=1:ppn=12

### Inform the scheduler of the amount of memory you expect to use.
### Use units of 'b', 'kb', 'mb', or 'gb'
#PBS -l mem=50gb

### Set the destination for your program's output.
#PBS -o $HOME/genotype_perm7-12+1.out
#PBS -e $HOME/genotype_perm7-12+1.err

module restore
cd /group/im-lab/nas40t2/tae/differentialeQTLs/muscle_skeletal/codes/interaction/

Rscript genotype_permutation.R 7 &
Rscript genotype_permutation.R 8 &
Rscript genotype_permutation.R 9 &
Rscript genotype_permutation.R 10 &
Rscript genotype_permutation.R 11 &
Rscript genotype_permutation.R 12 &
Rscript genotype_permutation.R 1 


