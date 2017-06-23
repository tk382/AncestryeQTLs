#!/bin/bash

########################
#                      #
# Scheduler Directives #
#                      #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N lampforall

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=10:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate four cores on a single node.
#PBS -l nodes=1:ppn=4

### Inform the scheduler of the amount of memory you expect to use.
### Use units of 'b', 'kb', 'mb', or 'gb'
#PBS -l mem=100gb

### Set the destination for your program's output.
#PBS -o $HOME/lampforall.out
#PBS -e $HOME/lampforall.err
module load intel/2017
module load R/3.3.2
module load gcc/6.2.0
module load parallel
function runlamp {
name="$1"
for i in {1..22}
do
	cd /group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$name/runlamp/bychr/chr${i}
	/group/im-lab/nas40t2/tae/lamproot/bin/lamp config.txt
	/group/im-lab/nas40t2/tae/lamproot/bin/generategraph.sh finalAncestry.txt
	chmod 770 index.html
	chmod 770 tmp
	cd tmp
	chmod 770 *
	echo ${i}
done

}

export -f runlamp

parallel runlamp :::: /group/im-lab/nas40t2/tae/differentialeQTLs/data/tissues
