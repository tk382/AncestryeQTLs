#!/bin/bash

########################
#                      #
# Scheduler Directives #
#                      #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N cleanDataForTests_all

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
#PBS -o $HOME/cleanDataForTests_all.out
#PBS -e $HOME/cleanDataForTests_all.err

module load gcc/6.2.0
module load parallel

cd /group/im-lab/nas40t2/tae/differentialeQTLs/codes/dataclean
filename="/group/im-lab/nas40t2/tae/differentialeQTLs/data/tissues"
while read -r line
do
	name="$line"
	for i in {1..22}
	do
		cd /group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$name/data
		chmod 750 global
	done
	echo $name
done < "$filename"
