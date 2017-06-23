#!/bin/bash

########################
#                      #
# Scheduler Directives #
#                      #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N pca

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=100:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate four cores on a single node.
#PBS -l procs=10

### Inform the scheduler of the amount of memory you expect to use.
### Use units of 'b', 'kb', 'mb', or 'gb'
#PBS -l mem=100gb

### Set the destination for your program's output.
#PBS -o $HOME/pca.out
#PBS -e $HOME/pca.err

module restore
module load gcc/6.2.0
module load parallel

cd /group/im-lab/nas40t2/tae/differentialeQTLs/codes/dataclean
function dopca {
	name="$1"
	var1="/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$name/runlamp/intermediateData/AAsubjects"
	var2="/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$name/runlamp/intermediateData/EAsubjects"
        var3="/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$name/runlamp/bychr/chr1/CEUgeno.shortandfat.final"
        var4="/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$name/runlamp/bychr/chr1/CHBgeno.shortandfat.final"
        var5="/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$name/runlamp/bychr/chr1/YRIgeno.shortandfat.final"
	naa=$(wc -l < "$var1")
	nea=$(wc -l < "$var2")
	nceu=$(wc -l < "$var3")
	nchb=$(wc -l < "$var4")
	nyri=$(wc -l < "$var5")
	Rscript pca.R $name $naa $nea $nceu $nchb $nyri
	echo $name
}
export -f dopca

#parallel -j5 dopca :::: /group/im-lab/nas40t2/tae/differentialeQTLs/data/tissues


filename="/group/im-lab/nas40t2/tae/differentialeQTLs/data/tissues"
while read -r line
do
       name="$line"
       dopca $name
done < "$filename"
