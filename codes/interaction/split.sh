#!/bin/bash
########################
#                      #
# Scheduler Directives #
#                      #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N split_bcor

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=500:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate four cores on a single node.
#PBS -l nodes=1:ppn=8

### Inform the scheduler of the amount of memory you expect to use.
### Use units of 'b', 'kb', 'mb', or 'gb'
#PBS -l mem=100gb

### Set the destination for your program's output.
#PBS -o $HOME/split_bcor.out
#PBS -e $HOME/split_bcor.err

start=`date +%s`

module load intel/2017
module load gcc/6.2.0
module load R/3.3.2
module load parallel

function do_this_tissue {
	if test $# -eq 3
	then
		TISSUE=$1
		NAA=$2
		NEA=$3
	else
		echo 'usage: ./split.sh [tissue_name]'
		exit
	fi
	echo starting $TISSUE
	N=12
	for CHR in {1..22}
	do	
		echo starting chr $CHR
		cd /group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$TISSUE/data/bychr/chr$CHR
		OLINES="$(wc -l exp | awk '{print $1}')"
		LINES="$((${OLINES}-1))"
		OSTEP=$((${LINES}/${N}))
		for ((i=1; i<$N; i+=1))
		do
			OSKIP=$((1 + (${i}-1)*${OSTEP}))
			Rscript /group/im-lab/nas40t2/tae/differentialeQTLs/codes/interaction/interaction.R $TISSUE $CHR $OSKIP $OSTEP $i $NAA $NEA &
		done
		OSKIP=$((1+ (${i}-1)*$OSTEP))
		Rscript /group/im-lab/nas40t2/tae/differentialeQTLs/codes/interaction/interaction.R $TISSUE $CHR $OSKIP -1 $N $NAA $NEA &
	done
	echo waiting...
	PIDS="$(pgrep -P $$)"
	for pid in $PIDS
	do
		wait $pid
		echo $pid done!
	done
	for CHR in {1..22}
	do
		Rscript /group/im-lab/nas40t2/tae/differentialeQTLs/codes/interaction/addup.R $TISSUE $CHR $N &
		rm /group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$TISSUE/interaction/chr${CHR}result*.txt
	done
}

export -f do_this_tissue

name=Brain_Cortex
var1="/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$name/runlamp/intermediateData/AAsubjects"
var2="/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$name/runlamp/intermediateData/EAsubjects"
naa=$(wc -l < "$var1")
nea=$(wc -l < "$var2")
do_this_tissue $name $naa $nea

end=`date +%s`

runtime=$((end-start))
echo runtime


