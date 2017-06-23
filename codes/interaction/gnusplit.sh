#!/bin/bash
########################
#                      #
# Scheduler Directives #
#                      #
########################

### Set the name of the job, where jobname is a unique name for your job
#PBS -N gnusplit_liver

### Select the shell you would like the script to execute within
#PBS -S /bin/bash

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=500:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate four cores on a single node.
#PBS -l nodes=5:ppn=10

### Inform the scheduler of the amount of memory you expect to use.
### Use units of 'b', 'kb', 'mb', or 'gb'
#PBS -l mem=50gb

### Set the destination for your program's output.
#PBS -o $HOME/gnusplit_liver.out
#PBS -e $HOME/gnusplit_liver.err
start=`date +%s`

module load intel/2017
module load gcc/6.2.0
module load R/3.3.2
module load parallel

function do_this_chr {
	if test $# -eq 4
	then
		TISSUE=$1
		CHR=$2
		NAA=$3
		NEA=$4
	else
		echo 'usage: ./gnusplit.sh [tissue_name, chr, naa, nea]'
		exit
	fi
	
	echo starting $TISSUE $CHR
	N=12
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

	PIDS="$(pgrep -P $$)"
	for pid in $PIDS
	do
		wait $pid
	done

        Rscript /group/im-lab/nas40t2/tae/differentialeQTLs/codes/interaction/addup.R $TISSUE $CHR $N
        rm /group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$TISSUE/interaction/chr${CHR}result*.txt        
}

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

	PIDS="$(pgrep -P $$)"
	for pid in $PIDS
	do
		wait $pid
	done

	for CHR in {1..22}
	do
		Rscript /group/im-lab/nas40t2/tae/differentialeQTLs/codes/interaction/addup.R $TISSUE $CHR $N
		rm /group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$TISSUE/interaction/chr${CHR}result*.txt
	done
}

function run {
	chr=$1
	name="Liver"
	var1="/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$name/runlamp/intermediateData/AAsubjects"
	var2="/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/$name/runlamp/intermediateData/EAsubjects"
	naa=$(wc -l < "$var1")
	nea=$(wc -l < "$var2")
	do_this_chr $name $chr $naa $nea
}

export -f do_this_chr
export -f do_this_tissue
export -f run

parallel run ::: {1..22}
end=`date +%s`
runtime=$((end-start))
echo runtime
