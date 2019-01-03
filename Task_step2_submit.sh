#!/bin/bash





workDir=~/compute/STT_new
scriptDir=${workDir}/Scripts
slurmDir=${workDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/sttN2_${time}

mkdir -p $outDir

cd $workDir
for i in s*; do

	# keep s1949 intermediates
	[ $i == s1949 ]; test=$?

    sbatch \
    -o ${outDir}/output_sttN2_${i}.txt \
    -e ${outDir}/error_sttN2_${i}.txt \
    ${scriptDir}/Task_step2_sbatch_regress_condensed.sh $i $test

    sleep 1
done
