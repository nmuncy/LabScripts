#!/bin/bash


# stderr and stdout are written to ${outDir}/error_* and ${outDir}/output_* for troubleshooting.
# job submission output are time stamped for troubleshooting

workDir=~/compute/STT_new   ###??? update this
scriptDir=${workDir}/Scripts
slurmDir=${workDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/sttN1_${time}

mkdir -p $outDir

cd $workDir
for i in s*; do

    sbatch \
    -o ${outDir}/output_sttN1_${i}.txt \
    -e ${outDir}/error_sttN1_${i}.txt \
    ${scriptDir}/Task_step1_sbatch_preproc.sh $i

    sleep 1
done
