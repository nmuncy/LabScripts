#!/bin/bash


# stderr and stdout are written to ${outDir}/error_* and ${outDir}/output_* for troubleshooting.
# job submission output are time stamped for troubleshooting

workDir=~/compute/Context_bids   ###??? update this
scriptDir=${workDir}/code
slurmDir=${workDir}/derivatives/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/con1_${time}

mkdir -p $outDir

cd $workDir
for i in sub*; do

    sbatch \
    -o ${outDir}/output_con1_${i}.txt \
    -e ${outDir}/error_con1_${i}.txt \
    ${scriptDir}/Task_step1_sbatch_preproc.sh $i

    sleep 1
done
