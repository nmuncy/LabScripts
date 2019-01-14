#!/bin/bash





workDir=~/compute/STT_bids
slurmDir=${workDir}/derivatives/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/sttN4_${time}

mkdir -p $outDir


sbatch \
-o ${outDir}/output_sttN4.txt \
-e ${outDir}/error_sttN4.txt \
Task_step4_grpAnalysis.sh
