#!/bin/bash

#SBATCH --time=30:00:00   # walltime
#SBATCH --ntasks=6   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8gb   # memory per CPU core
#SBATCH -J "sttN2"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE




# Written by Nathan Muncy on 10/24/18

# -- Notes, in no particular order
#
# 1) will do deconvolutions (OLS, GLS), post files, and print out info
#
# 2) I highly recommend using the arrays I built in order to keep multiple phases/runs/sub-bricks straight
#
# 3) Behavioral vectors (beh.01.1D) need to be created before this script can run! This is not included in the pipeline
#
# 4) "###???" needs attention
#
# 5) use testMode to save all intermediates
#
# 6) big steps are marked "### --- foobar --- ###", annotations with "#"
#
# 7) Written to show that you can do multiple deconvoutions for multiple phases of the experiment
#		- Obviously, you won't have the same deconvolutions
#
# 8) regArr=(${deconList[@]:0:2}) -> pulls info from array "deconList" starting with position 0 and continuing for 2 positions (0,1) and assigns it to new array "regArr"
#		- Needed in case different phases have a different number of deconvolutions
#
# 9) deepClean to remove as much as possible
#
# 10) This will write multiple deconvolution scripts, and then run them as needed



subj=$1
testMode=$2


###??? update these
workDir=~/compute/STT_new/$subj
phaseArr=(STUDY TEST{1,2})  		# Each PHASE of experiment (e.g. STUDY, TEST1, TEST2). Phase string should match block string, e.g. phase (TEST1) consists of blocks (TEST11 and TEST12)
blockArr=(1 2 4)   					# number of blocks (runs) in each phase
phaseLen=${#phaseArr[@]}





cd $workDir


### --- Motion --- ###
#
# motion and censor files are constructed. Multiple motion files
# include mean and derivative of motion.


c=0; while [ $c -lt $phaseLen ]; do

	phase=${phaseArr[$c]}
	nruns=${blockArr[$c]}

	cat dfile.${phase}*.1D > dfile_rall_${phase}.1D

	if [ ! -s censor_${phase}_combined.1D ]; then

		# files: de-meaned, motion params (per phase)
		1d_tool.py -infile dfile_rall_${phase}.1D -set_nruns $nruns -demean -write motion_demean_${phase}.1D
		1d_tool.py -infile dfile_rall_${phase}.1D -set_nruns $nruns -derivative -demean -write motion_deriv_${phase}.1D
		1d_tool.py -infile motion_demean_${phase}.1D -set_nruns $nruns -split_into_pad_runs mot_demean_${phase}
		1d_tool.py -infile dfile_rall_${phase}.1D -set_nruns $nruns -show_censor_count -censor_prev_TR -censor_motion 0.3 motion_${phase}


		# determine censor
		cat out.cen.${phase}*.1D > outcount_censor_${phase}.1D
		1deval -a motion_${phase}_censor.1D -b outcount_censor_${phase}.1D -expr "a*b" > censor_${phase}_combined.1D
	fi

	let c=$[$c+1]
done




### --- Deconvolve --- ###
#
# A deconvolution (OLS) script is written (for review) and then run.
# Currently, behaviors are called beh_1, beh_2 etc due to practical
# constraints. I'll probably fix this in the future. A deconv script
# is written for each run of the experiment, and for different
# phases a conditional allows control over things like timing
# files and different block times. Timing files are not constructed.
# A REML (GLS) script is written, for future use.


# Function - write deconvolution script
GenDecon (){

	# assign vars for readability
	h_phase=$1
	h_block=$2
	h_tfile=$3
	h_trlen=$4
    h_input=$5
    h_out=$6


	# build motion list
	unset stimBase
    x=1; for ((r=1; r<=${h_block}; r++)); do
        for ((b=0; b<=5; b++)); do

            stimBase+="-stim_file $x mot_demean_${h_phase}.r0${r}.1D'[$b]' -stim_base $x -stim_label $x mot_$x "
            x=$[$x+1]
        done
    done


	# build behavior list
	unset stimBeh
    tBeh=`ls ${h_tfile}* | wc -l`
    for ((t=1; t<=$tBeh; t++)); do

        stimBeh+="-stim_times $x ${h_tfile}.0${t}.1D \"BLOCK(${h_trlen},1)\" -stim_label $x beh_$t "
        x=$[$x+1]
    done


	# num_stimts
    h_nstim=$(($x-1))


	# write script
    echo "3dDeconvolve \
    -input $h_input \
    -censor censor_${h_phase}_combined.1D \
    -polort A -float \
    -num_stimts $h_nstim \
    $stimBase \
    $stimBeh \
    -jobs 6 \
    -bout -fout -tout \
    -x1D X.${h_out}.xmat.1D \
    -xjpeg X.${h_out}.jpg \
    -x1D_uncensored X.${h_out}.nocensor.xmat.1D \
    -errts ${h_out}_errts \
    -bucket ${h_out}_stats" > ${h_out}_deconv.sh
}



### Run multiple decons for each phase
c=0; while [ $c -lt $phaseLen ]; do

    phase=${phaseArr[$c]}

	# create input list
	unset input
	for j in ${phase}*_scale+tlrc.HEAD; do
		input+="${j%.*} "
	done


	###??? update number of conditionals to equal $phaseLen
	if [ $phase == ${phaseArr[0]} ]; then

		# variables/arrays for various deconvolutions from STUDY phase (SpT1, SpT1pT2)						###??? update these in each block
		trialLen=3																							# stim length (to be modeled)
		tfile=(timing_files/{Study_pred_Test1_TF_4_behVect,Study_pred_Test1_pred_Test2_TF_4_behVect})		# array of timing files
		out=(SpT1 SpT1pT2)																					# array of prefix strings
		outLen=${#out[@]}

		d=0; while [ $d -lt $outLen ]; do

			# write script
			GenDecon $phase ${blockArr[$c]} ${tfile[$d]} $trialLen "$input" ${out[$d]}

			# run script
			if [ ! -f ${out}_stats+tlrc.HEAD ]; then
				source ${out}_deconv.sh
			fi
			let d=$[$d+1]
		done

	elif [ $phase == ${phaseArr[1]} ]; then

		trialLen=3
		tfile=(timing_files/{Test1_TF_4_behVect,Test1_pred_Test2_TF_4_behVect})
		out=(T1 T1pT2)
		outLen=${#out[@]}

		d=0; while [ $d -lt $outLen ]; do

			GenDecon $phase ${blockArr[$c]} ${tfile[$d]} $trialLen "$input" ${out[$d]}

			if [ ! -f ${out}_stats+tlrc.HEAD ]; then
				source ${out}_deconv.sh
			fi
			let d=$[$d+1]
		done

	elif [ $phase == ${phaseArr[2]} ]; then

		trialLen=4.5
		tfile=(timing_files/{Test2_TF_4_behVect,Test2_ref_Test1_TF_4_behVect})
		out=(T2 T2fT1)
		outLen=${#out[@]}

		d=0; while [ $d -lt $outLen ]; do

			GenDecon $phase ${blockArr[$c]} ${tfile[$d]} $trialLen "$input" ${out[$d]}

			if [ ! -f ${out}_stats+tlrc.HEAD ]; then
				source ${out}_deconv.sh
			fi
			let d=$[$d+1]
		done
	fi
	let c=$[$c+1]
done


# determine number of decons conducted
c=0; for i in *_stats+tlrc.HEAD; do
	deconList[$c]=${i%_*}
	let c=$[$c+1]
done




#### --- REML and Post Calcs --- ###
#
# REML deconvolution (GLS) is run, excluding WM signal. REML will
# probably become standard soon, so I'll get this working at some point.
# Global SNR and corr are calculated.


for i in ${phaseArr[@]}; do

	# align phase w/decon output
	if [ $i == ${phaseArr[0]} ]; then
		regArr=(${deconList[@]:0:2})  		###??? update which decons match for each phase
	elif [ $i == ${phaseArr[1]} ]; then
		regArr=(${deconList[@]:2:2})
	elif [ $i == ${phaseArr[2]} ]; then
		regArr=(${deconList[@]:4:2})
	fi


	# all runs signal
	count=`1d_tool.py -infile censor_${i}_combined.1D -show_trs_uncensored encoded`

	if [ ! -f ${regArr[0]}_TSNR+tlrc.HEAD ]; then

		3dTcat -prefix tmp_${i}_all_runs ${i}*_scale+tlrc.HEAD
		3dTstat -mean -prefix tmp_${i}_allSignal tmp_${i}_all_runs+tlrc"[${count}]"
	fi


	## timeseries of eroded WM - ##### Fix here
	#if [ ! -f tmp_allRuns_${i}_WMe+tlrc.HEAD ] && [ ! -f ${regArr[0]}_stats_REML+tlrc.HEAD ]; then

		#3dTcat -prefix tmp_allRuns_${i} ${i}*_volreg_clean+tlrc.HEAD
		#3dcalc -a tmp_allRuns_${i}+tlrc -b final_mask_WM_eroded+tlrc -expr "a*bool(b)" -datum float -prefix tmp_allRuns_${i}_WMe
		#3dmerge -1blur_fwhm 20 -doall -prefix ${i}_WMe_rall tmp_allRuns_${i}_WMe+tlrc
	#fi


	for j in ${regArr[@]}; do

		# kill if decon failed
		if [ ! -f ${j}_stats+tlrc.HEAD ]; then
			echo "Decon failed on $j ..." >&2
			exit 1
		fi


		###### to do: fix this section

		## REML
		#if [ ! -f ${j}_stats_REML+tlrc.HEAD ]; then
			#tcsh -x ${j}_stats.REML_cmd -dsort ${i}_WMe_rall+tlrc
		#fi


		## kill if REMl failed
		#if [ ! -f ${j}_stats_REML+tlrc.HEAD ]; then
			#echo "REML failed on $j ..." >$2
			#exit 2
		#fi


		## calc SNR, corr
		#if [ ! -f ${j}_TSNR+tlrc.HEAD ]; then

			#3dTstat -stdev -prefix tmp_${j}_allNoise ${j}_errts_REML+tlrc"[${count}]"

			#3dcalc -a tmp_${i}_allSignal+tlrc \
			#-b tmp_${j}_allNoise+tlrc \
			#-c full_mask+tlrc \
			#-expr 'c*a/b' -prefix ${j}_TSNR

			#3dTnorm -norm2 -prefix tmp_${j}_errts_unit ${j}_errts_REML+tlrc
			#3dmaskave -quiet -mask full_mask+tlrc tmp_${j}_errts_unit+tlrc > ${j}_gmean_errts_unit.1D
			#3dcalc -a tmp_${j}_errts_unit+tlrc -b ${j}_gmean_errts_unit.1D -expr 'a*b' -prefix tmp_${j}_DP
			#3dTstat -sum -prefix ${j}_corr_brain tmp_${j}_DP+tlrc
		#fi


		# detect pairwise cor
		1d_tool.py -show_cormat_warnings -infile X.${j}.xmat.1D | tee out.${j}.cormat_warn.txt
	done
done



for i in ${deconList[@]}; do

	# sum of regressors, stim only x-matrix
	if [ ! -s X.${i}.stim.xmat.1D ]; then

		reg_cols=`1d_tool.py -infile X.${i}.nocensor.xmat.1D -show_indices_interest`
		3dTstat -sum -prefix ${i}_sum_ideal.1D X.${i}.nocensor.xmat.1D"[$reg_cols]"
		1dcat X.${i}.nocensor.xmat.1D"[$reg_cols]" > X.${i}.stim.xmat.1D
	fi
done




#### --- Print out info, Clean --- ###
#
# Print out information about the data - of particulatr interest
# is the number of TRs censored, which step4 uses. This involves
# producing the needed files, generating a set of review scripts,
# and then running my favorite. Many intermediates get removed.


# organize files for what gen*py needs
3dcopy full_mask+tlrc full_mask.${subj}+tlrc

for i in ${phaseArr[@]}; do

	cat outcount.${i}*.1D > outcount_all_${i}.1D

	c=1; for j in ${i}*+orig.HEAD; do

		file=${j%.*}
		prefix=${j%%_*}
		3dcopy $file pb00.${subj}.r0${c}.tcat
		3dcopy ${prefix}_volreg_clean+tlrc pb02.${subj}.r0${c}.volreg

		let c=$[$c+1]
	done


	for k in ${deconList[@]}; do

		# a touch more organization (gen*py is very needy)
		dset=${k}_stats+tlrc
		cp X.${k}.xmat.1D X.xmat.1D
		3dcopy ${k}_errts+tlrc errts.${subj}+tlrc


		# generate script
		gen_ss_review_scripts.py \
		-subj ${subj} \
		-rm_trs 0 \
		-motion_dset dfile_rall_${i}.1D \
		-outlier_dset outcount_all_${i}.1D \
		-enorm_dset  motion_${i}_enorm.1D \
		-mot_limit 0.3 \
		-out_limit 0.1 \
		-xmat_regress X.${k}.xmat.1D \
		-xmat_uncensored X.${k}.nocensor.xmat.1D \
		-stats_dset ${dset} \
		-final_anat final_anat+tlrc \
		-final_view tlrc \
		-exit0


		# run script - write an output for e/analysis
		./\@ss_review_basic | tee out_summary_${k}.txt


		# clean
		rm errts.*
		rm X.xmat.1D
		rm pb0*
		rm *ss_review*
	done
done



# clean
if [ $testMode == 1 ]; then
	rm tmp_*
	rm -r a*
	rm final_mask_{CSF,GM}*
	rm *corr_brain*
	rm *gmean_errts*
	rm *volreg*
	rm Temp*
	rm *WMe_rall*
	rm full_mask.*
fi
