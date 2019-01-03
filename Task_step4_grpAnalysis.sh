#!/bin/bash

#SBATCH --time=40:00:00   # walltime
#SBATCH --ntasks=10   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8gb   # memory per CPU core
#SBATCH -J "sttN4"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE





# Written by Nathan Muncy on 11/2/18


###--- Notes, in no particular order
#
# 1) does the new ETAC multiple comparison method
#
# 2) constructs a gray matter intersection mask
#		depends on ACT priors
#
# 3) accepts multiple blurs and p-values (multiple should be used)
#
# 4) generates subj lists for 3dttest++, excluding subjects who moved too much (ref step3)
#
# 5) each etac run takes about 5 hours, so do 5*$etacLen for walltime
#
# 6) I'll update this in the future for use with REML files
#
# 7) Notice the number of processors - I recommend updating OMP_NUM_THREADS in .afnirc
#		ETAC has heavy memory and processor needs


### update by Nathan Muncy on 11/30/18
#
# 8) Updated to support multiple pairwise comparisons
#		- if arrays A-C exist, will do an ETAC for AvB, AvC, BvC.
#		- will keep things straight with compList, so if you have the following:
#
#			compList=(first second third)
#			arrA=(1 2 3)
#			arrB=(11 22 33)
#			arrC=(111 222 333)
#			listX=ABC
#
#			then for the "first" comparison, you would get 1v11, 1v111, and 11v111
#			and for the "second" comparison, you would get 2v22, 2v222, and 22v222, etc.
#
# 9) Will now write out the grpAnalysis scripts to $outDir for your review.
#
# 10) added notes to each section
#
# 11) Added a section that will run MVMs via the ACF multiple comparison method.
#
#			compList=(first second third)
#			arrA=(1 2 3)
#			arrB=(11 22 33)
#			arrC=(111 222 333)
#			listX=ABC
#			bsArr=(Con Aut)
#
#			then for the "third" comparison, you would get 2(Con, Aut) x 3(3,33,333) comparison







### --- Set up --- ###										###??? update variables/arrays
#
# This is where the script will orient itself.
# Notes are supplied, and is the only section
# that really needs to be changed for each
# experiment.


# General variables
workDir=~/compute/STT_new									# par dir of data
outDir=${workDir}/Analyses/test								# where output will be written (should match step3)
tempDir=~/bin/Templates/vold2_mni							# desired template
priorDir=${tempDir}/priors_ACT								# location of atropos priors
refFile=${workDir}/s1295/TEST24_scale.nii.gz				# reference file for dimensions etc
mask=Intersection_GM_mask+tlrc								# this will be made


# grpAnalysis
doETAC=0													# Toggle ETAC analysis
doACF=1														# MVM
runIt=1														# whether ETAC/MVM actually run (and not just written)

thr=0.3														# thresh value for Group_EPI_mask, ref Group_EPI_mean

compList=(SpT1 T1 T1pT2 )									# matches decon prefixes, and will be prefix of output files
compLen=${#compList[@]}

arrA=(33 53 59)												# setA beh sub-brik for compList. Must be same length as compList
arrB=(36 56 62)												# setB
arrC=(39 59 65)
listX=ABC													# list of arr? used, for building permutations

namA=(RpH Hit FpH)											# names of behaviors from arrA. Must be same length as arrA
namB=(RpF FA FpF)
namC=(RpCR CR MpH)


# ETAC arrs
blurX=({2..4})    											# blur multipliers (integer)
pval_list=(0.01 0.005 0.001)								# p-value thresholds


# MVM vars/arrs
bsArr=(hold)												# Bx-subject variables (groups)
bsList=hold													# Needed when bsArr > 1. List where col1 = subj identifier, col2 = group membership (e.g. s1295 Con)
blurM=2														# blur multiplier, float/int




### --- Functions --- ###

# search array for string
MatchString () {

	local e match="$1"

	shift
	for e; do
		[[ "$e" == "$match" ]] && return 0
	done
	return 1
}


# make perumtations of length 2
MakePerm () {

	local items=$1
	local string i j hArr

	for ((i=0; i<${#items}; i++)); do

		string=${items:$i+1}
		for ((j=0; j<${#string}; j++)); do

			hArr+="${items:$i:1}${string:$j:1} "
		done
	done
	echo $hArr
}




### --- Checks, Permutations --- ###

# check
if [ ${#arrA[@]} != ${#arrB[@]} ] || [ ${#arrA[@]} != $compLen ]; then
	echo "Replace user and try again - grpAnalysis variables incorrect" >&2
	exit 1
fi

if [ ${#bsArr[@]} -gt 1 ] && [ ! -s $bsList ]; then
	echo "Replace user and try again - MVM vars/arrs incorrect" >&2
	exit 1
fi


# make permutation lists
arr=(`MakePerm $listX`)
alpha=(`echo {A..Z}`)
wsList=(${alpha[@]:0:${#listX}})


for ((a=0; a<${#bsArr[@]}; a++)); do
	tmpList+=$a
done
arrBS=(`MakePerm $tmpList`)




### --- Create Masks --- ###
#
# This section will create a group mean intersection mask
# then threshold it at $thr to create a binary intersection mask.
# A gray matter mask will be constructed, and then the GM mask
# will be multiplied with the intersection mask to create a
# single GM intersection mask


cd $outDir

if [ $runIt == 1 ]; then

	# ref file
	if [ ! -f $refFile ]; then
		3dcopy ${refFile%.nii*}+tlrc $refFile							###??? check this line
	fi


	# intersection mask
	if [ ! -f Group_epi_mask.nii.gz ]; then

		for i in ${workDir}/s*; do                                    	###??? check this line
			list+="${i}/mask_epi_anat+tlrc "
		done

		3dMean -prefix ${outDir}/Group_epi_mean.nii.gz $list
		3dmask_tool -input $list -frac $thr -prefix ${outDir}/Group_epi_mask.nii.gz
	fi


	# make $mask
	if [ ! -f ${mask}.HEAD ]; then

		# GM mask
		c3d ${priorDir}/Prior2.nii.gz ${priorDir}/Prior4.nii.gz -add -o Prior_GM.nii.gz
		3dresample -master $refFile -rmode NN -input Prior_GM.nii.gz -prefix Template_GM_mask.nii.gz

		# combine GM and intersection mask
		c3d Template_GM_mask.nii.gz Group_epi_mask.nii.gz -multiply -o Intersection_GM_prob_mask.nii.gz
		c3d Intersection_GM_prob_mask.nii.gz -thresh 0.1 1 1 0 -o Intersection_GM_mask.nii.gz
		3dcopy Intersection_GM_mask.nii.gz $mask
	fi


	# get template
	if [ ! -f vold2_mni_brain+tlrc.HEAD ]; then
		cp ${tempDir}/vold2_mni_brain+tlrc* .
	fi
fi



### --- ETAC --- ###
#
# This section will generate needed arguments for ETAC, including
# the p-values, blur sizes, and z-scores. It will then create a
# subject list which will exclude participants who were excessively
# noisy. An ETAC script will then be written, ran, and then the
# output will be combined across blurs and p-values to yield
# FINALall_something files. The extract sig clusters and stitching
# may have a better solution.
#
# Currently, this will only do paired comparisons e.g. BehA vs BehB,
# but multiple t-tests are supported (AvB, AvC, BvC)


if [ $doETAC == 1 ]; then

	# generate z,plist
	c=0; for i in ${pval_list[@]}; do

		zval_list[$c]=`ccalc "qginv(${pval_list[$c]})"`
		pval_hold+="${pval_list[$c]},"
		let c=$[$c+1]
	done

	pval_all=${pval_hold%?}


	# gen $blur
	gridSize=`fslhd $refFile | grep "pixdim1" | awk '{print $2}'`
	int=`printf "%.0f" $gridSize`

	c=0; for i in ${blurX[@]}; do
		hold="$(($int * $i))"
		blur+="$hold "
		blurArr[$c]=$hold
		let c=$[$c+1]
	done


	# do ETAC for e/permutation of sub-briks, for each $compList
	c=0; while [ $c -lt $compLen ]; do

		# variables													###??? check these, especially if running multiple tests on same decon file
		pref=${compList[$c]}
		scan=${pref}_stats+tlrc
		outPre=${pref}_ETAC_Decon

		for a in ${arr[@]}; do

			# unpack sub-brik value/name for e/permutation set
			eval declare -a var1=(arr${a:0:1})
			eval declare -a var2=(arr${a:1:1})
			eval declare -a nam1=(nam${a:0:1})
			eval declare -a nam2=(nam${a:1:1})

			brik1=$(eval echo \${${var1}[$c]})
			brik2=$(eval echo \${${var2}[$c]})
			name1=$(eval echo \${${nam1}[$c]})
			name2=$(eval echo \${${nam2}[$c]})

			out=${outPre}_${name1}-${name2}


			# construct setA and setB of subjects not found in info_rmSubj.txt
			arrRem=(`cat ${outDir}/info_rmSubj_${pref}.txt`)
			unset ListA
			unset ListB

		    for i in ${workDir}/s*; do                      		###??? check this line

				subj=${i##*\/}
				MatchString "$subj" "${arrRem[@]}"

				if [ $? == 1 ]; then
					ListA+="$subj ${i}/${scan}[${brik1}] "
					ListB+="$subj ${i}/${scan}[${brik2}] "
				fi
			done


			# write script
			echo "3dttest++ \\
				-paired \\
				-mask $mask \\
				-prefix $out \\
				-prefix_clustsim ${out}_clustsim \\
				-ETAC \\
				-ETAC_blur $blur \\
				-ETAC_opt name=NN1:NN1:2sid:pthr=$pval_all \\
				-setA A $ListA \\
				-setB B $ListB" > ${outDir}/${out}.sh


			# Do ETAC
			if [ $runIt == 1 ]; then

				if [ ! -f ${out}.B${blur:0:1}.0.nii ]; then
					source ${outDir}/${out}.sh
				fi

	            # Extract sig clusters
				for j in ${blurArr[@]}; do
			        if [ ! -f ${out}_B${j}_allMask+tlrc.HEAD ]; then

			            3dMultiThresh \
			            -mthresh ${out}_clustsim.NN1.ETAC.mthresh.B${j}.0.5perc.nii \
			            -input ${out}.B${j}.0.nii \
			            -1tindex 1 \
			            -prefix ${out}_B${j}_allMask \
			            -allmask ${out}_B${j}_binMask
			        fi
			    done

			    # Stitch together
			    if [ ! -f FINALall_${out}+tlrc.HEAD ]; then
					3dMean -sum -prefix FINALall_${out} ${out}_B*_allMask+tlrc.HEAD
			    fi
		    fi
		done

		let c=$[$c+1]
	done
fi




### --- ACF --- ###
#
# This will blur both stats and errts files according to $blurM, and
# then the blurred errts files will be used to model noise with
# an auto-correlation function. MVM scripts will be written and run,
# and none of this will happen on participants who move too much.
# A variable number of bx/wi subj variables is accepted, but this
# will not run a t-test.
#
# Currently, MVM post-hoc comparisons are permutations of bx/wi-subj
# variables. I.e. Behaviors A B C for groups Aut Con will yield
# comparisons of Aut-Con A-B, Aut-Con A-C, Aut-Con B-C. I could build
# more comparisons in the future.


if [ $doACF == 1 ]; then

	if [ ${#listX} -lt 3 ] && [ ${#bsArr[@]} == 1 ]; then
		echo "Replace user and try again - don't use ACF for a pairwise comparison" >&2
		exit 1
	fi


	arrCount=0; while [ $arrCount -lt $compLen ]; do

		pref=${compList[$arrCount]}
		outPre=${pref}_MVM_Decon
		print=ACF_raw_${pref}.txt


		# make subj list
		unset subjList

		for j in ${workDir}/s*; do                      			###??? check this line

			arrRem=(`cat info_rmSubj_${pref}.txt`)
			subj=${j##*\/}
			MatchString "$subj" "${arrRem[@]}"

			if [ $? == 1 ]; then
				subjList+=("$subj ")
			fi
		done


		# blur, determine parameter estimate
		gridSize=`fslhd $refFile | grep "pixdim1" | awk '{print $2}'`
		blurH=`echo $gridSize*$blurM | bc`
		blurInt=`printf "%.0f" $blurH`
		scan=${pref}_stats_blur${blurInt}+tlrc

		if [ $runIt == 1 ]; then
			if [ ! -s $print ]; then
				for k in ${subjList[@]}; do
					for m in stats errts; do

						hold=${workDir}/${k}/${pref}_${m}
						if [ ! -f ${hold}_blur${blurInt}+tlrc.HEAD ]; then
							3dmerge -prefix ${hold}_blur${blurInt} -1blur_fwhm $blurInt -doall ${hold}+tlrc
						fi
					done

					file=${workDir}/${k}/${pref}_errts_blur${blurInt}+tlrc
					3dFWHMx -mask $mask -input $file -acf >> $print
				done
			fi


			# simulate noise, determine thresholds
			if [ ! -s ACF_MC_${pref}.txt ]; then

				sed '/ 0  0  0    0/d' $print > tmp

				xA=`awk '{ total += $1 } END { print total/NR }' tmp`
				xB=`awk '{ total += $2 } END { print total/NR }' tmp`
				xC=`awk '{ total += $3 } END { print total/NR }' tmp`

				3dClustSim -mask $mask -LOTS -iter 10000 -acf $xA $xB $xC > ACF_MC_${pref}.txt
				rm tmp
			fi
		fi


		# set up - determine/construct variables for script
		unset conVar gltCount dataFrame

		if [ ${#bsArr[@]} -gt 1 ]; then


			# header, bx-subj title
			bsVars="'BSVARS'"
			header="Subj $bsVars WSVARS InputFile"


			# make $conVar (post-hoc comparisons)
			for x in ${!arrBS[@]}; do

				h1=${arrBS[$x]:0:1}
				h2=${arrBS[$x]:1:1}

				bsCon="1*${bsArr[$h1]} -1*${bsArr[$h2]}"
				bsLab=${bsArr[$h1]}-${bsArr[$h2]}

				for y in ${!arr[@]}; do

					gltCount=$[$gltCount+1]
					ws1h=${arr[$y]:0:1}
					ws2h=${arr[$y]:1:1}

					eval declare -a nam1=(nam${ws1h})
					eval declare -a nam2=(nam${ws2h})
					name1=$(eval echo \${${nam1}[$arrCount]})
					name2=$(eval echo \${${nam2}[$arrCount]})

					conVar+="-gltLabel $gltCount ${bsLab}_${name1}-${name2} -gltCode $gltCount '${bsVars}: $bsCon WSVARS: 1*$name1 -1*$name2' "

				done
			done


			# determine group membership, write dataframe
			bsSubj=(`cat $bsList | awk '{print $1}'`)
			bsGroup=(`cat $bsList | awk '{print $2}'`)

			for m in ${subjList[@]}; do
				for n in ${!bsSubj[@]}; do
					if [ $m == ${bsSubj[$n]} ]; then
						for o in ${wsList[@]}; do

							brik=$(eval echo \${arr${o}[$arrCount]})
							name=$(eval echo \${nam${o}[$arrCount]})

							dataFrame+="$m ${bsGroup[$n]} $name ${workDir}/${m}/'${scan}[${brik}]' "
						done
					fi
				done
			done

		else
			bsVars=1
			header="Subj WSVARS InputFile"

			for y in ${!arr[@]}; do

				gltCount=$[$gltCount+1]
				ws1h=${arr[$y]:0:1}
				ws2h=${arr[$y]:1:1}

				eval declare -a nam1=(nam${ws1h})
				eval declare -a nam2=(nam${ws2h})
				name1=$(eval echo \${${nam1}[$arrCount]})
				name2=$(eval echo \${${nam2}[$arrCount]})

				conVar+="-gltLabel $gltCount ${name1}-${name2} -gltCode $gltCount 'WSVARS: 1*$name1 -1*$name2' "
			done

			for m in ${subjList[@]}; do
				for n in ${wsList[@]}; do

					brik=$(eval echo \${arr${n}[$arrCount]})
					name=$(eval echo \${nam${n}[$arrCount]})

					dataFrame+="$m $name ${workDir}/${m}/'${scan}[${brik}]' "
				done
			done
		fi


		# write script
		echo "module load r/3/5

			3dMVM -prefix $outPre \\
			-jobs 10 \\
			-mask $mask \\
			-bsVars $bsVars \\
			-wsVars 'WSVARS' \\
			-num_glt $gltCount \\
			$conVar \\
			-dataTable \\
			$header \\
			$dataFrame" > ${outDir}/${outPre}.sh


		# run MVM
		if [ $runIt == 1 ]; then
			if [ ! -f ${outPre}+tlrc.HEAD ]; then
				source ${outDir}/${outPre}.sh
			fi
		fi

		let arrCount=$[$arrCount+1]
	done
fi
