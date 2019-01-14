#!/bin/bash


# Written by Nathan Muncy on 11/28/18


###--- Notes, in no particular order
#
# 1) the script will split the cluster file into multiple masks, and pull betas from each participant.
#
# 2) assumes clusters from step4 output have been saved in Clust_$fileArr format
#		will use the comparisonString portion to keep different group analyses straight
#		comparisonString should match the decon prefix (step4 $etacList)
#
# 3) assumes that decon files exist locally (bring them back from the supercomputer)



# Variables
workDir=/Volumes/Yorick/STT_bids/derivatives						###??? Update this section
grpDir=${workDir}/Analyses/grpAnalysis
clustDir=${grpDir}/etac_clusters
outDir=${grpDir}/etac_betas
refFile=${workDir}/sub-1295/run-1_STUDY_scale+tlrc					# reference file for dimensions etc


fileArr=(SpT1 T1 T1pT2 T2fT1)										# decon files from which betas will be extracted - should match step4. Only include prefixes for files that have results
arrA=(33 53 59 97)													# sub-bricks corresponding to $fileArr
arrB=(36 56 62 103)
arrLen=${#arrA[@]}

blurX=2																# first mulitplier from step 4


# function - search array for string
MatchString (){
	local e match="$1"
	shift
	for e; do [[ "$e" == "$match" ]] && return 0; done
	return 1
}


# check
if [ ${#fileArr[@]} != $arrLen ] || [ ${#arrB[@]} != $arrLen ]; then
	echo "Replace user and try again - script set up incorreclty" >&2
	exit 1
fi


# determine blur
gridSize=`3dinfo -dk $refFile`
blurH=`echo $gridSize*$blurX | bc`
blurInt=`printf "%.0f" $blurH`


# organize
mkdir $outDir
mkdir $clustDir

mv ${grpDir}/Clust* $clustDir



### pull mean betas for e/cluster from e/comparison from e/subject
cd $clustDir

c=0; while [ $c -lt $arrLen ]; do

	hold=${fileArr[$c]}
	betas=${arrA[$c]},${arrB[$c]}


	# make subj list
	unset subjHold
	arrRem=(`cat ${grpDir}/info_rmSubj_${hold}.txt`)
	for i in ${workDir}/s*; do
		subj=${i##*\/}
		MatchString "$subj" "${arrRem[@]}"
		if [ $? == 1 ]; then
			subjHold+="$subj "
		fi
	done
	subjList=(${subjHold})


	# split clust masks
	if [ ! -f Clust_${hold}_c1+tlrc.HEAD ]; then

		3dcopy Clust_${hold}_mask+tlrc ${hold}.nii.gz
		num=`3dinfo Clust_${hold}_mask+tlrc | grep "At sub-brick #0 '#0' datum type is short" | sed 's/[^0-9]*//g' | sed 's/^...//'`

		for (( j=1; j<=$num; j++ )); do
			if [ ! -f Clust_${hold}_c${j}+tlrc.HEAD ]; then

				c3d ${hold}.nii.gz -thresh $j $j 1 0 -o ${hold}_${j}.nii.gz
				3dcopy ${hold}_${j}.nii.gz Clust_${hold}_c${j}+tlrc
			fi
		done
		rm *.nii.gz
	fi


	# pull betas
	for i in Clust_${hold}_c*+tlrc.HEAD; do

		tmp=${i##*_}
		cnum=${tmp%+*}
		print=${outDir}/Betas_${hold}_${cnum}.txt
		> $print

		for j in ${subjList[@]}; do
			subjDir=${workDir}/${j}

			# blur
			if [ ! -f ${subjDir}/${hold}_stats_blur${blurInt}+tlrc.HEAD ]; then
				3dmerge -prefix ${subjDir}/${hold}_stats_blur${blurInt} -1blur_fwhm $blurInt -doall ${subjDir}/${hold}_stats+tlrc
			fi

			file=${subjDir}/${hold}_stats_blur${blurInt}+tlrc
			stats=`3dROIstats -mask $i "${file}[${betas}]"`
			echo "$j $stats" >> $print
		done
	done
	let c=$[$c+1]
done



### organize prints into Master files
cd $outDir

# make master files
for i in Betas*_c1.txt; do

	tmp=${i#*_}
	string=${tmp%_*}
    > Master_${string}_betas.txt
done


# fill master files
for j in Betas*; do

    tmp1=${j#*_}
	tmp2=${j##*_}

	maskN=${tmp2%.*}
	string=${tmp1%_c*}
	print=Master_${string}_betas.txt

	echo "Mask $maskN" >> $print
	cat $j >> $print
	echo >> $print
done


# make list for R script
> Final_List.txt
for i in Master*; do
	tmp=${i#*_}
    echo ${tmp%.*} >> Final_List.txt
done


# clean
rm Betas*txt
