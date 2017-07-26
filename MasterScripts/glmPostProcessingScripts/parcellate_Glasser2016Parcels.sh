#!/bin/bash
#Script to create a parcellated CIFTI file, based on the Glasser2016 parcellation

basedir=/projects/IndivRITL/
datadir=${basedir}/data/
subjList="013 014 016 017 018 021 023 024 025 026 027 028 030 031 032 033 034 035 037 038 039 040 041 042 043 045 046 047 048 049 050 052 053 055 056 057 058 062"
#063 064 066 067 068 069 070 072 074 075 076 077 078 079 081 082 085 086 087 088 090 092 093 094 095 097 098 099 101 102 103 104 105 106 108 109 110 111 112 114 115 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 134 135 136 137 138 139 140 141"
runNameList="Rest1 Task1 Task2 Task3 Task4 Task5 Task6 Task7 Task8"
parcelFile=/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210.LR.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii

for subjNum in $subjList; do
  subjDir=${datadir}/${subjNum}/
  subjDirOutput=/projects2/ModalityControl2/data/${subjNum}
  if [ "`ls -d $subjDirOutput`" == "" ]; then mkdir $subjDirOutput; fi	#Making directory if it doesn't exist yet
  subjAnalysisDir=${subjDirOutput}/analysis/
  if [ "`ls -d $subjAnalysisDir`" == "" ]; then mkdir $subjAnalysisDir; fi	#Making directory if it doesn't exist yet

  for runName in $runNameList; do
    subjRunDir=${subjDir}/MNINonLinear/Results/${runName}/
    pushd ${subjRunDir}

    inputFileName=${runName}_Atlas.dtseries.nii

    #Run commands to parcellate data file
    # Note that the hemispheres for this file are flippd; however, since all subsequent analyses use this template, files will be flipped when converting back to cifti
    wb_command -cifti-parcellate ${inputFileName} ${parcelFile} COLUMN ${subjAnalysisDir}/${runName}_Atlas.LR.Glasser2016Parcels.32k_fs_LR.ptseries.nii -method MEAN

    popd
  done

done
