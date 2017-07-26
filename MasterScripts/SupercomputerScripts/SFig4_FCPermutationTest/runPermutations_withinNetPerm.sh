#!/bin/bash

##NM3 supercomputer (at Rutgers University-Newark) batch script
##
##Author: Taku Ito
##
## This script runs a MATLAB command on a large number of subjects using the supercomputer queuing system (Slurm).
##


##*** Modify these variables:
scriptDir="/work/ti61/projects/ModalityControl/docs/scripts_master/NM3Scripts/FC_PermutationTest/"
subjBatchScriptDir="/work/ti61/projects/ModalityControl/docs/scripts_master/NM3Scripts/FC_PermutationTest/subjbatch/"
jobNamePrefix="cycle_"
# Completed: "013 032 033"  
##Make and execute a batch script for each subject

for cycle in {0..10}
do

 	cd ${subjBatchScriptDir}
 	
	batchFilename=${cycle}_pythonBatch.sh
	
	echo "#!/bin/bash" > $batchFilename
	echo "#SBATCH --time=12:00:00" >> $batchFilename
	echo "#SBATCH --nodes=1" >> $batchFilename
	echo "#SBATCH --ntasks=1" >> $batchFilename
	echo "#SBATCH --partition=day-long,week-long,month-long" >> $batchFilename
	echo "#SBATCH --job-name=${jobNamePrefix}${cycle}" >> $batchFilename
	echo "#SBATCH --output=slurm.${jobNamePrefix}${cycle}.out" >> $batchFilename
	echo "#SBATCH --error=slurm.${jobNamePrefix}${cycle}.err" >> $batchFilename
	echo "#SBATCH --cpus-per-task=20" >> $batchFilename
	
	echo "#Run the python command" >> $batchFilename
	echo "cd $scriptDir" >> $batchFilename
        #echo "python ${scriptDir}/computeTaskRuleFC_GlasserParcels.py '${subjNum}'" >> $batchFilename
        echo "python -W ignore -c 'import runPermutations_withinNetPerm as run; run.runPermutations("${cycle}")'" >> $batchFilename

	#Submit the job
	sbatch $batchFilename
	
done
