#!/bin/bash

##NM3 supercomputer (at Rutgers University-Newark) batch script
##
##Author: Taku Ito
##
## This script runs a MATLAB command on a large number of subjects using the supercomputer queuing system (Slurm).
##


##*** Modify these variables:
scriptDir="/work/ti61/projects/ModalityControl/docs/scripts_master/NM3Scripts/"
subjBatchScriptDir="/work/ti61/projects/ModalityControl/docs/scripts_master/NM3Scripts/subjbatch/"
jobNamePrefix="ite_"
subjNums="013 032 033 037 038 039 045 014 016 017 018 021 023 024 025 026 027 031 035 046 042 028 048 053 040 049 057 062 050 030 047 034"
# Completed: ""  
# Total: "013 032 033 037 038 039 045 014 016 017 018 021 023 024 025 026 027 031 035 046 042 028 048 053 040 049 057 062 050 030 047 034" 
##Make and execute a batch script for each subject

for subjNum in ${subjNums}
do

 	cd ${subjBatchScriptDir}
 	
	batchFilename=${subjNum}_pythonBatch.sh
	
	echo "#!/bin/bash" > $batchFilename
	echo "#SBATCH --nodes=1" >> $batchFilename
	echo "#SBATCH --ntasks=1" >> $batchFilename
	echo "#SBATCH --partition=week-long" >> $batchFilename
	echo "#SBATCH --time=48:00:00" >> $batchFilename
	echo "#SBATCH --job-name=${jobNamePrefix}${subjNum}" >> $batchFilename
	echo "#SBATCH --output=slurm.${jobNamePrefix}${subjNum}.out" >> $batchFilename
	echo "#SBATCH --error=slurm.${jobNamePrefix}${subjNum}.err" >> $batchFilename
	echo "#SBATCH --cpus-per-task=10" >> $batchFilename
	
	echo "#Run the python command" >> $batchFilename
	echo "cd $scriptDir" >> $batchFilename
        echo "/home/apps/MATLAB/R2016a/bin/matlab -nodisplay -r \"ActFlow_ITE_DecodePerformance_LogRegression_v2('${subjNum}'), exit\"" >> $batchFilename

	#Submit the job
	sbatch $batchFilename
	
done
