PROJECTDIR="/hpc/pmc_kuiper/SingleCellWGS/"
PATIENT="P0625"
BATCH="${PATIENT}_includingbultumor"
VCF_PATH="${PROJECTDIR}/ANALYSIS/includebulktumor/${PATIENT}/GenoType/"
VCF="${PATIENT}.output.snv.vcf.gz"
EXTENSION=".vcf.gz"
BLACKLIST="${PROJECTDIR}/ANALYSIS/blacklist/PTA_blacklist_05032025.bed"
JOBDIR="${PROJECTDIR}/CODE/jobs/BlacklistFilter/${BATCH}/"

NUM_CPU=1

mkdir -p $JOBDIR

echo -e "created blacklist filter job for ${BATCH}"
echo -e """#!/bin/bash
#SBATCH --job-name=FilterVCFonBlacklist-${BATCH}.sh
#SBATCH --output=$JOBDIR/FilterVCFonBlacklist-${BATCH}.sh.out
#SBATCH --error=$JOBDIR/FilterVCFonBlacklist-${BATCH}.sh.err
#SBATCH --partition=cpu
#SBATCH --time=00:59:00
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task ${NUM_CPU}
#SBATCH --gres=tmpspace:10G
#SBATCH --nodes=1
#SBATCH --open-mode=append

set -e # exit if any subcommand or pipeline returns a non-zero status
set -u # exit if any uninitialised variable is used


startTime=\$(date +%s)
echo \"startTime: \$startTime\"

# Load required modules
module load bedtools/2.31.1

#Filter vcf on blacklist
bedtools intersect -a ${VCF_PATH}/${VCF} -b ${BLACKLIST} -wa -header -v > ${VCF_PATH}/$(basename ${VCF} ${EXTENSION}).blacklistfiltered.vcf

echo \"Finished filtering vcf file with blacklist\"

#Retrieve and check return code
returnCode=\$?
echo \"Return code \${returnCode}\"

if [ \"\${returnCode}\" -eq \"0\" ]
then
        
        echo -e \"Return code is zero, process was succesfull\n\n\"
        
else
  
        echo -e \"\nNon zero return code not making files final. Existing temp files are kept for debugging purposes\n\n\"
        #Return non zero return code
        exit 1
        
fi


#Write runtime of process to log file
endTime=\$(date +%s)
echo \"endTime: \$endTime\"

#Source: http://stackoverflow.com/questions/12199631/convert-seconds-to-hours-minutes-seconds-in-bash

num=\$endTime-\$startTime
min=0
hour=0
day=0
if((num>59));then
    ((sec=num%60))
    ((num=num/60))
    if((num>59));then
        ((min=num%60))
        ((num=num/60))
        if((num>23));then
            ((hour=num%24))
            ((day=num/24))
        else
            ((hour=num))
        fi
    else
        ((min=num))
    fi
else
    ((sec=num))
fi
echo \"Running time: \${day} days \${hour} hours \${min} mins \${sec} secs\"

""" > ${JOBDIR}/FilterVCFonBlacklist-${BATCH}.sh
