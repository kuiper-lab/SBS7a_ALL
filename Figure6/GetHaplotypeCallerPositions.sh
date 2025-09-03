PROJECTDIR="/hpc/pmc_kuiper/SingleCellWGS/"
BATCH="P0625"
UNFILTERED_VCF="${PROJECTDIR}/ANALYSIS/PTATO/Output_${BATCH}/ptato_vcfs/${BATCH}/${BATCH}.ptato.merged.vcf.gz"
FILTERED_VCFS="${PROJECTDIR}/ANALYSIS/PTATO/${BATCH}/*.vcf"
OUTPUTDIR="$PROJECTDIR/DATA/Treebuilding/${BATCH}/"
JOBDIR="$PROJECTDIR/CODE/jobs/CreateCellPhyInputVCF/"

NUM_CPU=1

mkdir -p $OUTPUTDIR
mkdir -p $JOBDIR

echo -e "created CreateCellPhyInputVCF job for ${BATCH}"
echo -e """#!/bin/bash
#SBATCH --job-name=CreateCellPhyInputVCF-${BATCH}.sh
#SBATCH --output=$JOBDIR/CreateCellPhyInputVCF-${BATCH}.sh.out
#SBATCH --error=$JOBDIR/CreateCellPhyInputVCF-${BATCH}.sh.err
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
module load bcftools/1.17

# Select all mutations after PTATO filtering
bedtools intersect -a ${UNFILTERED_VCF} -b ${FILTERED_VCFS} -wa -header -u > ${OUTPUTDIR}/${BATCH}_filtered.vcf
# Select only SNVs
bcftools view -m2 -M2 -v snps ${OUTPUTDIR}/${BATCH}_filtered.vcf > ${OUTPUTDIR}/${BATCH}_filtered_snv.vcf

echo \"Finished Generating Input VCF for CellPhy\"

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

""" > ${JOBDIR}/CreateCellPhyInputVCF-${BATCH}.sh
