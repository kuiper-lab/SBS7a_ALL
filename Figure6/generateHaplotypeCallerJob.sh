
###
BATCHDIR="P0625"
PROJECTDIR="/hpc/pmc_kuiper/SingleCellWGS/"
JOBDIR="${PROJECTDIR}/CODE/jobs/HaplotypeCaller/${BATCHDIR}/"
OUTPUTDIRHAPLO="${PROJECTDIR}/ANALYSIS/includebulktumor/${BATCHDIR}/HaplotypeCaller/"
SAMPLEFILE="${PROJECTDIR}/CODE/sampleLists/sampleListHaploGeno.${BATCHDIR}.txt"
VCF_FILE="${PROJECTDIR}/DATA/Treebuilding/${BATCHDIR}/MQandadjustedVAF025filtered/${BATCHDIR}_filtered_snv.vcf"
CRAM_EXT=".bam"
REFGENOME="/hpc/pmc_kuiper/References/Boxtel/homo_sapiens.GRCh38.GATK.illumina/Homo_sapiens_assembly38.fasta"


mkdir -p "${JOBDIR}"
mkdir -p "${OUTPUTDIRHAPLO}"

# Job settings
DATE=$(date --iso-8601=seconds)
NUMCPUS=16
WALLTIME="47:59:00"
MEM="30G"
JVMHEAP="30g"
TMPSPACE="700G"
NUM_CORES=8
NUM_FORKS=8

# Read samplefile list
while IFS='\t' read -r line || [[ $line ]]
do
    arr=($line)
    NORMAL="${arr[0]}"

    NORMAL_CRAM_FILE="${arr[1]}"
	
    echo -e "Generating job for sample: ${NORMAL}"

    echo "#!/bin/bash
#SBATCH --job-name=${NORMAL}.HaplotypeCaller.sh
#SBATCH --output=${JOBDIR}/${NORMAL}.HaplotypeCaller.sh.out
#SBATCH --error=${JOBDIR}/${NORMAL}.HaplotypeCaller.sh.err
#SBATCH --partition=cpu
#SBATCH --time=36:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=${NUM_CORES}
#SBATCH --gres=tmpspace:300G 
#SBATCH --nodes=1
#SBATCH --open-mode=append

set -e # exit if any subcommand or pipeline returns a non-zero status
set -u # exit if any uninitialised variable is used

startTime=\$(date +%s)
echo \"startTime: \$startTime\"


# load the required modules
module load gatk/4.2.0.0

# Haplotype Caller
java -Djava.io.tmpdir=\$TMPDIR -Xmx30G -jar \$GATK4 \\
HaplotypeCaller \\
 -R ${REFGENOME} \\
 -O ${OUTPUTDIRHAPLO}/${NORMAL}.g.vcf.gz \\
 -I ${NORMAL_CRAM_FILE} \\
 -L ${VCF_FILE} \\
 -ERC GVCF
# -L ${BEDFILE} \\

# Retrieve and check return code
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


# Write runtime of process to log file
endTime=\$(date +%s)
echo \"endTime: \$endTime\"


# Source: http://stackoverflow.com/questions/12199631/convert-seconds-to-hours-minutes-seconds-in-bash

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

" > $JOBDIR/${NORMAL}.HaplotypeCaller.sh

done < ${SAMPLEFILE}



