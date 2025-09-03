
###
BATCHDIR="P0625"
PROJECTDIR="/hpc/pmc_kuiper/SingleCellWGS/"
JOBDIR="${PROJECTDIR}/CODE/jobs/Genotype/${BATCHDIR}/"
OUTPUTDIR="${PROJECTDIR}/ANALYSIS/includebulktumor/${BATCHDIR}/GenoType/"
INPUTDIR="${PROJECTDIR}/ANALYSIS/includebulktumor/${BATCHDIR}/HaplotypeCaller/"
REFGENOME="/hpc/pmc_kuiper/References/Boxtel/homo_sapiens.GRCh38.GATK.illumina/Homo_sapiens_assembly38.fasta"


mkdir -p "${JOBDIR}"
mkdir -p "${OUTPUTDIR}"

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
echo -e "Generating job for sample: ${BATCHDIR}"

echo "#!/bin/bash
#SBATCH --job-name=${BATCHDIR}.Genotype.sh
#SBATCH --output=${JOBDIR}/${BATCHDIR}.Genotype.sh.out
#SBATCH --error=${JOBDIR}/${BATCHDIR}.Genotype.sh.err
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
module load bcftools/1.17

#Combine gvcfs
java -Djava.io.tmpdir=\$TMPDIR -Xmx30G -jar \$GATK4 \\
CombineGVCFs \\
 -O ${OUTPUTDIR}/${BATCHDIR}.g.vcf.gz \\" > $JOBDIR/${BATCHDIR}.Genotype.sh

for i in ${INPUTDIR}/*.g.vcf.gz
do
echo " -V $i \\" >> $JOBDIR/${BATCHDIR}.Genotype.sh
done

echo " -R ${REFGENOME}

# Genotype Caller 
java -Djava.io.tmpdir=\$TMPDIR -Xmx30G -jar \$GATK4 \\
GenotypeGVCFs \\
 -R ${REFGENOME} \\
 -V ${OUTPUTDIR}/${BATCHDIR}.g.vcf.gz \\
 -O ${OUTPUTDIR}/${BATCHDIR}.output.vcf.gz \\

#Get snvs
bcftools view -m2 -M2 -v snps ${OUTPUTDIR}/${BATCHDIR}.output.vcf.gz > ${OUTPUTDIR}/${BATCHDIR}.output.snv.vcf.gz

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

" >> $JOBDIR/${BATCHDIR}.Genotype.sh



