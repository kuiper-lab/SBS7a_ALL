#Set variables
BATCH="UVRisingVariants_24032023"
REFERENCE_GENOME="/Users/liannesuurenbroek/Projects/hypermutated_ALL/REFERENCES/Homo_sapiens_assembly38.fasta"
OUTPUT_FOLDER="/Users/liannesuurenbroek/Projects/hypermutated_ALL/ANALYSES/deepSequencing/bam-readcount/readcountFiles/$BATCH"
INPUT_DIR="/Users/liannesuurenbroek/Projects/hypermutated_ALL/DATA/iSeq/bam_with_duplicates/"
TARGET_FILE="/Users/liannesuurenbroek/Projects/hypermutated_ALL/ANALYSES/ddPCR/VariantIdentification/variant_positions_UVRisingVariants_bamreadcount.bed"

mkdir -p $OUTPUT_FOLDER
conda activate bamreadcount

#Run bam-readcount
for BAM in $(ls ${INPUT_DIR}/*.bam)
do
    SAMPLE_ID=$(basename ${BAM} _WXS.bam)
    BAM_PATH="${BAM}"
    echo $SAMPLE_ID
    echo $BAM_PATH

    # run commando
    bam-readcount \
    -f $REFERENCE_GENOME \
    -l $TARGET_FILE \
    -w 10 \
    $BAM_PATH >> $OUTPUT_FOLDER/$SAMPLE_ID.bam-readcount.variants.txt
done
