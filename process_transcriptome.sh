#!/bin/bash

## USER INPUT -- CHANGE INPUT PARAMETERS HERE

## scratch folder for intermediate files
TMP_FOLDER=/tmp
## input files
INPUT_FILE_R1=sample1_R1.fastq.gz
INPUT_FILE_R2=sample1_R2.fastq.gz
## output folder for result files
OUTPUT_FOLDER=output
## resource limit
NUM_CORES=5
MAX_MEMORY=20G



## Ensembl reference genome fasta and annotation files -- change to human if necessary
## For mouse:
## Downloaded and gunzipped from https://ftp.ensembl.org/pub/release-97/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
## and https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
## For human:
## Downloaded and gunzipped from https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
## and https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
REF_GENOME_FILE=Mus_musculus.GRCm38.dna.primary_assembly.fa
REF_GENOME_ANNOT_FILE=Mus_musculus.GRCm38.102.gtf


STAR_OVERHANG=99
REF_GENOME_BUILD_LOCATION=${TMP_FOLDER}/ref_genome_GRCm38_sj${STAR_OVERHANG}


## PROCESSING SCRIPT START -- DON'T CHANGE

set -e

## extract filename
BASENAME_R1=$(basename "${INPUT_FILE_R1}")
BASENAME_R2=$(basename "${INPUT_FILE_R2}")
echo $BASENAME_R1
# Check if the file ends with .fastq.gz or .fq.gz
if [[ ! "$BASENAME_R1" =~ \.fastq\.gz$|\.fq\.gz$ ]]; then
    echo "Error: input file must end with .fastq.gz or .fq.gz"
    exit 1
fi

# Extract sample name (remove .fastq or  .fq at the end)
SAMPLE_NAME_R1=$(echo "${BASENAME_R1}" | sed -E 's/(\.fastq\.gz|\.fq\/gz)$//')
SAMPLE_NAME_R2=$(echo "${BASENAME_R2}" | sed -E 's/(\.fastq\.gz|\.fq\.gz)$//')
SAMPLE_NAME="${SAMPLE_NAME_R1//R1/SAMPLE}"
echo "Sample names (R1, R2): ${SAMPLE_NAME_R1}, ${SAMPLE_NAME_R2}"


INTERMEDIATE_FILE1="${SAMPLE_NAME_R1}_trimmed.fastq"
INTERMEDIATE_FILE2="${SAMPLE_NAME_R2}_trimmed.fastq"

cutadapt -j 4 \
	--trim-n \
	--match-read-wildcards \
  -a CTGTCTCTTATACACATCT -a AAAAAAAAAAA \
	-A CTGTCTCTTATACACATCT -A AAAAAAAAAAA \
  --nextseq-trim 20 -m 10 \
  -o "$TMP_FOLDER/$INTERMEDIATE_FILE1" -p "$TMP_FOLDER/$INTERMEDIATE_FILE2" \
  "$INPUT_FILE_R1" "$INPUT_FILE_R2"
        
gzip "$TMP_FOLDER/$INTERMEDIATE_FILE1"
gzip "$TMP_FOLDER/$INTERMEDIATE_FILE2"


## Step 5: Align to genome using STAR
## https://github.com/alexdobin/STAR

echo "Aligning with STAR..."

## Build index on first run
if ! [ -d "${REF_GENOME_BUILD_LOCATION}/" ]; then
    echo "Star index not found. Building index."
    mkdir "${REF_GENOME_BUILD_LOCATION}"
    # parameters recommendation: https://github.com/alexdobin/STAR/issues/255
    STAR  \
    	--runMode genomeGenerate --runThreadN $NUM_CORES --genomeDir $REF_GENOME_BUILD_LOCATION \
    	--genomeFastaFiles $REF_GENOME_FILE \
    	--sjdbGTFfile $REF_GENOME_ANNOT_FILE --sjdbOverhang $STAR_OVERHANG 
fi

STAR \
       	--runThreadN $NUM_CORES --readFilesCommand zcat --genomeDir $REF_GENOME_BUILD_LOCATION \
       	--sjdbGTFfile $REF_GENOME_ANNOT_FILE  \
	--sjdbOverhang $STAR_OVERHANG --readFilesIn "$TMP_FOLDER/$INTERMEDIATE_FILE1.gz" "$TMP_FOLDER/$INTERMEDIATE_FILE2.gz" \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM GeneCounts \
	--outFileNamePrefix "$TMP_FOLDER/$SAMPLE_NAME"
	

SAMFILE_GENOME="$TMP_FOLDER/${SAMPLE_NAME}Aligned.sortedByCoord.out.bam"
FEATURE_COUNTS_OUTFILE="${SAMPLE_NAME}_feature_counts.txt"



featureCounts \
	-O `#Assign reads to all their overlapping meta-features.` \
	-t exon  `#Specify feature type in the GTF/GFF annotation to summarise the counts.` \
	-g gene_id  `#Specify attribute type in GTF/GFF annotation. This GTF/GFF determines the name of the features.` \
	-T $NUM_CORES `#number of cores` \
	-p  --countReadPairs \
	-a "$REF_GENOME_ANNOT_FILE" `#annotation` \
	-o "${OUTPUT_FOLDER}/${FEATURE_COUNTS_OUTFILE}" `#output` \
	"$SAMFILE_GENOME" `#input`
	
	
	



