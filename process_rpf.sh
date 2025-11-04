#!/bin/bash

## scratch folder for intermediate files
TMP_FOLDER=/tmp
## input file
INPUT_FILE=RPF_sample1.fastq.gz
## output folder for result files
OUTPUT_FOLDER=output
## resource limit
NUM_CORES=5
MAX_MEMORY=20G

## default is mouse RRNA mouse_rrna.fasta (downloaded from  https://www.ncbi.nlm.nih.gov/nuccore/BK000964)
## if you want to analyze human data, use human_rrna.fasta (downloaded from https://www.ncbi.nlm.nih.gov/nuccore/U13369.1)
RRNA_LOCATION=mouse_rrna.fasta
RRNA_LIBRARY_PATH=rrna_library/mouse_rrna
## These are the mouse chromosomes, change to 1..22 for human 
CHROMOSOMES=({1..19} X Y MT)


## Ensembl reference genome fasta and annotation files -- change to human if necessary
## For mouse:
## Downloaded and gunzipped from https://ftp.ensembl.org/pub/release-97/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
## and https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
## For human:
## Downloaded and gunzipped from https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
## and https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
REF_GENOME_FILE=Mus_musculus.GRCm38.dna.primary_assembly.fa
REF_GENOME_ANNOT_FILE=Mus_musculus.GRCm38.102.gtf



STAR_OVERHANG=36
REF_GENOME_BUILD_LOCATION=${TMP_FOLDER}/ref_genome_GRCm38_sj${STAR_OVERHANG}

## downloaded from RSeQC and corrected only to change chromosome names to remove the 'chr': e.g. 'chr3' -> '3'
## awk 'BEGIN{OFS="\t"} {if($1=="chrM") $1="MT"; else gsub(/^chr/, "", $1); print}' mm10_GENCODE_VM25_basic.bed > mm10_GENCODE_VM25_basic_corr.bed
REF_GENOME_BED=mm10_GENCODE_VM25_basic_corr.bed



## PROCESSING SCRIPT START -- DON'T CHANGE

set -e

## extract filename
BASENAME=$(basename "$INPUT_FILE")

# Check if the file ends with .fastq, .fq, or .fq.gz
if [[ ! "$BASENAME" =~ \.fastq\.gz$|\.fq\.gz$ ]]; then
    echo "Error: input file must end with .fastq.gz or .fq.gz"
    exit 1
fi

# Extract sample name (remove .fastq or  .fq at the end)
SAMPLE_NAME=$(echo "$BASENAME" | sed -E 's/(\.fastq\.gz|\.fq\.gz)$//')

echo "Sample name: $SAMPLE_NAME"


## Step 1: copy UMI from sequence to header
INTERMEDIATE_FILE1="${SAMPLE_NAME}_copy_umi.fastq.gz"
echo "Copying UMI with fumi_tools..."
fumi_tools copy_umi --threads $NUM_CORES --umi-length 12 -i "$INPUT_FILE" -o "$TMP_FOLDER/$INTERMEDIATE_FILE1"

# Step 2: Trim adapters

INTERMEDIATE_FILE2="${SAMPLE_NAME}_trimmed.fastq.gz"
echo "Trimming adapters with cutadapt..."
cutadapt --trim-n -j $NUM_CORES \
--match-read-wildcards -u 16 \
-n 4 -a AGATCGGAAGAGCACACGTCTG -a AAAAAAAA \
-a GAACTCCAGTCAC -e 0.166666 --nextseq-trim 20 -m 18 \
-o "$TMP_FOLDER/$INTERMEDIATE_FILE2" "$TMP_FOLDER/$INTERMEDIATE_FILE1"


## Step 3: remove ribosomal RNA sequences
INTERMEDIATE_FILE3="${SAMPLE_NAME}_filtered.fastq"

echo "Removing ribosomal RNA sequences..."

## need to build with bowtie on first usage
if ! [ -f "${RRNA_LIBRARY_PATH}.1.bt2" ]; then
    echo "Bowtie2 index not found. Building index..."
    bowtie2-build -f "$RRNA_LOCATION" "$RRNA_LIBRARY_PATH"
fi

bowtie2 \
	--seedlen=23 \
	-x "$RRNA_LIBRARY_PATH" \
	"$TMP_FOLDER/$INTERMEDIATE_FILE2" --un "$TMP_FOLDER/$INTERMEDIATE_FILE3" >/dev/null 

gzip $TMP_FOLDER/$INTERMEDIATE_FILE3

## Step 4: check read lengths -- for later analysis
READ_LENGTH_FILE="${SAMPLE_NAME}_read_length.txt"
echo "Saving read lengths to ${READ_LENGTH_FILE}..."
zcat "$TMP_FOLDER/$INTERMEDIATE_FILE3.gz" | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > "$OUTPUT_FOLDER/$READ_LENGTH_FILE"


## Step 5: Align to genome using STAR
## https://github.com/alexdobin/STAR

echo "Aligning with STAR..."

## Build index on first run
if ! [ -d "${REF_GENOME_BUILD_LOCATION}/" ]; then
    echo "Star index not found. Building index."
    mkdir  "${REF_GENOME_BUILD_LOCATION}/"
    # parameters recommendation: https://github.com/alexdobin/STAR/issues/255
    STAR  \
    	--runMode genomeGenerate --runThreadN $NUM_CORES --genomeDir $REF_GENOME_BUILD_LOCATION \
    	--genomeFastaFiles $REF_GENOME_FILE \
    	--sjdbGTFfile $REF_GENOME_ANNOT_FILE --sjdbOverhang $STAR_OVERHANG 
fi


# parameters from RiboLite paper and from https://github.com/alexdobin/STAR/issues/255 -- optimized to short reads
STAR \
       	--runThreadN $NUM_CORES --readFilesCommand zcat --genomeDir $REF_GENOME_BUILD_LOCATION \
       	--sjdbGTFfile $REF_GENOME_ANNOT_FILE  \
	--outFilterMismatchNmax 2 --outFilterMultimapNmax 20 --outFilterMatchNmin 16 --alignEndsType EndToEnd \
	--sjdbOverhang $STAR_OVERHANG --readFilesIn "$TMP_FOLDER/$INTERMEDIATE_FILE3.gz" \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM GeneCounts \
	--outFileNamePrefix "$TMP_FOLDER/$SAMPLE_NAME"
	

## Step 6: filter and deduplicate. 

echo "Filtering and deduplicating..."

SAMFILE_GENOME="$TMP_FOLDER/${SAMPLE_NAME}Aligned.sortedByCoord.out.bam"
SAMFILE_GENOME_FILTERED="$TMP_FOLDER/${SAMPLE_NAME}Aligned.sortedByCoord_filtered.out.bam"
SAMFILE_GENOME_FILTERED_DEDUP="$TMP_FOLDER/${SAMPLE_NAME}Aligned.sortedByCoord_filtered_dedup.out.bam"

samtools index "$SAMFILE_GENOME"
samtools view -b "$SAMFILE_GENOME" ${CHROMOSOMES[@]} > "$SAMFILE_GENOME_FILTERED"

fumi_tools dedup --threads $NUM_CORES --memory $MAX_MEMORY \
-i  "$SAMFILE_GENOME_FILTERED" \
-o "$SAMFILE_GENOME_FILTERED_DEDUP"



SAMFILE_TRANSCRIPTOME="$TMP_FOLDER/${SAMPLE_NAME}Aligned.toTranscriptome.out.bam"
SAMFILE_TRANSCRIPTOME_SORTED="$TMP_FOLDER/${SAMPLE_NAME}Aligned.toTranscriptome_sorted.out.bam"
SAMFILE_TRANSCRIPTOME_SORTED_DEDUP="$TMP_FOLDER/${SAMPLE_NAME}Aligned.toTranscriptome_sorted_dedup.out.bam"

samtools sort -@ 10 -o "$SAMFILE_TRANSCRIPTOME_SORTED" $SAMFILE_TRANSCRIPTOME

fumi_tools dedup --threads $NUM_CORES --memory $MAX_MEMORY \
-i  "$SAMFILE_TRANSCRIPTOME_SORTED" \
-o "$SAMFILE_TRANSCRIPTOME_SORTED_DEDUP"




## Step 7: Write out read distribution with RSeQC
echo "Exporting read distribution..."

if ! command -v read_distribution.py &> /dev/null; then
    echo "rseqc not found. Installing via pip..."
    pip install --user RSeQC || { echo "pip install failed"; exit 1; }
else
    echo "rseqc is already installed."
fi



read_distribution.py -i "$SAMFILE_GENOME_FILTERED_DEDUP" -r ${REF_GENOME_BED}> "$OUTPUT_FOLDER/${SAMPLE_NAME}_read_distribution.txt"


## Step 8: metaplots
echo "Generating metaplots"
REF_GENOME_ANNOT_FILE_UPDATED=${REF_GENOME_ANNOT_FILE/.gtf/_updated.gtf}

if ! command -v GTFupdate &> /dev/null; then
    echo "RiboCode not found. Installing via pip..."
    pip install --user RiboCode || { echo "pip install failed"; exit 1; }
else
    echo "RiboCode is already installed."
fi


if ! [ -f "${REF_GENOME_ANNOT_FILE_UPDATED}" ]; then
  echo "GTF file not yet updated. Updating..."
  GTFupdate "${REF_GENOME_ANNOT_FILE}" > "${REF_GENOME_ANNOT_FILE_UPDATED}"
fi

if ! [ -d "${TMP_FOLDER}/ribocode_annot" ]; then 
  echo "Ribocode annot not yet built. Building..."
  prepare_transcripts -g "${REF_GENOME_ANNOT_FILE_UPDATED}" -f "${REF_GENOME_FILE}" -o "$TMP_FOLDER/ribocode_annot"

fi

## ribo metaplots
ls -1 "${SAMFILE_TRANSCRIPTOME_SORTED_DEDUP}"  > "$TMP_FOLDER/file_${SAMPLE_NAME}.txt"
#echo metaplots -a "$TMP_FOLDER/ribocode_annot" --minimum-length 18 -i "$TMP_FOLDER/file_${SAMPLE_NAME}.txt"
metaplots -o "${OUTPUT_FOLDER}/metaplots_${SAMPLE_NAME}" -a "$TMP_FOLDER/ribocode_annot" --minimum-length 18 -i "$TMP_FOLDER/file_${SAMPLE_NAME}.txt"


#Final Step: counting features
echo "Run featureCounts..."
FEATURE_COUNTS_OUTFILE="${SAMPLE_NAME}_feature_counts.txt"

featureCounts \
	-O `#Assign reads to all their overlapping meta-features.` \
	-s 1 -M \
	-t exon  `#Specify feature type in the GTF/GFF annotation to summarise the counts.` \
	-g gene_id  `#Specify attribute type in GTF/GFF annotation. This GTF/GFF determines the name of the features.` \
	-T $NUM_CORES `#number of cores` \
	-a "$REF_GENOME_ANNOT_FILE" `#annotation` \
	-o "${OUTPUT_FOLDER}/${FEATURE_COUNTS_OUTFILE}" `#output` \
	"$SAMFILE_GENOME_FILTERED_DEDUP" `#input`


echo "Done. Intermediate files were left in the temporary folder for your convenience. Don't forget to clean up!"



