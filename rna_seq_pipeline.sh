#!/bin/bash

# RNAseq Pipeline - A streamlined pipeline for RNA-seq data analysis
# This script automates RNA-seq data processing from raw FASTQ files to aligned BAMs
# Version: 2.0

set -e  # Exit on error

# Default parameters
FASTQ_DIR="."
GENOME="hg19"
LEADING=0
TRAILING=0
CROP=50
OUTPUT_DIR="results"
THREADS=$(nproc --ignore=1)  # Use all but one thread by default

# Display help information
show_help() {
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  -i, --input DIR       Directory containing FASTQ files (default: current directory)"
    echo "  -g, --genome STRING   Reference genome to use: hg19 or hg38 (default: hg19)"
    echo "  -o, --output DIR      Output directory (default: results)"
    echo "  -c, --crop INT        Crop reads to this length (default: 50)"
    echo "  -l, --leading INT     Remove leading low quality bases (default: 0)"
    echo "  -t, --trailing INT    Remove trailing low quality bases (default: 0)"
    echo "  -p, --threads INT     Number of threads to use (default: all but one)"
    echo "  -h, --help            Display this help message and exit"
    echo
    echo "Example: $0 -i fastq_data -g hg38 -o my_results -c 60 -p 8"
    echo
    exit 0
}

# Parse command line arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input)
                FASTQ_DIR="$2"
                shift 2
                ;;
            -g|--genome)
                GENOME="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -c|--crop)
                CROP="$2"
                shift 2
                ;;
            -l|--leading)
                LEADING="$2"
                shift 2
                ;;
            -t|--trailing)
                TRAILING="$2"
                shift 2
                ;;
            -p|--threads)
                THREADS="$2"
                shift 2
                ;;
            -h|--help)
                show_help
                ;;
            *)
                echo "Unknown option: $1"
                show_help
                ;;
        esac
    done

    # Remove trailing slash if present
    FASTQ_DIR=$(echo "$FASTQ_DIR" | sed 's/\/$//')
    
    # Validate input
    if [ ! -d "$FASTQ_DIR" ]; then
        echo "Error: Input directory $FASTQ_DIR does not exist"
        exit 1
    fi
    
    if [[ "$GENOME" != "hg19" && "$GENOME" != "hg38" ]]; then
        echo "Error: Genome must be either hg19 or hg38"
        exit 1
    fi
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    
    # Print configuration
    echo "=== RNA-seq Pipeline Configuration ==="
    echo "Input directory: $FASTQ_DIR"
    echo "Reference genome: $GENOME"
    echo "Output directory: $OUTPUT_DIR"
    echo "Crop length: $CROP"
    echo "Leading quality: $LEADING"
    echo "Trailing quality: $TRAILING" 
    echo "Threads: $THREADS"
    echo "==================================="
}

# Setup required directories
setup_directories() {
    echo "Setting up reference directories..."
    
    GENOMEDIR="$HOME/GENOMEDIR"
    
    # Directories to create
    DIRECTORIES=(
        "$GENOMEDIR"
        "$GENOMEDIR/hg19_index"
        "$GENOMEDIR/hg38_index"
        "$GENOMEDIR/picard"
        "$GENOMEDIR/adapter"
        "$GENOMEDIR/RefFlat"
    )

    # Check and create directories
    for DIR in "${DIRECTORIES[@]}"; do
        if [ ! -d "$DIR" ]; then
            mkdir -p "$DIR"
            echo "Created directory: $DIR"
        fi
    done
}

# Download necessary genome files
download_genome_files() {
    echo "Downloading genome files for $GENOME..."
    
    if [ "$GENOME" = "hg19" ]; then
        FILES=(
            "hg19.fa.gz:https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
            "hg19.refGene.gtf.gz:https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz"
        )
    else
        FILES=(
            "hg38.fa.gz:https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz"
            "hg38.refGene.gtf.gz:https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz"
        )
    fi
    
    for FILE_URL in "${FILES[@]}"; do
        FILE="${FILE_URL%%:*}"
        URL="${FILE_URL#*:}"
        
        if [ ! -e "$GENOMEDIR/$FILE" ]; then
            echo "Downloading $FILE..."
            wget -q --show-progress "$URL" -O "$GENOMEDIR/$FILE"
            echo "Downloaded $FILE"
        else
            echo "$FILE already exists, skipping download"
        fi
    done
}

# Download tools and adapter sequences
download_tools() {
    echo "Downloading necessary tools and sequences..."
    
    # Download Picard.jar
    if [ ! -e "$GENOMEDIR/picard/picard.jar" ]; then
        echo "Downloading picard.jar..."
        wget -q --show-progress -O $GENOMEDIR/picard/picard.jar https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar
        chmod +x $GENOMEDIR/picard/picard.jar
        echo "Downloaded picard.jar"
    fi
    
    # Download adapter sequences
    ADAPTERS=("TruSeq3-SE.fa" "TruSeq3-PE.fa")
    for ADAPTER in "${ADAPTERS[@]}"; do
        if [ ! -e "$GENOMEDIR/adapter/$ADAPTER" ]; then
            echo "Downloading $ADAPTER..."
            wget -q --show-progress -O "$GENOMEDIR/adapter/$ADAPTER" "https://github.com/usadellab/Trimmomatic/raw/main/adapters/$ADAPTER"
            echo "Downloaded $ADAPTER"
        fi
    done
    
    # Download gtfToGenePred
    if [ ! -f "$GENOMEDIR/gtfToGenePred" ]; then
        echo "Downloading gtfToGenePred..."
        wget -q --show-progress https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred -O "$GENOMEDIR/gtfToGenePred"
        chmod +x "$GENOMEDIR/gtfToGenePred"
        echo "Downloaded gtfToGenePred"
    fi
}

# Function to decompress files
decompress_file() {
    local file="$1"
    local decompressed_file="${file%.gz}"

    if [[ -f "$decompressed_file" ]]; then
        echo "Decompressed file $decompressed_file already exists"
    elif [[ -f "$file" && "$file" =~ \.gz$ ]]; then
        echo "Decompressing $file..."
        gunzip -kf "$file"
        echo "Decompression completed"
    fi
}

# Get RefFlat file
get_refflat_file() {
    local genome="$1"
    local refFlatFile="$GENOMEDIR/RefFlat/${genome}_ref_flat.txt"

    if [ ! -f "$refFlatFile" ]; then
        echo "Downloading RefFlat file for $genome..."
        wget -q --show-progress "https://hgdownload.cse.ucsc.edu/goldenPath/$genome/database/refFlat.txt.gz" -O "$refFlatFile.gz"
        decompress_file "$refFlatFile.gz"
        echo "Downloaded and decompressed RefFlat file"
    else
        echo "RefFlat file for $genome already exists"
    fi
    
    return "$refFlatFile"
}

# Perform genome indexing
perform_indexing() {
    local genome_dir="$1"
    local genome_fasta="$2"
    local sjdb_gtf="$3"
    
    # STAR parameters
    local sjdboverhang=59  # splice junction database overhang = max(ReadLength)-1
    local genomesaindexnbases=14  # min(14, log2(GenomeLength)/2 - 1)

    # Check if the index already exists
    if [ ! -e "$genome_dir/SAindex" ]; then
        echo "Indexing genome for $GENOME..."
        STAR --runMode genomeGenerate \
             --genomeDir "$genome_dir" \
             --genomeFastaFiles "$genome_fasta" \
             --sjdbGTFfile "$sjdb_gtf" \
             --sjdbOverhang $sjdboverhang \
             --genomeSAindexNbases $genomesaindexnbases \
             --runThreadN $THREADS
        echo "Genome indexing complete"
    else
        echo "Genome index already exists for $GENOME"
    fi
}

# Setup genome files and indexing
setup_genome() {
    # Set paths for genome files
    if [ "$GENOME" = "hg19" ]; then
        genome_file="$GENOMEDIR/hg19.fa.gz"
        gtf_file="$GENOMEDIR/hg19.refGene.gtf.gz"
        genome_dir="$GENOMEDIR/hg19_index"
        refFlat="$GENOMEDIR/RefFlat/hg19_ref_flat.txt"
    else
        genome_file="$GENOMEDIR/hg38.fa.gz"
        gtf_file="$GENOMEDIR/hg38.refGene.gtf.gz"
        genome_dir="$GENOMEDIR/hg38_index"
        refFlat="$GENOMEDIR/RefFlat/hg38_ref_flat.txt"
    fi
    
    # Decompress files
    decompress_file "$genome_file"
    decompress_file "$gtf_file"
    
    # Update paths to decompressed files
    genome_file="${genome_file%.gz}"
    gtf_file="${gtf_file%.gz}"
    
    # Perform indexing
    perform_indexing "$genome_dir" "$genome_file" "$gtf_file"
    
    # Get RefFlat file
    get_refflat_file "$GENOME"
    
    # Return values needed for processing
    GENOME_FASTA="$genome_file"
    GTF_FILE="$gtf_file"
    GENOME_DIR="$genome_dir"
    REFFLAT="$refFlat"
}

# Process a single sample
process_sample() {
    local r1_file="$1"
    local basename=$(basename "$r1_file")
    local sample_name="${basename%%_R1*.fastq*}"
    local sample_dir="$OUTPUT_DIR/$sample_name"
    local r2_file="${r1_file/_R1/_R2}"
    local is_paired=false
    
    # Create sample directory
    mkdir -p "$sample_dir"
    
    echo "==== Processing sample: $sample_name ===="
    
    # Check if paired-end or single-end
    if [ -e "$r2_file" ]; then
        is_paired=true
        echo "Detected paired-end sequencing"
    else
        echo "Detected single-end sequencing"
    fi
    
    # Run FastQC
    echo "Running FastQC..."
    fastqc "$r1_file" -o "$sample_dir"
    if [ "$is_paired" = true ]; then
        fastqc "$r2_file" -o "$sample_dir"
    fi
    
    # Trimming step
    echo "Performing quality trimming..."
    if [ "$is_paired" = true ]; then
        trimmomatic PE -threads $THREADS -phred33 "$r1_file" "$r2_file" \
            "$sample_dir/${sample_name}_trimmed_R1.fastq" "$sample_dir/${sample_name}_untrimmed_R1.fastq" \
            "$sample_dir/${sample_name}_trimmed_R2.fastq" "$sample_dir/${sample_name}_untrimmed_R2.fastq" \
            ILLUMINACLIP:"$GENOMEDIR/adapter/TruSeq3-PE.fa":2:30:10 LEADING:$LEADING TRAILING:$TRAILING CROP:$CROP SLIDINGWINDOW:4:15 MINLEN:15
        
        read_files=("$sample_dir/${sample_name}_trimmed_R1.fastq" "$sample_dir/${sample_name}_trimmed_R2.fastq")
        prefix="$sample_dir/${sample_name}_PE_"
    else
        trimmomatic SE -threads $THREADS -phred33 "$r1_file" "$sample_dir/${sample_name}_trimmed.fastq" \
            ILLUMINACLIP:"$GENOMEDIR/adapter/TruSeq3-SE.fa":2:30:10 LEADING:$LEADING TRAILING:$TRAILING \
            CROP:$CROP SLIDINGWINDOW:4:15 MINLEN:10
            
        read_files=("$sample_dir/${sample_name}_trimmed.fastq")
        prefix="$sample_dir/${sample_name}_SE_"
    fi
    
    # Alignment step
    echo "Performing STAR alignment..."
    STAR --runThreadN "$THREADS" \
         --genomeDir "$GENOME_DIR" \
         --readFilesIn "${read_files[@]}" \
         --outFileNamePrefix "$prefix" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
         --sjdbGTFfile "$GTF_FILE"
    
    # Set input BAM file path
    if [ "$is_paired" = true ]; then
        input_bam="${prefix}Aligned.sortedByCoord.out.bam"
        gene_count_file="${prefix}ReadsPerGene.out.tab"
    else
        input_bam="${prefix}Aligned.sortedByCoord.out.bam"
        gene_count_file="${prefix}ReadsPerGene.out.tab"
    fi
    
    # Run Picard metrics
    echo "Collecting RNA-seq metrics..."
    java -jar "$GENOMEDIR/picard/picard.jar" CollectRnaSeqMetrics \
        I="$input_bam" \
        O="$sample_dir/${sample_name}_RNA_Metrics.txt" \
        REF_FLAT="$REFFLAT" \
        STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
        MINIMUM_LENGTH=10

    echo "Collecting alignment metrics..."
    java -jar "$GENOMEDIR/picard/picard.jar" CollectAlignmentSummaryMetrics \
        R="$GENOME_FASTA" \
        I="$input_bam" \
        O="$sample_dir/${sample_name}_mapping_metrics.txt"

    echo "Marking duplicates..."
    java -jar "$GENOMEDIR/picard/picard.jar" MarkDuplicates \
        I="$input_bam" \
        O="$sample_dir/${sample_name}_marked_duplicates.bam" \
        M="$sample_dir/${sample_name}_duplicate_metrics.txt"
    
    # Run MultiQC
    echo "Running MultiQC..."
    multiqc -f "$sample_dir" -o "$sample_dir"
    
    echo "Sample $sample_name processing complete"
    echo "===============================\n"
}

# Main function
main() {
    # Parse command line arguments
    parse_args "$@"
    
    # Setup directories
    setup_directories
    
    # Download tools and adapter sequences
    download_tools
    
    # Download genome files
    download_genome_files
    
    # Setup genome and get paths
    setup_genome
    
    echo "\n===============================\n"
    echo "Starting RNA-seq analysis pipeline..."
    
    # Count FASTQ files
    fastq_count=$(find "$FASTQ_DIR" -name "*_R1*.fastq*" | wc -l)
    if [ "$fastq_count" -eq 0 ]; then
        echo "Error: No FASTQ files found in $FASTQ_DIR"
        exit 1
    fi
    echo "Found $fastq_count sample(s) to process"
    
    # Process each sample
    sample_num=1
    for r1_file in "$FASTQ_DIR"/*_R1*.fastq*; do
        echo "Processing sample $sample_num of $fastq_count"
        process_sample "$r1_file"
        sample_num=$((sample_num + 1))
    done
    
    echo "All samples processed successfully!"
    echo "Results are available in: $OUTPUT_DIR"
}

# Execute main function
main "$@"
