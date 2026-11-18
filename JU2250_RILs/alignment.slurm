#!/bin/bash
#SBATCH --job-name=JUalign_NILs
#SBATCH --output=JUalign_NILs_%j.out
#SBATCH --error=JUalign_NILs_%j.err
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

cd /home/yinchen/scratch/merged_JUvcf_breakpoints

# ======== CONFIGURATION ========
REF="20250625_c_elegans_WS283.genome.fa"

# Load modules (adjust for your cluster environment)
module load bwa
module load samtools

# ======== STEP 1: Index reference if not yet done ========
if [ ! -f "${REF}.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index $REF
fi

# ======== STEP 2: Align all paired-end reads ========
for R1 in *_R1_001.fastq.gz
do
    SAMPLE=$(basename $R1 _R1_001.fastq.gz)
    R2=${SAMPLE}_R2_001.fastq.gz
    OUTBAM=${SAMPLE}.sorted.bam

    if [ ! -f "$R2" ]; then
        echo "Warning: missing R2 for $SAMPLE, skipping."
        continue
    fi

    echo "Processing sample: $SAMPLE"

    # Align, convert to BAM, sort
    bwa mem -t 8 $REF $R1 $R2 | samtools sort -@ 4 -o $OUTBAM -

    # Index BAM
    samtools index $OUTBAM

    echo "Finished $SAMPLE"
done

echo "✅ All alignments completed."

sbatch variant_calling.slurm
