#!/bin/bash
#SBATCH --job-name=JUjoint_variant_calling
#SBATCH --output=JUjoint_variant_calling_%j.out
#SBATCH --error=JUjoint_variant_calling_%j.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

cd /home/yinchen/scratch/merged_JUvcf_breakpoints
# ======== CONFIGURATION ========
REF="20250625_c_elegans_WS283.genome.fa"
OUTPREFIX="NILs_joint"
THREADS=8

# Load required modules
module load bcftools
module load samtools

# ======== STEP 1: List all BAM files ========
BAMLIST=$(ls *.sorted.bam | tr '\n' ' ')

echo "Samples included in variant calling:"
echo $BAMLIST

# ======== STEP 2: Generate pileup and call variants ========
bcftools mpileup -f $REF -Ou -a FORMAT/AD,FORMAT/DP,INFO/AD \
    --threads $THREADS $BAMLIST | \
bcftools call -mv -Oz -o ${OUTPREFIX}.raw.vcf.gz

# ======== STEP 3: Index the VCF ========
bcftools index ${OUTPREFIX}.raw.vcf.gz

# ======== STEP 4: Basic filtering (optional but recommended) ========
# This filters out low-quality or low-depth sites
bcftools filter -s LOWQUAL -e 'QUAL<30 || INFO/DP<10' \
    ${OUTPREFIX}.raw.vcf.gz -Oz -o ${OUTPREFIX}.filtered.vcf.gz


bcftools index ${OUTPREFIX}.filtered.vcf.gz

echo "✅ Joint variant calling complete!"
echo "Output VCF: ${OUTPREFIX}.filtered.vcf.gz"

