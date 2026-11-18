#!/bin/bash
#SBATCH --job-name=JU2522_RILs_align_call
#SBATCH --output=JU2522_RILs_align_call_%j.out
#SBATCH --error=JU2522_RILs_align_call_%j.err
#SBATCH --time=25:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

cd /home/yinchen/scratch/JU2522_RILs

# ======== CONFIGURATION ========
REF="20250625_c_elegans_WS283.genome.fa"
OUTPREFIX="JU2522_RILs_joint"
THREADS=8

# Load modules (adjust for your cluster environment)
module load bwa
module load samtools
module load bcftools

# Stop the job if any command (or any part of a pipe) fails,
# so variant calling never runs on missing or truncated BAMs
set -eo pipefail

# ================================================================
#                        PART 1: ALIGNMENT
# ================================================================

# ======== STEP 1: Index reference if not yet done ========
if [ ! -f "${REF}.bwt" ]; then
    echo "Indexing reference genome (bwa)..."
    bwa index $REF
fi

if [ ! -f "${REF}.fai" ]; then
    echo "Indexing reference genome (samtools faidx)..."
    samtools faidx $REF
fi

# ======== STEP 2: Align each lane, then merge lanes per sample ========
# Files look like: JU2522_RIL_J11_S6_L001_R1_001.fastq.gz
#   PREFIX = JU2522_RIL_J11_S6   (everything before the lane)
#   SAMPLE = JU2522_RIL_J11      (PREFIX minus the _S# sample-sheet number)
#   LANE   = L001

PREFIXES=$(ls *_L[0-9][0-9][0-9]_R1_001.fastq.gz \
    | sed -E 's/_L[0-9]{3}_R1_001\.fastq\.gz$//' | sort -u)

BAMS=()

for PREFIX in $PREFIXES
do
    SAMPLE=$(echo $PREFIX | sed -E 's/_S[0-9]+$//')
    echo "Processing sample: $SAMPLE"

    LANE_BAMS=()

    for R1 in ${PREFIX}_L[0-9][0-9][0-9]_R1_001.fastq.gz
    do
        LANE=$(echo $R1 | sed -E 's/.*_(L[0-9]{3})_R1_001\.fastq\.gz$/\1/')
        R2=${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}
        LANE_BAM=${SAMPLE}_${LANE}.sorted.lane.bam

        if [ ! -f "$R2" ]; then
            echo "Warning: missing R2 for $SAMPLE $LANE, skipping this lane."
            continue
        fi

        echo "  Aligning $LANE"

        # Read group: same SM for all lanes so they count as one sample,
        # unique ID per lane
        RG="@RG\tID:${SAMPLE}.${LANE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA\tPU:${LANE}"

        bwa mem -t $THREADS -R "$RG" $REF $R1 $R2 \
            | samtools sort -@ 4 -o $LANE_BAM -

        LANE_BAMS+=("$LANE_BAM")
    done

    if [ ${#LANE_BAMS[@]} -eq 0 ]; then
        echo "Warning: no complete lanes for $SAMPLE, skipping sample."
        continue
    fi

    # Merge lanes into one BAM per sample
    OUTBAM=${SAMPLE}.sorted.bam
    echo "  Merging ${#LANE_BAMS[@]} lane(s) into $OUTBAM"
    samtools merge -f -@ $THREADS $OUTBAM "${LANE_BAMS[@]}"
    samtools index $OUTBAM

    # Remove per-lane intermediates
    rm -f "${LANE_BAMS[@]}"

    BAMS+=("$OUTBAM")
    echo "Finished $SAMPLE"
done

echo "✅ All alignments completed."

# ================================================================
#                   PART 2: JOINT VARIANT CALLING
# ================================================================

# ======== STEP 3: List the BAM files made in this run ========
echo "Samples included in variant calling:"
printf '  %s\n' "${BAMS[@]}"

# ======== STEP 4: Generate pileup and call variants ========
bcftools mpileup -f $REF -Ou -a FORMAT/AD,FORMAT/DP,INFO/AD \
    --threads $THREADS "${BAMS[@]}" | \
bcftools call -mv -Oz -o ${OUTPREFIX}.raw.vcf.gz

# ======== STEP 5: Index the VCF ========
bcftools index ${OUTPREFIX}.raw.vcf.gz

# ======== STEP 6: Basic filtering (optional but recommended) ========
# This filters out low-quality or low-depth sites
bcftools filter -s LOWQUAL -e 'QUAL<30 || INFO/DP<10' \
    ${OUTPREFIX}.raw.vcf.gz -Oz -o ${OUTPREFIX}.filtered.vcf.gz

bcftools index ${OUTPREFIX}.filtered.vcf.gz

echo "✅ Joint variant calling complete!"
echo "Output VCF: ${OUTPREFIX}.filtered.vcf.gz"

