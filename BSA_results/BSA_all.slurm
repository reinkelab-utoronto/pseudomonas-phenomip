#!/bin/bash
#SBATCH --job-name=BSA_all
#SBATCH --array=0-5
#SBATCH --output=BSA_%A_%a.out
#SBATCH --error=BSA_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --nodes=1

# ==============================================================================
# Bulk segregant analysis: reads -> allele depths, one SLURM array task per
# population, all populations processed with the same pipeline:
#
#   merge lanes -> fastp -> bwa mem (with read groups) -> sort -> Picard
#   MarkDuplicates -> bcftools mpileup (-a FORMAT/AD,FORMAT/DP, -d 2000)
#   -> bcftools call -mv -> filter (QUAL >= 20, per-sample DP >= 10)
#   -> allele depths
#
# Populations (array index):
#   0 JU2610  MYb114  uninfected / 33% ("low") / 100% ("full")
#   1 ED3005  MYb114  uninfected / 33% ("low") / 100% ("full")
#   2 JU2522  PA14    G6: uninfected / 10%;  G8: uninfected / 50%
#   3 WN2050  PA14    G6: uninfected / 10%;  G8: uninfected / 50%
#   4 JU2250  MYb395  uninfected / 100% ("m")
#   5 EG4946  MYb395  uninfected / 100% ("m")
#
# Usage (from the final_BSA folder):
#   sbatch BSA_all.slurm              # all six populations
#   sbatch --array=4,5 BSA_all.slurm  # only some of them
#
# Inputs are read from the existing population folders and are never
# modified. Everything is written to rerun/<population>/.
#
# Output per population (same format as before, so the R analysis scripts
# work unchanged):
#   <POP>_allele_depths.tsv   CHROM, POS, then AD (ref,alt) per sample, no header,
#                             samples in the order listed in <POP>_samples.txt
#   <POP>_samples.txt         sample order of the AD columns
#   <POP>_samples.vcf.gz      filtered calls (+ raw .bcf, unfiltered)
#   <POP>_versions.txt        tool versions used for this run
#   <POP>_BSA_results.tar.gz  the files to download (see the end of the script)
# ==============================================================================

set -euo pipefail

BASE=/home/yinchen/scratch/final_BSA
REF=${BASE}/reference/20250625_c_elegans_WS283.genome.fa
THREADS=${SLURM_CPUS_PER_TASK:-16}

MAX_DEPTH=2000   # per-sample depth cap for mpileup (default 250 is too low for pools)
MIN_QUAL=20      # minimum variant QUAL
MIN_DP=10        # minimum per-sample depth

module load fastp
module load bwa
module load samtools
module load bcftools
module load picard

# ------------------------------------------------------------------------------
# Population settings. Each sample is "label:input_prefix"; the input FASTQs are
# <input_prefix>_L00*_R1_001.fastq.gz and _R2_ (all lanes are merged). Samples
# are listed in the column order of the allele-depth file used previously.
# ------------------------------------------------------------------------------
POPS=(JU2610 ED3005 JU2522 WN2050 JU2250 EG4946)
POP=${POPS[$SLURM_ARRAY_TASK_ID]}

case $POP in
  JU2610)
    INDIR=${BASE}/JU2610
    SAMPLES=(full1:JU2610_full1_S213 full2:JU2610_full2_S214
             low1:JU2610_low1_S211   low2:JU2610_low2_S212
             uni1:JU2610_uni1_S209   uni2:JU2610_uni2_S210) ;;
  ED3005)
    INDIR=${BASE}/ED3005
    SAMPLES=(full1:ED3005_full1_S207 full2:ED3005_full2_S208
             low1:ED3005_low1_S205   low2:ED3005_low2_S206
             uni1:ED3005_uni1_S204   uni2:ED3005_uni2_S215) ;;
  JU2522)
    INDIR=${BASE}/JU2522_BSA
    SAMPLES=(G6_10percentPA1:JU2522_BSA_G6_10percentPA1_S14
             G6_10percentPA2:JU2522_BSA_G6_10percentPA2_S15
             G6_uni1:JU2522_BSA_G6_uni1_S1
             G6_uni2:JU2522_BSA_G6_uni2_S13
             G8_50percentPA1:JU2522_BSA_G8_50percentPA1_S4
             G8_50percentPA2:JU2522_BSA_G8_50percentPA2_S5
             G8_uni1:JU2522_BSA_G8_uni1_S2
             G8_uni2:JU2522_BSA_G8_uni2_S3) ;;
  WN2050)
    INDIR=${BASE}/WN2050_BSA
    SAMPLES=(G6_10percentPA1:WN2050_BSA_G6_10percentPA1_S8
             G6_10percentPA2:WN2050_BSA_G6_10percentPA2_S9
             G6_uni1:WN2050_BSA_G6_uni1_S6
             G6_uni2:WN2050_BSA_G6_uni2_S7
             G8_50percentPA1:WN2050_BSA_G8_50percentPA1_S11
             G8_50percentPA2:WN2050_BSA_G8_50percentPA2_S12
             G8_uni1:WN2050_BSA_G8_uni1_S16
             G8_uni2:WN2050_BSA_G8_uni2_S10) ;;
  JU2250)
    INDIR=${BASE}/JU2250_BSA
    SAMPLES=(m1:JU_BSA_m1_S98     m2:JU_BSA_m2_S99
             uni1:JU_BSA_uni1_S96 uni2:JU_BSA_uni2_S97) ;;
  EG4946)
    INDIR=${BASE}/EG_BSA
    SAMPLES=(m1:EG_BSA_m1_S94     m2:EG_BSA_m2_S95
             uni1:EG_BSA_uni1_S92 uni2:EG_BSA_uni2_S93) ;;
esac

OUTDIR=${BASE}/rerun/${POP}
mkdir -p ${OUTDIR}/merged
cd ${OUTDIR}
echo "=== ${POP}: reading from ${INDIR}, writing to ${OUTDIR}"

# ------------------------------------------------------------------------------
# 0. Checks: reference index and input FASTQs
# ------------------------------------------------------------------------------
for ext in "" .amb .ann .bwt .pac .sa .fai; do
  if [ ! -f ${REF}${ext} ]; then
    echo "ERROR: missing ${REF}${ext}. Set up the reference folder first (see README)."
    exit 1
  fi
done

for s in "${SAMPLES[@]}"; do
  prefix=${s#*:}
  for read in R1 R2; do
    if ! ls ${INDIR}/${prefix}_L00*_${read}_001.fastq.gz >/dev/null 2>&1; then
      echo "ERROR: no ${INDIR}/${prefix}_L00*_${read}_001.fastq.gz"
      exit 1
    fi
  done
done

# Record tool versions
{
  echo "population: ${POP}"; date
  echo "fastp:    $(fastp --version 2>&1)"
  echo "bwa:      $(bwa 2>&1 | grep -i '^Version')"
  echo "samtools: $(samtools --version | head -1)"
  echo "bcftools: $(bcftools --version | head -1)"
  echo "picard:   $(java -jar $EBROOTPICARD/picard.jar MarkDuplicates --version 2>&1 | tail -1)"
  echo "parameters: MAX_DEPTH=${MAX_DEPTH} MIN_QUAL=${MIN_QUAL} MIN_DP=${MIN_DP}"
} > ${POP}_versions.txt

BAMLIST=${POP}_bams.txt
: > ${BAMLIST}

for s in "${SAMPLES[@]}"; do
  label=${s%%:*}
  prefix=${s#*:}
  echo "=== ${POP} ${label} (${prefix})"

  # ----------------------------------------------------------------------------
  # 1. Merge lanes (sorted, so lanes are always concatenated in the same order)
  # ----------------------------------------------------------------------------
  for read in R1 R2; do
    out=merged/${label}_${read}.fastq.gz
    if [ ! -f ${out} ]; then
      cat $(ls ${INDIR}/${prefix}_L00*_${read}_001.fastq.gz | sort) > ${out}
    fi
  done

  # ----------------------------------------------------------------------------
  # 2. Trim
  # ----------------------------------------------------------------------------
  fastp \
    -i merged/${label}_R1.fastq.gz -I merged/${label}_R2.fastq.gz \
    -o ${label}_R1.trimmed.fastq.gz -O ${label}_R2.trimmed.fastq.gz \
    -h ${label}_fastp.html -j ${label}_fastp.json \
    --thread ${THREADS}

  # ----------------------------------------------------------------------------
  # 3. Align with a read group, so the sample name in the VCF is the label
  # ----------------------------------------------------------------------------
  RG="@RG\tID:${label}\tSM:${label}\tLB:${label}\tPL:ILLUMINA"
  bwa mem -t ${THREADS} -R "${RG}" ${REF} \
      ${label}_R1.trimmed.fastq.gz ${label}_R2.trimmed.fastq.gz \
    | samtools sort -@ ${THREADS} -o ${label}.sorted.bam -
  samtools index ${label}.sorted.bam

  # ----------------------------------------------------------------------------
  # 4. Mark PCR duplicates (mpileup skips reads flagged as duplicates)
  # ----------------------------------------------------------------------------
  java -Xmx16g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I=${label}.sorted.bam O=${label}.dedup.bam M=${label}.dup_metrics.txt
  samtools index ${label}.dedup.bam

  echo ${label}.dedup.bam >> ${BAMLIST}
  rm -f ${label}_R1.trimmed.fastq.gz ${label}_R2.trimmed.fastq.gz
done

# ------------------------------------------------------------------------------
# 5. Joint variant calling, in the sample order of ${BAMLIST}
# ------------------------------------------------------------------------------
bcftools mpileup -Ou -a FORMAT/AD,FORMAT/DP -d ${MAX_DEPTH} -f ${REF} -b ${BAMLIST} \
  | bcftools call -mv -Ob -o ${POP}_samples.raw.bcf
bcftools index ${POP}_samples.raw.bcf

# ------------------------------------------------------------------------------
# 6. Filter: drop sites with QUAL < MIN_QUAL or any sample with DP < MIN_DP
# ------------------------------------------------------------------------------
bcftools filter -e "QUAL<${MIN_QUAL} || FMT/DP<${MIN_DP}" -Ob -o ${POP}_samples.bcf ${POP}_samples.raw.bcf
bcftools index ${POP}_samples.bcf
bcftools view ${POP}_samples.bcf -Oz -o ${POP}_samples.vcf.gz
bcftools index ${POP}_samples.vcf.gz

# ------------------------------------------------------------------------------
# 7. Allele depths (same format as before) and the sample order
# ------------------------------------------------------------------------------
bcftools query -f '%CHROM\t%POS[\t%AD]\n' ${POP}_samples.vcf.gz > ${POP}_allele_depths.tsv
bcftools query -l ${POP}_samples.vcf.gz > ${POP}_samples.txt

echo "Sample order of the AD columns:"; cat ${POP}_samples.txt
echo "Sites in ${POP}_allele_depths.tsv: $(wc -l < ${POP}_allele_depths.tsv)"

# ------------------------------------------------------------------------------
# 8. Package the files to download
# ------------------------------------------------------------------------------
cp ${BASE}/BSA_all.slurm .
LOGS=$(ls ${BASE}/BSA_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${BASE}/BSA_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err 2>/dev/null || true)   # logs so far
# Snapshot the logs: tar fails if it reads a log the job is still writing to
LOGCOPIES=""
for f in ${LOGS}; do cp "$f" . && LOGCOPIES="${LOGCOPIES} $(basename "$f")"; done
tar -czf ${POP}_BSA_results.tar.gz \
  ${POP}_allele_depths.tsv ${POP}_samples.txt ${POP}_versions.txt \
  ${POP}_samples.raw.bcf ${POP}_samples.bcf ${POP}_samples.vcf.gz \
  *_fastp.json *.dup_metrics.txt BSA_all.slurm ${LOGCOPIES}

echo "=== ${POP} finished: ${OUTDIR}/${POP}_BSA_results.tar.gz"

