export PATH=/../tools/htslib-1.16/bin:$PATH
source /../miniconda3/etc/profile.d/conda.sh
conda activate IMMerge

CHR=$1

LIST=($(awk '{print $2}' /../files/imputed/chr$CHR.list))

OUTDIR=/../vcfs/imputed_gt

VCF=${LIST[$SLURM_ARRAY_TASK_ID]}
ETH=$(dirname $VCF | xargs dirname | xargs basename)
COHORT=$(dirname $VCF | xargs basename)

bcftools=/../tools/bcftools-1.19/bcftools

mkdir -p $OUTDIR/$ETH/$COHORT
OUT=$OUTDIR/$ETH/$COHORT/chr$CHR.gt.vcf.gz
TMP=$OUTDIR/$ETH/$COHORT/chr$CHR.tmp

$bcftools view --threads 2 $VCF \
        | $bcftools annotate -x ^FORMAT/GT -Ou \
        | $bcftools sort -Oz --temp-dir $TMP -o $OUT
$bcftools index -t $OUT

$bcftools stats $OUT > ${OUT%%.gz}.stats

python3 /../tools/IMMerge/src/IMMerge/make_info.py \
        --output_dir $OUTDIR/$ETH/$COHORT \
        --input $OUT
