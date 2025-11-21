CHR=$1
chunk_idx=$2

OUTDIR=/../vcfs/IMMerge/chunks
mkdir -p $OUTDIR
CHUNK_FILE=/../files/chunks/chr$CHR.chunks.txt

MY_LIST=($(awk '{print $1}'  /../files/imputed_gt/chr$CHR.used.list))
MY_INFO_LIST=($(awk '{print $1}'  /../files/imputed_gt/chr$CHR.used.info.list))

CHUNKS=($(awk '{print $1}' $CHUNK_FILE))
REGION=${CHUNKS[$chunk_idx]}

NUM=${#MY_LIST[@]}

SAMPLE_LIST=""
for (( i = 0; i < $NUM; i++ ));do SAMPLE_LIST="$SAMPLE_LIST ${MY_LIST[i]}";done
echo $SAMPLE_LIST

VCF=${MY_LIST[$SLURM_ARRAY_TASK_ID]}

ETH=$(dirname $VCF | xargs dirname | xargs basename)
COHORT=$(dirname $VCF | xargs basename)

mkdir -p $OUTDIR/$ETH/$COHORT
OUT=$OUTDIR/$ETH/$COHORT/chr$CHR.chunk$chunk_idx.vcf.gz
OUT_INFO=$OUTDIR/$ETH/$COHORT/
bcftools=/qnap-wlee2/wanpingleelab/chengp/tools/bcftools-1.19/bcftools

# chunk vcf file
$bcftools view --threads 2 -r $REGION $VCF \
        | $bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o $OUT
$bcftools index -t $OUT
$bcftools stats $OUT > ${OUT%%.gz}.stats

# split_info
split_info=/../scripts/split_info_by_region.py

if [[ -f $OUT_INFO/chr$CHR.chunk$chunk_idx.inf.gz ]]; then
        rm  $OUT_INFO/chr$CHR.chunk$chunk_idx.inf.gz
fi
INFO=${MY_INFO_LIST[$SLURM_ARRAY_TASK_ID]}
python3 $split_info -i $INFO -o $OUT_INFO --name-style chunk -c chunk$chunk_idx -r $REGION
