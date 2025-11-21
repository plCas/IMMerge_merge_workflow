CHR=$1

export PATH=/../htslib-1.16/bin:$PATH
source /../miniconda3/etc/profile.d/conda.sh
conda activate IMMerge

target_chunk_num=$SLURM_ARRAY_TASK_ID
#target_chunk_num=$2

DIR=/../vcfs/IMMerge/chunks
mkdir -p $DIR
LISTDIR=/../vcfs/IMMerge/chunks
mkdir -p $LISTDIR
FILES=($(find $LISTDIR -mindepth 2 -type f -iname "chr$CHR.*"|grep .vcf.gz$))
FILES_INFO=($(find $LISTDIR -mindepth 2 -type f -iname "chr$CHR.*"|grep .info.gz$))

SAMPLE_LIST=""
for i in ${FILES[*]};do
    chunk_num=$(echo $i| grep -oP 'chunk\K[0-9]+')
    if [ $chunk_num == $target_chunk_num ];then
        SAMPLE_LIST="$SAMPLE_LIST $i"
    fi
done

echo $SAMPLE_LIST

SAMPLE_INFO_LIST=""
for i in ${FILES_INFO[*]};do
    chunk_num=$(echo $i| grep -oP 'chunk\K[0-9]+')
    if [ $chunk_num == $target_chunk_num ];then
        SAMPLE_INFO_LIST="$SAMPLE_INFO_LIST $i"
    fi
done

echo $SAMPLE_INFO_LIST

python3 /../tools/IMMerge/src/IMMerge/merge_files.py \
    --thread 10 \
    --input $SAMPLE_LIST \
    --info $SAMPLE_INFO_LIST \
    --output $DIR/chr$CHR.chunk$target_chunk_num \
    --check_duplicate_id true \
    --missing 0

python3 /../tools/IMMerge/src/IMMerge/make_info.py \
    --thread 10 \
    --output_dir $DIR \
    --output_fn chr$CHR.chunk$target_chunk_num.info.gz \
    --input $DIR/chr$CHR.chunk$target_chunk_num.vcf.gz

bcftools=/../tools/bcftools-1.19/bcftools
$bcftools index -t $DIR/chr$CHR.chunk$target_chunk_num.vcf.gz
$bcftools stats $DIR/chr$CHR.chunk$target_chunk_num.vcf.gz  > $DIR/chr$CHR.chunk$target_chunk_num.vcf.stats
