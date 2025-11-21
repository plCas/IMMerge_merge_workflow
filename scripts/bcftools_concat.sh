CHR=$1
threads=12

DIR=/../vcfs/IMMerge/chunks
OUT=/../vcfs/IMMerge
mkdir -p $OUT

VCFS=$(ls -v $DIR/chr$CHR.chunk*|grep vcf.gz$)

bcftools=/../tools/bcftools-1.19/bcftools
$bcftools concat --threads $threads --allow-overlaps --remove-duplicates $VCFS -Oz -o $OUT/chr$CHR.vcf.gz
$bcftools index --threads $threads -t $OUT/chr$CHR.vcf.gz
$bcftools stats $OUT/chr$CHR.vcf.gz > $OUT/chr$CHR.vcf.stats
