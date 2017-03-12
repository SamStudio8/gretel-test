#@CHITIN_META name tree_makereads
#@CHITIN_INPUT 1 read_size
#@CHITIN_INPUT 2 read_uuid
#@CHITIN_INPUT 3 haps_uuid
#@CHITIN_INPUT 4 coverage
#@CHITIN_INPUT 5 exp_uuid

mkdir $5/data/
mkdir $5/data/$2

python ~/Projects/Packages/gretel-test/shredder4.py --cover $4 --sam $5/data/$2/out.sam $1 $5/haps/$3.fa > $5/data/$2/reads.fq
samtools view -bS $5/data/$2/out.sam | samtools sort > $5/data/$2/out.sort.bam
samtools index $5/data/$2/out.sort.bam

python ~/Projects/Packages/gretel-test/snpper.py $5/data/$2/out.sort.bam 'HOOT' 3000 0 > $5/data/$2/raw.vcf

cp $5/data/$2/raw.vcf $5/data/$2/raw.raw.vcf
bgzip $5/data/$2/raw.vcf
tabix $5/data/$2/raw.vcf.gz

