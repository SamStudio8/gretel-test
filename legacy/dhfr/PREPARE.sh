#@CHITIN_META name prepare
#@CHITIN_INPUT 1 read_size
#@CHITIN_INPUT 2 run_path
#@CHITIN_INPUT 3 coverage

python ../shredder.py --error 1 --cover $3 $1 genes.aln.fa > $2/reads.fq

bowtie2 --local -D 20 -R 3 -L 3 -N 1 -p 8 --gbar 1 --mp 3 -x master.bti -U $2/reads.fq --un $2/unaligned.fq --no-unal -S $2/out.sam

#@CHITIN_START_BLOCK
samtools view -bS $2/out.sam > $2/out.bam
samtools sort $2/out.bam > $2/out.sort.bam
samtools index $2/out.sort.bam
#@CHITIN_END_BLOCK

python ../snpper.py $2/out.sort.bam 'gi|161089328|gb|EU145592.1|' 564 > $2/raw.vcf

#@CHITIN_START_BLOCK
cp $2/raw.vcf $2/raw.raw.vcf
bgzip $2/raw.vcf
tabix $2/raw.vcf.gz
#@CHITIN_END_BLOCK
