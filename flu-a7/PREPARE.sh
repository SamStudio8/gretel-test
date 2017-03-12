#@CHITIN_META name prepare
#@CHITIN_INPUT 1 read_size
#@CHITIN_INPUT 2 run_path
#@CHITIN_INPUT 3 coverage
#@CHITIN_INPUT 4 fasta

python ../shredder.py --cover $3 $1 $4 > $2/reads.fq

bowtie2 --local -D 20 -R 3 -L 7 -N 1 -p 8 --gbar 1 --mp 2 -x master_982_24.bti -U $2/reads.fq --un $2/unaligned.fq --no-unal -S $2/out.sam

#@CHITIN_START_BLOCK
samtools view -bS $2/out.sam | samtools sort > $2/out.sort.bam
samtools index $2/out.sort.bam
#@CHITIN_END_BLOCK

python ../snpper.py $2/out.sort.bam 'KP404423.1' 982 > $2/raw.vcf

#@CHITIN_START_BLOCK
cp $2/raw.vcf $2/raw.raw.vcf
bgzip $2/raw.vcf
tabix $2/raw.vcf.gz
#@CHITIN_END_BLOCK
