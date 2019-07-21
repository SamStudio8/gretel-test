#/home/ubuntu/ware/pb/smrtlink/install/smrtlink-release_6.0.0.47841/bundles/smrttools/install/smrttools-release_6.0.0.47835/smrtcmds/bin/ccs p1.bam ../../../pub33/pbio3/r54118_20190408_081407/3_C01/m54118_190409_113909.subreads.bam --minReadScore=0.65 --minLength=100 --minPasses=20
#samtools view pb_ccs_1.bam | awk '{print ">"$1"\n"$10}' > pb_ccs_1.fa
#python count_passes.py > passes.ls

python mask_qual.py > pb_ccs_1.fa

~/git/minimap2/minimap2 -ax map-pb ../ref-seq.contigs.g.fa pb_ccs_1.fa --secondary=no --sam-hit-only > pbp1-on-g2019.mm2.sam

samtools sort pbp1-on-g2019.mm2.sam -o pbp1-on-g2019.mm2.sorted.bam; samtools index pbp1-on-g2019.mm2.sorted.bam

python drop_inserts.py pbp1-on-g2019.mm2.sorted.bam ../ref-seq.contigs.g.fa 'R123' snps/123.snps.3 --onlysnps --padstart --bed gretel_zmw.bed | seqkit rename -n | seqkit rmdup -s -D 123.compressed.rmdup.count > 123.compressed.rmdup.fa
python rename_ccs.py 123.compressed.rmdup.fa 123.compressed.rmdup.count 1 > 123.compressed.rmdup.count.fa

python drop_inserts.py pbp1-on-g2019.mm2.sorted.bam ../ref-seq.contigs.g.fa 'R90' snps/90.snps.3 --onlysnps --padstart  --bed gretel_zmw.bed  | seqkit rename -n | seqkit rmdup -s -D 90.compressed.rmdup.count > 90.compressed.rmdup.fa
python rename_ccs.py 90.compressed.rmdup.fa 90.compressed.rmdup.count 1 > 90.compressed.rmdup.count.fa

python drop_inserts.py pbp1-on-g2019.mm2.sorted.bam ../ref-seq.contigs.g.fa 'R31' snps/31.snps.3 --onlysnps --padstart  --bed gretel_zmw.bed  | seqkit rename -n | seqkit rmdup -s -D 31.compressed.rmdup.count > 31.compressed.rmdup.fa
python rename_ccs.py 31.compressed.rmdup.fa 31.compressed.rmdup.count 1 > 31.compressed.rmdup.count.fa

rm all.compressed.rmdup.count.fa
rm all.compressed.rmdup.count
cat *.compressed.rmdup.count.fa > all.compressed.rmdup.count.fa
cat *.compressed.rmdup.count > all.compressed.rmdup.count
samtools faidx all.compressed.rmdup.count.fa
samtools faidx ../$1.snp_hap.fasta

/cephfs/lomanlabz/sam/bin/ncbi-blast-2.9.0+-src/c++/ReleaseMT/bin/makeblastdb -dbtype nucl -in ../$1.snp_hap.fasta
/cephfs/lomanlabz/sam/bin/ncbi-blast-2.9.0+-src/c++/ReleaseMT/bin/blastn -db ../$1.snp_hap.fasta -query all.compressed.rmdup.count.fa  -num_threads 8 -outfmt 6 -num_alignments 1 > $1.snp_hap.b6
sed 's,/,\t,' $1.snp_hap.b6 | sed 's,__,\t,' | sed 's,/compressed_,\t,' | sort -k1,1 -k5,5nr > $1.snp_hap.b6.tsv

/cephfs/lomanlabz/sam/bin/ncbi-blast-2.7.1+/bin/blastn -db ../$1.snp_hap.fasta -query all.compressed.rmdup.count.fa  -num_threads 72 -num_alignments 1 > $1.snp_hap.b1

python results_table.py | sort -k1,1 -k8,8nr | column -t > result
grep 'G123_0' result | awk '{print $5}' > result.g123-0.repzmw.ls
seqkit grep -f result.g123-0.repzmw.ls pb_ccs_1.fa > zmw-check/result.g123-0.repzmw.fa
