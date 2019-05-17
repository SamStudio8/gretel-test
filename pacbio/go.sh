python drop_inserts.py pbp1-on-g2019.mm2.sorted.bam 2019.g.out.top.fa 'G123_0__16.65' snps/123.snps.3 | seqkit rename -n | seqkit rmdup -s -D 123.compressed.rmdup.count > 123.compressed.rmdup.fa
python rename_ccs.py 123.compressed.rmdup.fa 123.compressed.rmdup.count 10 > 123.compressed.rmdup.count.fa
