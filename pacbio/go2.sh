# Ensure index
samtools faidx ../$1.full_hap.fasta

# Get FASTQ from CCS BAM
#bamtools convert -format fastq -in ../../../p1.bam -out pb_ccs_1.fq
#seqkit fq2fa pb_ccs_1.fq > pb_ccs_1.fa

# Align ZMW FASTQ to haplotypes
minimap2 -ax map-pb ../$1.full_hap.fasta pb_ccs_1.fq --secondary=no --MD | samtools sort - -o $1.pbp1-on-g2019.mm2.sorted.bam; samtools index $1.pbp1-on-g2019.mm2.sorted.bam


# Count ZMW
#python count_passes.py > passes.ls

# Parse and grade haplotypes
python haplify_zmw.py --bam $1.pbp1-on-g2019.mm2.sorted.bam --snps snps/$1.all.snps --bed gretel_zmw.bed --infasta ../$1.full_hap.fasta --counts passes.ls --outfasta $1.pb_ccs_1.hap.fasta --onlysnps --onlybed --padstart --qual 40
