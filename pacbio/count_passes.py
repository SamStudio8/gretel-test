import pysam
bam_fh = pysam.AlignmentFile("../../../p1.bam", check_sq=False)
for read in bam_fh.fetch(until_eof=True):
    print(read.query_name, read.get_tag("np"))
