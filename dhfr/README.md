gretel-test: DHFR
=================

## Prep

    bowtie2-build master.fa master.bti

## Generate New Data

    python chitinise.py
    python upload.py

## Original Data

    tar xvf data_fbc8eada-a2a7-4ede-b1b5-757e3640fa93.tar
    tar xvf result_fbc8eada-a2a7-4ede-b1b5-757e3640fa93.tar

## Results

    python ../fasta_snpper.py muscle-I20170117-180402-0567-32986293-oy.clustalw.fa > muscle-I20170117-180402-0567-32986293-oy.clustalw.vcf
    python ../collate_hamming.py result/fbc8eada-a2a7-4ede-b1b5-757e3640fa93_final/ muscle-I20170117-180402-0567-32986293-oy.clustalw.fa muscle-I20170117-180402-0567-32986293-oy.clustalw.vcf genes_aln.bed data/fbc8eada-a2a7-4ede-b1b5-757e3640fa93/ > fbc_hamming_wpd.txt
    python ../correspond_uuid.py data/fbc8eada-a2a7-4ede-b1b5-757e3640fa93/fbc8eada-a2a7-4ede-b1b5-757e3640fa93.manifest fbc_hamming_wpd.txt > fbc_hamming_wpd_wmeta.txt
