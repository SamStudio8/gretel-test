gretel-test: DHFR
=================

## Prep

    bowtie2-build master.fa master.bti

## Generate New Data

    python chitinise.py
    python upload.py

## Original Data

Original results can be [downloaded from Google Drive](https://drive.google.com/open?id=0B4t7QqgmvVLHSzlpbXVmMUlmVWs).

## Results

    python ../fasta_snpper.py muscle-I20170117-180402-0567-32986293-oy.clustalw.fa > muscle-I20170117-180402-0567-32986293-oy.clustalw.vcf
    python ../collate_hamming.py result/fbc8eada-a2a7-4ede-b1b5-757e3640fa93_final/ muscle-I20170117-180402-0567-32986293-oy.clustalw.fa muscle-I20170117-180402-0567-32986293-oy.clustalw.vcf genes_aln.bed data/fbc8eada-a2a7-4ede-b1b5-757e3640fa93/ > fbc_hamming_wpd.txt
    python ../correspond_uuid.py data/fbc8eada-a2a7-4ede-b1b5-757e3640fa93/fbc8eada-a2a7-4ede-b1b5-757e3640fa93.manifest fbc_hamming_wpd.txt > fbc_hamming_wpd_wmeta.txt

    cd result/fbc8eada-a2a7-4ede-b1b5-757e3640fa93_final/
    for crumb in `ls *.crumbs`; do echo $crumb | sed 's,.crumbs,,'; tail -n+2 $crumb | cut -f3 | sed 's,-,,' | sort -n | awk 'NR == 1 { print }END{ if (NR==0) print "0\t0\n"; else print; }'; done | paste - - - > ../../fbc8eada-a2a7-4ede-b1b5-757e3640fa93_final.minmaxlikl.ls # add add a nice header
    python ../../correspond_uuid.py fbc8eada-a2a7-4ede-b1b5-757e3640fa93_final.minmaxlikl.ls fbc_hamming_wpd_wmeta.txt > fbc_hamming_wpd_wmeta.withbestworstlikl.txt

    cut -f3 fbc_hamming_wpd_wmeta.txt | sed 's,.*__-,,' | sed 's,None,0,' > fbc_hamming_wpd_wmeta.likl.ls #and change header to abslik
    paste fbc_hamming_wpd_wmeta.withbestworstlikl.txt fbc_hamming_wpd_wmeta.likl.ls > fbc_hamming_wpd_wmeta.withbestworstlikl.withabslikl.txt

