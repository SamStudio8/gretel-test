gretel-test: Flu-A7
===================

## Prep

    bowtie2-build master_982_24.fa master_982_24.bti

## Generate New Data

    python chitinise.py segment7.no99dups.nomaster.only982.24.fa
    python upload.py

## Original Data

Original data and results can be [downloaded from Google Drive](https://drive.google.com/open?id=0B4t7QqgmvVLHSzlpbXVmMUlmVWs).

## Results

    python ../collate_hamming.py result/09d1e503-3e17-455c-bd47-dbaea3f96b3e_final/  segment7.no99dups.nomaster.only982.24.clustalw.alnfa segment7.no99dups.nomaster.only982.24.clustalw.aln.vcf 0 data/09d1e503-3e17-455c-bd47-dbaea3f96b3e/ > collated_hamming.txt
    python ../correspond_uuid.py data/09d1e503-3e17-455c-bd47-dbaea3f96b3e/09d1e503-3e17-455c-bd47-dbaea3f96b3e.manifest collated_hamming.txt > 091_hamming_wmeta.txt
