gretel-test: 5VM
================

## Data
`SRR961514`

## Result Processing with BLAST

    makeblastdb -in strains.fa -out strains.blast -dbtype nucl
    python process_hiv.py <gene> <start_clip> <end_clip> > <gene>/out.minus.fasta
    blastn -db strains.blast -query <gene>/out.minus.fasta -outfmt 6 -max_target_seqs 1 -max_hsps 1 > <gene>.list
    python parse_hiv_blast.py <gene>.list <gene> >> hiv.tab
    cut -f1,3,13,14 hiv.tab > hiv.c.tab

