# gretel-test: Treeviomes (Synthetic Metahaplomes)

## Generate New Data

    python chitinise.py

## Original Data

### Data
Full data is unfortunately only available on request due to its storage impracticality (>300 GB)... But, you can generate your own data using our protocol with `chitin`!
We do however provide the data tar that includes the trees, haplotypes and masters, [on our Google Drive](https://drive.google.com/open?id=0B4t7QqgmvVLHSzlpbXVmMUlmVWs).

### Results

Original results can be [downloaded from Google Drive](https://drive.google.com/open?id=0B4t7QqgmvVLHSzlpbXVmMUlmVWs).

## Results
### Hamming

    python collate4.py results/710dd289-621e-456c-9274-e080b4a5b5fd 710dd289-621e-456c-9274-e080b4a5b5fd/710dd289-621e-456c-9274-e080b4a5b5fd.manifest.chitin.FIXED 710dd289-621e-456c-9274-e080b4a5b5fd > collate_wmuscle.txt


### Dists

    python result_chitinise.py results/710dd289-621e-456c-9274-e080b4a5b5fd 710dd289-621e-456c-9274-e080b4a5b5fd/710dd289-621e-456c-9274-e080b4a5b5fd.manifest.chitin.FIXED 710dd289-621e-456c-9274-e080b4a5b5fd
