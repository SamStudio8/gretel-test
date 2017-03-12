# gretel-test: Treeviomes (Synthetic Metahaplomes)

## Generate New Data

    python chitinise.py

## Original Data

### Data
Data tar includes the trees, haplotypes and masters.
Full data is available on request due to its storage impracticality (>300 GB)... But, you can generate your own data using our protocol with `chitin`!

    tar xvf 710dd289-621e-456c-9274-e080b4a5b5fd.tar

### Results

    cat result_710dd289-621e-456c-9274-e080b4a5b5fd.tar.* > result_710dd289-621e-456c-9274-e080b4a5b5fd.tar
    tar xvf result_710dd289-621e-456c-9274-e080b4a5b5fd.tar

## Results
### Hamming

    python collate4.py results/710dd289-621e-456c-9274-e080b4a5b5fd 710dd289-621e-456c-9274-e080b4a5b5fd/710dd289-621e-456c-9274-e080b4a5b5fd.manifest.chitin.FIXED 710dd289-621e-456c-9274-e080b4a5b5fd > collate_wmuscle.txt


### Dists

    python result_chitinise.py results/710dd289-621e-456c-9274-e080b4a5b5fd 710dd289-621e-456c-9274-e080b4a5b5fd/710dd289-621e-456c-9274-e080b4a5b5fd.manifest.chitin.FIXED 710dd289-621e-456c-9274-e080b4a5b5fd
