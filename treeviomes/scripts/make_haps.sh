#@CHITIN_META name tree_makehaps
#@CHITIN_INPUT 1 num_species
#@CHITIN_INPUT 2 tree_name
#@CHITIN_INPUT 3 hap_path
#@CHITIN_INPUT 4 exp_uuid

#@CHITIN_START_BLOCK
mkdir $4/haps
FA_LINES=$(($1 * 2))
head -n $FA_LINES $4/trees/$2/tree.fa > $4/haps/$3.fa
#@CHITIN_END_BLOCK
