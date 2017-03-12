#@CHITIN_META name raxml_tree
#@CHITIN_INPUT 1 res_dir
#@CHITIN_INPUT 2 res_uuid
#@CHITIN_INPUT 3 haps_fa

cat $1/$2.out.fasta > $1/$2.mix.fasta
cat $3 >> $1/$2.mix.fasta
mkdir $1/raxml/$2
/home/sam/ware/standard-RAxML/raxmlHPC -w $1/raxml/$2/ -m GTRCAT -n $2 -s $1/$2.mix.fasta -p $RANDOM
/home/sam/Projects/Packages/gretel-test/treedist_all $1/raxml/$2/RAxML_result.$2 vector
cat $1/raxml/$2/RAxML_result.$2.dist | awk '{print $ 1"\t"$ 2"\t"$ 3; print $ 2"\t"$ 1"\t"$ 3}' | grep "^[A-Z]" | grep "__" | sed 's/__/	/g' | sort -k3,2 -k4,4n | awk 'BEGIN{prev=""};{if($ 2!=prev){prev=$ 2;print $ 0}}' > $1/dists/$2
