#@CHITIN_META name tree_harness
#@CHITIN_INPUT 1 prefix
#@CHITIN_INPUT 2 max_species
#@CHITIN_INPUT 3 divergence
#@CHITIN_INPUT 4 tree_name
#@CHITIN_INPUT 5 exp_uuid

#@CHITIN_START_BLOCK TREEIFY
TREE_END=$((67 + $2 - 3))

TREE_TREE="(A:$3,B:$3)"
for i in `seq 67 $TREE_END`; do CURR=`printf "\x$(printf %x $i)"`; TREE_TREE="($TREE_TREE:0,$CURR:$3)"; done
TREE_TREE="$TREE_TREE;"

mkdir $5
mkdir $5/masters/
mkdir $5/trees/
mkdir $5/trees/$4
cp treefile $5/trees/$4/tree.treefile
echo $TREE_TREE >> $5/trees/$4/tree.treefile
#@CHITIN_END_BLOCK

/home/sam/ware/Seq-Gen.v1.3.3/source/seq-gen -k 1 -mGTR < $5/trees/$4/tree.treefile | sed -e '1d' -e 's,^,>,' -e 's,\s\+,\n,' > $5/trees/$4/tree.fa
head -n3 $5/trees/$4/tree.treefile | sed -e '1d' -e 's,^,>,' -e 's,\s\+,\n,' > $5/masters/$4.master.fa
