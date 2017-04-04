#!/bin/bash  

if [ $# -lt 2 ]; then
   echo USAGE: $0 "[input fragments file] [output prefix] [optional: alignment subset size] [optional: placement subset size]"
   exit 1;
fi

# Should point to a (semi) permanent tmp for keeping the important parts of the results 
tmp=/oasis/scratch/$USER/temp_project/sepp/$PBS_JOBID
tmp=`mktemp -d sepp-temp-XXXX`

# Should point to a fast (hopefully ssd) tmp location which may be removed after the run
tmpssd=/scratch/$USER/$PBS_JOBID
tmpssd=`mktemp -d sepp-tempssd-XXXX`

# Input sequence file
f=$1

# Name of the output file prefix 
name=$2

# SEPP placement and alignment subset sizes
p=5000
a=1000
test $# -lt 3 || a=$3
test $# -lt 4 || p=$4

# Leave blank to choose all available CPUs. Otherwise set to e.g., -x 16
#cpus=-x 16
cpus=""

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Reference tree
t="$DIR/ref/reference-gg-raxml-bl-rooted-relabelled.tre"
# Reference alignment
alg="$DIR/ref/gg_13_5_ssu_align_99_pfiltered.fasta"
# RAxML info file generated when creating the reference RAxML tree
rxi="$DIR/ref/RAxML_info-reference-gg-raxml-bl.info"

set -e

python3 $DIR/sepp/run_sepp.py -P $p -A $a -t $t -a $alg -r $rxi -f $f $cpus -cp $tmpssd/chpoint-$name -o $name -d $tmp/ -p $tmpssd 1>sepp-$name-out.log 2>sepp-$name-err.log

tail sepp-$name-*

cp -r $tmpssd $tmp;

cp $tmp/${name}_placement.json .
cp $tmp/${name}_rename-json.py .

$DIR/sepp/.sepp/bundled-v3.2/guppy tog ${name}_placement.json

cat ${name}_placement.tog.tre | python3 ${name}_rename-json.py > ${name}_placement.tog.relabelled.tre

$DIR/sepp/.sepp/bundled-v3.2/guppy tog --xml ${name}_placement.json

cat ${name}_placement.tog.xml | python3 ${name}_rename-json.py > ${name}_placement.tog.relabelled.xml

echo output files are at ${name}_placement.* and more files are at $tmp
