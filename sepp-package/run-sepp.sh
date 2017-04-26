#!/bin/bash  

if [ $# -lt 2 ]; then
   echo USAGE: $0 "[input fragments file] [output prefix] [optional: -x number-of-cores ] [optional: -A alignment subset size] [optional: -P placement subset size] [optional: any other SEPP argument]
   Optional commands need not be in order. Any SEPP option can also be passed. For example, use
   -x 8
   to make SEPP us 8 threads"
   exit 1;
fi

# Should point to a (semi) permanent tmp for keeping the important parts of the results 
tmp=/oasis/scratch/$USER/temp_project/sepp/$PBS_JOBID
tmp=`mktemp -t sepp-temp-XXXX`

# Should point to a fast (hopefully ssd) tmp location which may be removed after the run
tmpssd=/scratch/$USER/$PBS_JOBID
tmpssd=`mktemp -t sepp-tempssd-XXXX`

# Input sequence file
f=$1
shift

# Name of the output file prefix 
name=$1
shift

# SEPP placement and alignment subset sizes
p=5000
a=1000

opts=""

while [[ $# -gt 1 ]]
do
	key="$1"

	case $key in
		-A|--alignmentSize)
			a="$2"
			shift # past argument
			;;
		-P|--placementSize)
			p="$2"
			shift # past argument
			;;
		*)
			opts="$opts"" ""$key"" ""$2"
			shift # past argument
			;;
	esac
	shift # past argument or value
done


DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Reference tree
t="$DIR/ref/reference-gg-raxml-bl-rooted-relabelled.tre"
# Reference alignment
alg="$DIR/ref/gg_13_5_ssu_align_99_pfiltered.fasta"
# RAxML info file generated when creating the reference RAxML tree
rxi="$DIR/ref/RAxML_info-reference-gg-raxml-bl.info"

set -e

python $DIR/sepp/run_sepp.py -P $p -A $a -t $t -a $alg -r $rxi -f $f -cp $tmpssd/chpoint-$name -o $name $opts -d $tmp/ -p $tmpssd 1>sepp-$name-out.log 2>sepp-$name-err.log

tail sepp-$name-*

cp -r $tmpssd $tmp;

cp $tmp/${name}_placement.json .
cp $tmp/${name}_rename-json.py .

$DIR/sepp/.sepp/bundled-v3.2/guppy tog ${name}_placement.json

cat ${name}_placement.tog.tre | python ${name}_rename-json.py > ${name}_placement.tog.relabelled.tre

$DIR/sepp/.sepp/bundled-v3.2/guppy tog --xml ${name}_placement.json

cat ${name}_placement.tog.xml | python ${name}_rename-json.py > ${name}_placement.tog.relabelled.xml

echo output files are at ${name}_placement.* and more files are at $tmp . Consider removing $tmp if its files are not needed. 
