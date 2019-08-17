#!/bin/bash  

if [ $# -lt 2 ]; then
   echo USAGE: $0 "[input fragments file] [output prefix] [optional: -x number-of-cores ] [optional: -A alignment subset size] [optional: -P placement subset size] [optional: any other SEPP argument] [optional: -t filename reference phylogeny] [optional: -a filename reference alignment] [optional: -r filename reference RAxML info file] [optional: -n 1 = no tree-, just placements- computation] [optional: -b 1 = report debugging information ]
   Optional commands need not be in order. Any SEPP option can also be passed. For example, use
   -x 8
   to make SEPP us 8 threads"
   exit 1;
fi

tmpbase=`mktemp -d -t 'tmp.XXXXXXXXXX'`
if [ $? -ne 0 ]; 
then
	echo "$0: Can't create temp directory, exiting..."
	exit 1
fi
export TMPDIR=$tmpbase

# Should point to a (semi) permanent tmp for keeping the important parts of the results 
tmp=`mktemp -d -t 'sepp-tmp-XXXXXXXXXX'`

if [ -z ${TMPDIRSSD+x} ];
then
    # if $TMPDIRSSD does not exist, just create under the main tmp
	tmpssd=`mktemp -d -t sepp-tempssd-XXXXXXXXXX`
else
	# Should point to a fast (hopefully ssd) tmp location which may be removed after the run
    export TMPDIR=${TMPDIRSSD}
	tmpssd=`mktemp -d -t sepp-tempssd-XXXXXXXXXX`

    # make sure the regular temp is reset so other consumers of the pipeline can use it
    export TMPDIR=$tmpbase
fi

# from http://stackoverflow.com/a/2130323
function cleanup {
  exitcode=`echo $?`
  if [ $exitcode != 0 ] && [ ! -z "$printDebug" ];
  then
    echo "========= Execution of SEPP failed with exit code $exitcode =================";
    echo "temporary working directories are NOT deleted for further inspection:";
    echo "  \$tmp = $tmp";
    echo "  \$tmpssd = $tmpssd";
    echo "--------- Content of STDOUT -----------------------------------------";
    cat sepp-$name-out.log
    echo "--------- Content of STDERR -----------------------------------------";
    cat sepp-$name-err.log
    echo "=====================================================================";
  else
    echo "Removing $tmp";
    rm -r $tmp
    rm -r $tmpssd
  fi
  unset TMPDIR
}
trap cleanup EXIT

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
		-a|--referenceAlignment)
			alg="$2"
			shift # past argument
			;;
		-t|--referencePhylogeny)
			t="$2"
			shift # past argument
			;;
		-r|--referenceInfofile)
			rxi="$2"
			shift # past argument
			;;
		-n|--noTreeComputation)
			noTree="$2"
			shift # past argument
			;;
		-b|--debugInformation)
			printDebug="$2"
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
if [ -z $t ]; then
	# resort to the default reference tree, if not specified via command line argument
	t="$DIR/ref/reference-gg-raxml-bl-rooted-relabelled.tre"
fi;
# Reference alignment
if [ -z $alg ]; then
	# resort to the default reference alignment, if not specified via command line argument
	alg="$DIR/ref/gg_13_5_ssu_align_99_pfiltered.fasta"
fi;
# RAxML info file generated when creating the reference RAxML tree
if [ -z $rxi ]; then
	# resort to the default reference RAxML info file, if not specified via command line argument
	rxi="$DIR/ref/RAxML_info-reference-gg-raxml-bl.info"
fi;

set -e

if [ ! -z "$printDebug" ]; then
	export SEPP_DEBUG=True
fi;
python $DIR/sepp/run_sepp.py -P $p -A $a -t $t -a $alg -r $rxi -f $f -o $name $opts -d $tmp/ -p $tmpssd 1>sepp-$name-out.log 2>sepp-$name-err.log

tail sepp-$name-*

cp -r $tmpssd $tmp;

cp $tmp/${name}_placement.json .
cp $tmp/${name}_rename-json.py .

# we might want to split computation in two parts: a) obtaining placements and b) creation of an insertion tree.
# If -n set to something, we stop after a) and leave it to the user to compute b) afterwards.
if [ -z ${noTree+x} ]; then
	gbin=$( dirname `grep -A1 "pplacer" $DIR/sepp/.sepp/main.config |grep path|sed -e "s/^path=//g"` )

	$gbin/guppy tog ${name}_placement.json

	cat ${name}_placement.tog.tre | python ${name}_rename-json.py > ${name}_placement.tog.relabelled.tre

	$gbin/guppy tog --xml ${name}_placement.json

	cat ${name}_placement.tog.xml | python ${name}_rename-json.py > ${name}_placement.tog.relabelled.xml
else
	echo "User requested skipping of insertion tree computation. Only placements are returned.";
fi;

echo output files are at ${name}_placement.* and more files are at $tmp . Consider removing $tmp if its files are not needed. 
