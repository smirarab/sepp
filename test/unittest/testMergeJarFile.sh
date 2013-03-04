#!/bin/sh
set -x

java -jar ../../tools/merge/seppJsonMerger.jar - - data/tmp/test$2.json -s -r 4 -p 0 -t data/mock/pyrg/all_taxon.taxonomy -m data/mock/pyrg/species.mapping -c data/tmp/.t1 $1 < data/tmp/test.merge.stdin

sort data/tmp/.t1 > data/tmp/test$2.classification.txt
