#!/bin/bash

# $1 : archea ids
# $2 : input
# $3 : output

nw_reroot $2 4480201 | nw_reroot - `cat $1|tr '\n' ' '` > $3

paste <( nw_ed -n $3 '1 == a' s |nw_stats -|grep leaves ) <( wc -l $1 )
