#!/bin/bash

for i in `ls | grep test | xargs`
do
TIME=`cat "${i}/out-pml.txt" | grep 'Solve time' | tail -1 | sed 's/ msec.*//g ; s/.* \t//g'`
DIFF=`cat "${i}/out-pml.txt" | grep 'Diff (L2)' | sed 's/.* //g' | xargs`
CHI=`echo "${i}" | sed 's/[^(]*(// ; s/)[^(]*(/ / ; s/,/ /g ; s/)//g'`
echo "${CHI} ${DIFF} ${TIME}"
done

