#!/bin/bash

for i in `ls | grep test | xargs`
do
TIME=`cat "${i}/out-pml.txt" | grep 'Solve time' | tail -1 | sed 's/ msec.*//g ; s/.* \t//g'`
DIFF=`cat "${i}/out-pml.txt" | grep 'Diff (L2)' | sed 's/.* //g' | xargs`
RESID=`cat "${i}/out-pml.txt" | grep 'V-Cycle Residual' | tail -1 | sed 's/.*\t//'`
CHI=`echo "${i}" | sed 's/[^(]*(// ; s/)[^(]*(/ / ; s/,/ /g ; s/)//g'`
echo "${CHI} ${RESID} ${DIFF} ${TIME}"
done

