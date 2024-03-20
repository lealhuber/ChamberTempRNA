#!/bin/bash

while read matrixFileNames2
do
    echo $matrixFileNames2
    awk '$2 > 0 {n=n+1} END {print n}' $matrixFileNames2 ;
done < matrixFileNames2
exit 0