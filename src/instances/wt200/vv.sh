#!/bin/sh
s=26
n=1
nj=200
for rdd in 0.2 0.4 0.6 0.8 1.0; do
  for tf in 0.2 0.4 0.6 0.8 1.0; do
    ../../../problem/prob -nj $nj -nm 1 -np 5 -tf $tf -rdd $rdd -seed $s
    s=$(expr $s + 1)
    for f in 001 002 003 004 005; do
      nn=$(printf "%03d" $n)
      echo " $nj" > wt${nj}_${nn}.dat
      tail +7 ${nj}_*_${f}.dat >> wt${nj}_${nn}.dat
      rm -f ${nj}_*_${f}.dat
      n=$(expr $n + 1)
    done
  done
done
