#!/bin/sh

cat /dev/null > POTFILES.in

for f in ../lib/src/*.c ; do
   echo $f | sed s+\.\./++ >> POTFILES.in
done

for f in ../cli/*.c ; do
   echo $f | sed s+\.\./++ >> POTFILES.in
done   

for f in `find ../gui -maxdepth 1 -name "*.c" -type f` ; do
   echo $f | sed s+\.\./++ >> POTFILES.in
done 

for f in ../gui2/*.c ; do
   echo $f | sed s+\.\./++ >> POTFILES.in
done

for f in ../plugin/*.c ; do
   echo $f | sed s+\.\./++ >> POTFILES.in
done 
