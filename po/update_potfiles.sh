#!/bin/sh

cat /dev/null > POTFILES.in

for f in ../lib/src/*.c ; do
   echo $f | sed s+\.\./++ >> POTFILES.in
done

for f in ../cli/*.c ; do
   echo $f | sed s+\.\./++ >> POTFILES.in
done   

for f in ../gui/*.c ; do
   if test -f $f ; then
      echo $f | sed s+\.\./++ >> POTFILES.in
   fi
done 

for f in ../gui2/*.c ; do
   echo $f | sed s+\.\./++ >> POTFILES.in
done

for f in ../plugin/*.c ; do
   echo $f | sed s+\.\./++ >> POTFILES.in
done 
