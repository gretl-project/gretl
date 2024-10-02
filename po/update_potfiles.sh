#!/bin/sh

cat /dev/null > POTFILES.in

for f in ../lib/src/*.c ; do
   echo $f | sed s+\.\./++ >> POTFILES.in
done

for f in ../cli/*.c ; do
   echo $f | sed s+\.\./++ >> POTFILES.in
done

for f in ../gui/*.c ; do
   echo $f | sed s+\.\./++ >> POTFILES.in
done

for f in ../plugin/*.c ; do
   echo $f | sed s+\.\./++ >> POTFILES.in
done

for f in ../doc/commands/*.c ; do
   echo $f | sed s+\.\./++ >> POTFILES.in
done

# we should ignore the other .c files in addons
echo "addons/addons-i18n.c" >> POTFILES.in
