#!/bin/sh

for f in *.html ; do
  if grep figures_it $f >/dev/null 2>&1 ; then
     sed s/figures_it/figures/g < $f > tmp.html
     mv tmp.html $f
  fi
done 
