#!/bin/sh

DATADIR="http://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA"

# Data file example:
# http://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Misra1a.dat

for datafile in `cat datalist` ; do
  if [ -f ./$datafile ] ; then
     echo "$datafile already present"
  else 
    echo "getting $datafile"
    wget $DATADIR/$datafile
  fi
done
  
