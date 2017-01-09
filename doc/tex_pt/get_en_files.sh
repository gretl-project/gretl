#!/usr/bin/sh

# Use this script to link untranslated files to original English ones.
# To start translating delete the linked file and hard copy it.
# Remove from list any file you started to work on.

untranslated="
appendices.tex 
cheatsheet.tex 
criteria.tex 
datafiles.tex 
datatypes.tex 
df.tex 
discrete.tex 
dpanel.tex 
forecast.tex 
functions.tex 
genr.tex 
gmm.tex 
graphs.tex 
gretlOctave.tex 
gretlOx.tex 
gretlR.tex 
gretltex.tex 
introduction.tex 
kalman.tex 
looping.tex 
matrices.tex 
mle.tex 
modes.tex 
nls.tex 
odbc.tex 
panel.tex 
persistent.tex 
probit.tex 
quantreg.tex 
robust_vcv.tex 
sampling.tex 
starting.tex 
timeseries.tex 
trouble.tex 
var.tex 
vecm.tex
"

for i in $untranslated; do
    if [ `diff ../tex/$i $i 2>/dev/null; echo $?` -eq 2 ];
    then
       ln -sf ../tex/$i $i
    fi
done

