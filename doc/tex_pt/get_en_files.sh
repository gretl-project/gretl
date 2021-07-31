#!/usr/bin/sh

# Use this script to copy untranslated files from original English ones.
# Remove from the "untranslated" list any file you started to work on.

untranslated="
calendar.tex
complex_matrices.tex
cmdtopics.tex
funcbody.tex
gretl-guide-a4.tex
gretl-guide.tex
gretlhds.sty
gretlJulia.tex
gretl-keys-a4.tex
gretl-keys.tex
gretl-lite.sty
gretl-mpi-a4.tex
gretl-mpi.tex
gretlPython.tex
gretl-ref-a4.tex
gretl-ref.tex
gretlStata.tex
gretl-svm.tex
hansl-primer-a4.tex
hp-bundles.tex
hp-ctrlflow.tex
hp-dataset.tex
hp-estimate.tex
hp-functions.tex
hp-further.tex
hp-greeks.tex
hp-hello.tex
hp-matrices.tex
hp-numerical.tex
hp-output.tex
hp-reference.tex
hp-series.tex
join.tex
lucidabr.sty
midas-gretl.tex
midas.tex
mixfreq.tex
nonparam.tex
numerical.tex
optshort.tex
pkgbook-a4.tex
pkgbook.tex
realtime.tex
refbody.tex
string_series.tex
system.tex
tsfilter.tex
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
tdisagg.tex
timeseries.tex 
trouble.tex 
var.tex 
vecm.tex
"

for i in $untranslated; do
    a=`diff -q ../tex/$i ./$i 2>/dev/null`
    code=`echo $?`
    if [ $code -eq 2 ];
    then
       cp ../tex/$i $i
    fi
done



