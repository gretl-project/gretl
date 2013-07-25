#!/bin/sh

# Shell script to update an installation of gretl on win32, by
# copying files from a build tree.

######### CONFIGURE HERE #########
PREFIX="c:/progra~1/gretl"
TOPDIR=".."
WINFILES=/mingw/winbuild
HLPDIR=${WINFILES}/help
DOCDIR=${WINFILES}/doc
install="cp -a"
##################################

# put binaries in place
${install} gretlcli.exe $PREFIX
${install} gretl.exe $PREFIX

# gretl dynamic libs
for f in dlls/*.dll ; do
    ${install} $f $PREFIX
done

# gretl plugins
for f in plugins/*.dll ; do
    ${install} $f $PREFIX/plugins
done

# auxiliary dlls
for f in \
    libgmp-3.dll libxml2.dll readline5.dll history5.dll zlib1.dll \
    libblas.dll liblapack.dll libgtksourceview.dll \
    iconv.dll intl.dll libfftw3-3.dll libpng12-0.dll
do
    ${install} ${WINFILES}/misc-dll/$f $PREFIX
done

# XML UI-definition files
for f in ${TOPDIR}/gui2/*.xml ; do
  ${install} $f $PREFIX/ui
done

# gretl lang files for gtksourceview
SPECDIR=$PREFIX/share/gtksourceview-1.0/language-specs
mkdir -p $SPECDIR
${install} gretl.lang $SPECDIR
${install} ${TOPDIR}/gui2/gnuplot.lang $SPECDIR
${install} ${TOPDIR}/gui2/R.lang $SPECDIR

# help files, license, logo
${install} ${HLPDIR}/gretlgui.hlp $PREFIX/gretlgui_hlp.txt
${install} ${HLPDIR}/gretlcmd.hlp $PREFIX/gretlcmd_hlp.txt
${install} ${HLPDIR}/gretlcli.hlp $PREFIX/gretlcli_hlp.txt
${install} ${HLPDIR}/genrcli.hlp $PREFIX/genrcli.hlp
${install} ${HLPDIR}/genrgui.hlp $PREFIX/genrgui.hlp
${install} ${HLPDIR}/genrcli.hlp.it $PREFIX/genrcli.hlp.it
${install} ${HLPDIR}/genrgui.hlp.it $PREFIX/genrgui.hlp.it
for f in ${TOPDIR}/share/texfigs/*.png ; do
   ${install} $f $PREFIX/helpfigs
done
for lang in es it ; do
   ${install} ${HLPDIR}/gretlgui.hlp.${lang} \
       $PREFIX/gretlgui_hlp_${lang}.txt
   ${install} ${HLPDIR}/gretlcmd.hlp.${lang} \
       $PREFIX/gretlcmd_hlp_${lang}.txt
   ${install} ${HLPDIR}/gretlcli.hlp.${lang} \
       $PREFIX/gretlcli_hlp_${lang}.txt
done
${install} ${TOPDIR}/COPYING $PREFIX
${install} ${TOPDIR}/pixmaps/gretl-logo.xpm $PREFIX
${install} ${DOCDIR}/gretl-guide.pdf $PREFIX/doc
${install} ${DOCDIR}/gretl-ref.pdf $PREFIX/doc

# Ramanathan data files
DATADIR=${TOPDIR}/share/data
for f in ${DATADIR}/*.gdt ; do
    ${install} $f $PREFIX/data
done
for f in ${DATADIR}/*.dtd ; do
    ${install} $f $PREFIX/data
done
${install} ${DATADIR}/descriptions $PREFIX/data

# Greene data files
DATADIR=${TOPDIR}/share/data/greene
for f in ${DATADIR}/*.gdt ; do
    ${install} $f $PREFIX/data/greene
done
${install} ${DATADIR}/wg_descriptions $PREFIX/data/greene

# Gretl data files
DATADIR=${TOPDIR}/share/data/misc
for f in ${DATADIR}/*.gdt ; do
    ${install} $f $PREFIX/data/misc
done
${install} ${DATADIR}/descriptions $PREFIX/data/misc

# NIST test data files
for f in ${TOPDIR}/tests/*.dat ; do
    ${install} $f $PREFIX/data/nist
done

# top-level script files
for f in ${TOPDIR}/share/scripts/*.inp ; do
    ${install} $f $PREFIX/scripts
done
${install} ${TOPDIR}/share/scripts/ps_descriptions $PREFIX/scripts
${install} ${TOPDIR}/share/scripts/wg_ps_descriptions $PREFIX/scripts

# "misc" script files
for f in ${TOPDIR}/share/scripts/misc/*.inp ; do
    ${install} $f $PREFIX/scripts/misc
done
${install} ${TOPDIR}/share/scripts/misc/ps_descriptions $PREFIX/scripts/misc

# database files
make -C db
${install} db/fedstl.bin $PREFIX/db
${install} ${TOPDIR}/share/bcih/fedstl.idx $PREFIX/db

# gretl translations (make sure they're up to date first)
make -C mo
LANGS=`cat ${TOPDIR}/po/LINGUAS | grep -v ^#`
for lang in $LANGS ; do
   mkdir -p $PREFIX/locale/${lang}/LC_MESSAGES
   ${install} mo/$lang.mo $PREFIX/locale/$lang/LC_MESSAGES/gretl.mo
done

# misc files
${install} ${TOPDIR}/plugin/data/urcdata.gz $PREFIX/plugins/data
${install} ${TOPDIR}/plugin/data/dwdata.gz $PREFIX/plugins/data
${install} gretl_website.url $PREFIX/gretl_website.url
${install} updater/gretl_updater.exe $PREFIX
date > gretl.stamp
${install} gretl.stamp $PREFIX


