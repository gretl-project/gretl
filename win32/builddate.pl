#!/usr/bin/perl

$_=`date +%m/%d/%Y`;
chomp;
print "#define BUILD_DATE \"build date ", $_, "\\n\"\n";
