#!/usr/bin/perl

# script to abbreviate NBER macrohistory data descriptions

my $line;

use strict;

while ($line = <STDIN>) {

    $line =~ s/MILLIONS OF DOLLARS/\$MILL./;
    $line =~ s/BILLIONS OF DOLLARS/\$BN./;
    $line =~ s/ AND / & /;
    $line =~ s/UNITED STATES/U.S./;
    $line =~ s/INCLUDING/INC./;
    $line =~ s/EXCLUDING/EXCL./;
    $line =~ s/MANUFACTURING/MFG./;
    $line =~ s/PERCENTAGE/%/;
    $line =~ s/PERCENT/%/;
    $line =~ s/PER CENT/%/;
    $line =~ s/FEDERAL RESERVE/FED./;
    $line =~ s/DEPARTMENT/DEPT./;
    $line =~ s/NUMBER/NO./;
    $line =~ s/ IN / /;
    $line =~ s/MANUFACTURERS/MFRS./;
    $line =~ s/ PLUS / + /;
    $line =~ s/MILLIONS OF/MILL./;
    $line =~ s/AVERAGE/AVG/;
    $line =~ s/WEIGHTED/WTD/;
    $line =~ s/THOUSANDS/000s/;

    print "$line";
}
