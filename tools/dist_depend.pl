#!/usr/bin/perl -w

use strict;

my $line;
my $got_depend;
my @newlines;

open (MAK, "<Makefile") || die "Can't open Makefile";
open (MAKIN, "<Makefile.in") || die "Can't open Makefile.in";

while ($line = <MAKIN>) {
    if ($line =~ /^# DO NOT DELETE/) {
	last;
    }
    push (@newlines, $line);
}
close (MAKIN);

$got_depend = 0;
while ($line = <MAK>) {
    if ($line =~ /^# DO NOT DELETE/) {
	$got_depend = 1;
    }
    if ($got_depend) {
	push (@newlines, $line);
    }
}
close (MAK);

if ($got_depend == 0) {
    exit 0;
}

open (MAKOUT, ">Makefile.in") || die "Can't write to Makefile.in";

foreach $line (@newlines) {
    print MAKOUT $line;
}

