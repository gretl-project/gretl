#!/usr/bin/perl -w

use strict;

my $line;
my $got_depend;
my @newlines;

open (MAKDEP, "<Make.windep") || die "Can't open Make.windep";
open (MAK, "<Makefile") || die "Can't open Makefile";

while ($line = <MAK>) {
    if ($line =~ /^# DO NOT DELETE/) {
	last;
    }
    push (@newlines, $line);
}
close (MAK);

$got_depend = 0;
while ($line = <MAKDEP>) {
    if ($line =~ /^# DO NOT DELETE/) {
	$got_depend = 1;
    }
    if ($got_depend) {
	$line =~ s+^../gui2/++;
	$line =~ s+^../cli/++;
	$line =~ s+^../plugin/++;
	push (@newlines, $line);
    }
}
close (MAKDEP);

if ($got_depend == 0) {
    exit 0;
}

open (MAKOUT, ">Makefile") || die "Can't write to Makefile";

foreach $line (@newlines) {
    print MAKOUT $line;
}

