#!/usr/bin/perl

use strict;

my($gotfunc, $header, $cfile, $func, $line, $used);
my(@fields, @funclist, @used);

open (ULIST, ">unused_functions") || die "Can't open output file";

foreach $header (<./*.h>) {

    print "***************\n";
    print "Listing functions in $header...\n";

    open (HDR, "$header") || die "Can't open $header";

    @funclist = ();
    $gotfunc = 0;
    while (<HDR>) {
	if (/functions follow/) { 
	    $gotfunc = 1; 
	    next;
	}
	if (!$gotfunc) {
	    next;
	}
	if (/^[a-zA-Z]/) {
	    @fields = split;
	    $fields[1] =~ s/^\*//;
	    print "   $fields[1]\n";
	    push(@funclist, $fields[1]);
	}
    }
    if (!$gotfunc) {
	print "No functions declared in $header\n";
    }
    close (HDR);

    foreach $func (@funclist) {
	$used = 0;
	print "checking for use of $func...\n";
	foreach $cfile (<../../*/*.c>) {
	    # print "checking $cfile\n";
	    open (SRC, "$cfile") || die "Can't open $cfile";
	    while ($line = <SRC>) {
		if ($line =~ /$func\(/) {
		    print "$cfile uses $func\n";
		    $used = 1;
		    last;
		}
	    }
	    close (SRC);
	    if ($used) {
		last;
	    }
	}
	if (!$used) {
	    print "$func not used in ../../*/*.c\n";
	    print ULIST "$header: $func\n";
	}
    }
    print "****************\n\n";
}
close (ULIST);


