#!/usr/bin/perl

use strict;

my $line;
my $manual = "./manual.tex";
my $textmp = "./mantmp.tex";
my $inmath = 0;

open (MAN, "<$manual") || die "Can't read $manual";
open (TMP, ">$textmp") || die "Can't write to $textmp";

sub unescape {
    $line =~ s/\\char92{}/\\/g;
    $line =~ s/\\char94{}/^/g;
    $line =~ s/\\char95{}/_/g;
    $line =~ s/\\{/{/g;
    $line =~ s/\\}/}/g;
    $line =~ s/\\\$/\$/g;
}

while ($line = <MAN>) {
    if ($line =~ /BEGINTEXMATH/) {
        $inmath = 1;
        $line =~ s/BEGINTEXMATH//;
    }
    if ($line =~ /ENDTEXMATH/) {
        $inmath = 0;
        $line =~ s/ENDTEXMATH//;
    }    
    if ($inmath) {
	unescape();
    }
    print TMP "$line";
}

close (MAN);
close (TMP);
system("cp $textmp $manual");
system("rm -f $textmp");






