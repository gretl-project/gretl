#!/usr/bin/perl

use strict;

my ($line, $texline, $foo, $figfile);
my $manual = "./manual.xml";
my $textmp = "./eqntmp";

open (MAN, "<$manual") || die "Can't open $manual";

sub printtex {
    open (TEX, ">$textmp.tex") || die "Can't open $textmp";
    print TEX "\\documentclass[12pt]{article}\n";
    print TEX "\\usepackage{mathtime\n}";
    print TEX "\\pagestyle{empty}\n";
    print TEX "\\begin{document}\n";
    print TEX "$texline\n";
    print TEX "\\end{document}\n";
    close (TEX);
    system ("latex $textmp");
    system ("dvips -E -o $textmp.eps $textmp");
    system ("convert -density 100x100 $textmp.eps $figfile.png");
    system ("rm -f $textmp.*");
}

sub grabtex {
    my $line = <MAN>; # grab line with literal tex
    $line =~ s/\s+/ /g;
    print "read TeX: $line";
    $texline += $line;
}

while ($line = <MAN>) {
    chomp($line);
    if ($line =~ /equation\>/) {
        $line = <MAN>;
	($foo, $figfile) = split(/\"/, $line);
        print "Got figfile: $figfile\n";
	$line = <MAN>; # skip <texmath>
	grabtex();
	while ($texline !~ /\$ +$/) {
	    grabtex();
	}
	print "creating $figfile.png\n";
	printtex();
    }
}




