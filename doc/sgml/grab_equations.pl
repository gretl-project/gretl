#!/usr/bin/perl

use strict;

my ($line, $foo, $figfile);
my $manual = "./manual.xml";
my $textmp = "./eqntmp";

open (MAN, "<$manual") || die "Can't open $manual";

sub printtex {
    open (TEX, ">$textmp.tex") || die "Can't open $textmp";
    print TEX "\\documentclass[12pt]{article}\n";
    print TEX "\\usepackage{mathtime\n}";
    print TEX "\\pagestyle{empty}\n";
    print TEX "\\begin{document}\n";
    print TEX "$line\n";
    print TEX "\\end{document}\n";
    close (TEX);
    system ("latex $textmp");
    system ("dvips -E -o $textmp.eps $textmp");
    system ("convert -density 100x100 $textmp.eps $figfile.png");
    system ("rm -f $textmp.*");
}

while ($line = <MAN>) {
    chomp($line);
    if ($line =~ /\<texmath/) {
	($foo, $figfile) = split(/\"/, $line);
        print "Got figfile: $figfile\n";
        $line = <MAN>;
        chomp($line);
        $line =~ s/\s+//g;
        print "read TeX: $line\n";
	print "creating $figfile.png\n";
	printtex();
    }
}




