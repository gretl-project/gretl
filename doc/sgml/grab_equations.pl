#!/usr/bin/perl

use strict;

my ($line, $texline, $foo, $figfile);
my $manual = "./manual.xml";
my $textmp = "./eqntmp";
my $i = 0;

open (MAN, "<$manual") || die "Can't open $manual";

sub printtex {
    open (TEX, ">$textmp.tex") || die "Can't open $textmp";
    print TEX "\\documentclass[12pt]{article}\n";
    print TEX "\\usepackage{mathtime}\n";
    print TEX "\\pagestyle{empty}\n";
    print TEX "\\begin{document}\n";
    print TEX "$texline";
    print TEX "\\end{document}\n";
    close (TEX);
    # print "printed tex file $i\n"; $i++;
    system ("latex $textmp");
    system ("dvips -o $textmp.eps $textmp -E");
    system ("convert -density 100x100 $textmp.eps $figfile.png");
    system ("cp $textmp.tex $textmp$i.tex"); $i++;
    system ("rm -f $textmp.*");
}

while ($line = <MAN>) {
    if ($line =~ /\<informalequation\>/) {
        $texline = '';
        $line = <MAN>;
	($foo, $figfile) = split(/\"/, $line);
        print "Got figfile: $figfile\n";
	$line = <MAN>; # skip opening <texmath>
	$line = <MAN>;
	while ($line !~ m+\</texmath\>+) {
	    print "Got tex line: $line";
	    $texline = $texline . $line;
	    $line = <MAN>;
	}
	print "creating $figfile.png\n";
	$texline =~ s/^\S+//g;
	print "Here's texline:\n";
	print "$texline";
	printtex();
    }
}





