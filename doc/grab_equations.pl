#!/usr/bin/perl

use strict;

my ($line, $eqn, $foo, $figfile);
my $manual = "./manual.sgml";
my $textmp = "./eqntmp";

open (MAN, "<$manual") || die "Can't open $manual";

sub printtex {
    open (TEX, ">$textmp.tex") || die "Can't open $textmp";
    print TEX "\\documentclass[12pt]{article}\n";
    print TEX "\\usepackage{mathtime}\n";
    print TEX "\\pagestyle{empty}\n";
    print TEX "\\begin{document}\n";
    print TEX "$eqn";
    print TEX "\\end{document}\n";
    close (TEX);
    system ("latex $textmp");
    system ("dvips -o $textmp.eps $textmp -E");
    system ("convert -density 100x100 $textmp.eps $figfile");
    system ("rm -f $textmp.*");
}

while ($line = <MAN>) {
    if ($line =~ /\<(informal|inline|)equation\>/) { 
	$eqn = $line;
	while ($line !~ /\<\/(informal|inline|)equation\>/) {
	    $line = <MAN>;
	    $eqn = $eqn . $line;
	}
	if ($eqn =~ /fileref="([a-zA-Z0-9_\/\.]+)"/) { 
	    print "got fileref $1\n";
	    $figfile = $1;
	}
	if ($eqn =~ /\<texmath\>((?:.|\s)*)\<\/texmath\>/) {
	    print "got texmath $1\n";
	    $eqn = $1;
	}
	printtex();
    }
}





