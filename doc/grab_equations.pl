#!/usr/bin/perl

use strict;

# configurable options:
my $texopt = "12pt";     
my $texpackage = "\\usepackage{mathtime}";
# end configurable options

my ($line, $eqn, $foo, $figfile);
my $textmp = "./eqntmp";

sub usage
{
    die <<"EndUsage";
usage: grab_equations.pl sgmlfile

grab_equations.pl -- A program for making png images of equations from
                     an SGML source file using dbtexmath mark-up.

EndUsage
}

if (@ARGV == 0) { &usage; }
my $doc = $ARGV[0];

open (DOC, "<$doc") || die "Can't open $doc";

sub printtex {
    open (TEX, ">$textmp.tex") || die "Can't open $textmp";
    print TEX "\\documentclass[$texopt]{article}\n";
    print TEX "$texpackage\n";
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

while ($line = <DOC>) {
    if ($line =~ /\<(informal|inline|)equation\>/) { 
	$eqn = $line;
	while ($line !~ /\<\/(informal|inline|)equation\>/) {
	    $line = <DOC>;
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

close (DOC);





