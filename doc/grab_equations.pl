#!/usr/bin/perl

use strict;

# configurable options:
my $texopt = "12pt";     
my $texpackage = "\\usepackage{mathtime}";
# end configurable options

my @figfiles;
my ($line, $eqn, $foo, $figfile);
my ($i, $n, $match);
my $textmp = "./eqntmp";

sub usage
{
    die <<"EndUsage";
usage: grab_equations.pl sgmlfile

grab_equations.pl -- A program for making png images of equations from
                     an SGML source file using dbtexmath mark-up.

EndUsage
}

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

if (@ARGV == 0) { &usage; }
my $doc = $ARGV[0];

open (DOC, "<$doc") || die "Can't open $doc";

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
	$n = @figfiles;
	$match = 0;
	for ($i = 0; $i < $n && !$match; $i++) {
	    if ($figfiles[$i] eq $figfile) { 
		$match = 1;
	    }
	}
	if ($match) {
	    print "$figfile: already done, skipping\n";
	    next; 
	}
	push(@figfiles, $figfile);
	if ($eqn =~ /\<alt\>((?:.|\s)*)\<\/alt\>/) {
	    print "got texmath $1\n";
	    $eqn = $1;
	}
	printtex();
    }
}

close (DOC);





