eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q'
  if 0;
use strict;
$^W=1; # turn warning on

my @figfiles;
my ($line, $eqn, $foo, $figfile);
my ($i, $n, $match);
my ($opt, $den, $pkg);
my $textmp = "./eqntmp";

sub usage
{
    die <<"EndUsage";
usage: texmath2png.pl sgmlfile

texmath2png.pl -- A program for making png images of equations from
                  an SGML source file using dbtexmath mark-up.

EndUsage
}

sub unescape {
    $eqn =~ s/&#38;/&/g;
    $eqn =~ s/&#62;/\>/g;
    $eqn =~ s/&#60;/\</g;
}

sub printtex {
    open (TEX, ">$textmp.tex") || die "Can't open $textmp";
    if ($opt) {
	print TEX "\\documentclass[$opt]{article}\n";
    } else {
	print TEX "\\documentclass{article}\n";
    }	
    if ($pkg) {
        print TEX "\\usepackage{$pkg}\n";
    }
    print TEX "\\pagestyle{empty}\n";
    print TEX "\\begin{document}\n";
    print TEX "$eqn";
    print TEX "\\end{document}\n";
    close (TEX);
    system ("latex $textmp");
    system ("dvips -o $textmp.eps $textmp -E");
    system ("convert -density $den $textmp.eps $figfile");
    system ("rm -f $textmp.*");
}
sub get_latexopt {
    $opt = $line;
    $line = <DOC>;
    $opt = $opt . $line;
    if ($opt =~ /\<latexopt\s*\>\s*((?:.|\s)*)<\/latexopt/) {
	$opt = $1;
	print "Got LaTeX option: '" . $opt . "'\n";
    }
}

sub get_density {
    $den = $line;
    $line = <DOC>;
    $den = $den . $line;
    if ($den =~ /\<density\s*\>\s*((?:.|\s)*)<\/density/) {
	$den = $1;
	print "Got image density spec: '" . $den . "'\n";
    }
}

sub get_package {
    $pkg = $line;
    $line = <DOC>;
    $pkg = $pkg . $line;
    if ($pkg =~ /\<usepackage\s*\>\s*((?:.|\s)*)<\/usepackage/) {
	$pkg = $1;
	print "Got package call: '" . $pkg . "'\n";
    }
}

if (@ARGV == 0) { &usage; }
my $doc = $ARGV[0];

open (DOC, "<$doc") || die "Can't open $doc";

$den = "96x96";
$opt = undef;
$pkg = undef;

while ($line = <DOC>) {
    if ($line =~ /\<latexopt/) { get_latexopt(); }
    if ($line =~ /\<density/) { get_density(); }
    if ($line =~ /\<usepackage/) { get_package(); }
    if ($line =~ /\<texequation/) { 
	$eqn = $line;
	while ($line !~ /\<\/texequation/) {
	    $line = <DOC>;
	    $eqn = $eqn . $line;
	}
	if ($eqn =~ s/\s*fileref="([a-zA-Z0-9_\/\.]+)"//) { 
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
	if ($eqn =~ /\<texequation\s*\>\s*((?:.|\s)*)<\/texequation/) {
	    $eqn = $1;
	    unescape();
	    print "got tex math $eqn\n";
	}
	printtex();
    }
}

close (DOC);





