eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q'
  if 0;
use strict;
$^W=1; # turn warning on

my @figfiles;
my ($line, $eqn, $figfile);
my ($i, $n, $match);
my ($opt, $den, $pkg);
my $den_dflt = "96x96";
my $textmp = "./eqntmp";
my $figdir = ".";

sub usage
{
    die <<"EndUsage";
usage: texmath2png.pl xmlfile

texmath2png.pl -- A program for making png images of equations from
                  an XML source file using dbtexmath mark-up.

EndUsage
}

sub add_png_ext {
    if ($figfile !~ /\.png$/) {
        $figfile = $figfile . ".png"
    }
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
    system ("convert -density $den $textmp.eps $figdir/$figfile");
    system ("rm -f $textmp.*");
}

if (@ARGV == 0) { &usage; }
my $doc = $ARGV[0];
if (@ARGV > 1) {
  $figdir = $ARGV[1];
}

open (DOC, "<$doc") || die "Can't open $doc";

$den = $den_dflt;
$opt = undef;
$pkg = undef;

while ($line = <DOC>) {
    if ($line =~ /latexopt="(\S+)"/) { 
	$opt = $1; 
	print "got LaTeX document option $opt\n";
    }
    if ($line =~ /density="(\S+)"/) { 
	$den = $1;
	print "got density specification $den\n";
    }
    if ($line =~ /usepackage="(\S+)"/) { 
	$pkg = $1; 
	print "got usepackage line $pkg\n";
    }
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
        add_png_ext();
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





