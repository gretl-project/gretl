#!/usr/bin/perl

# script to write gretl database using data from a given chapter
# (01-16) of the NBER Macrohistory Database.
# See http://www.nber.org/databases/macrohistory/contents/index.html
# Requires that get_nber_data.sh be run first.

use strict;

my($dbfile, $idxtitle);
my($line, $linebak, $linebak2, $sername, $descrip, $area, $units, $i);
my($pd, $pdstr, $stobs, $endobs, $obs, $nobs, $n_2, $ndiff, $sa, $yr, $per);
my(@fields, @obss);
my @chapters = (
"Production of Commodities",
"Construction",
"Transportation and Public Utilities",
"Prices",
"Stocks of Commodities",
"Distribution of Commodities",
"Foreign Trade",
"Income and Employment",
"Financial Status of Business",
"Savings and Investment",
"Security Markets",
"Volume of Transactions",
"Interest Rates",
"Money and Banking",
"Government and Finance",
"Leading Indicators"
);

my $id = $ARGV[0];

# takes args: endobs, frequency, stobs
sub dateton {
    my ($endobs, $pd, $stobs) = @_;
    my ($maj, $min, $start_maj, $start_min, $n);

    if ($endobs =~ /\./) {
      ($maj, $min) = split ('\.', $endobs);
    } else {
      $maj = $endobs;
      $min = 0;
    }
    if ($stobs =~ /\./) {
      ($start_maj, $start_min) = split ('\.', $stobs);
    } else {
      $start_maj = $stobs;
      $start_min = 0;
    }
    $n = $pd * ($maj - $start_maj);
    $n += $min - $start_min;
   
    return $n;
}

sub top_n_tail {
    my ($str) = @_;
    $str =~ s/^\s+//;
    $str =~ s/^\"\s+//;
    $str =~ s/\s+$//;
    $str =~ s/\s+\"$//;
    return $str;
}

sub do_missval {
    if ($obs =~ /0.10*E-/ || $obs =~ /1.e-37/ ||
	$obs eq "M" || $obs eq "m") {
	$obs = -999.0;
    }    
}

sub is_missing {
    my ($obs) = @_;
    if ($obs == -999.0) { 
	return 1; 
    }
    return 0;
}

sub decrement {
    my ($endobs, $pd, $diff) = @_;
    my ($yr, $per) = split(/\./, $endobs);
    if ($pd == 1) {
	$yr -= $diff;
	$endobs = $yr;
    } else {
	$per -= $diff;
	while ($per <= 0) {
	    $yr -= 1;
	    $per += $pd;
	}
	$endobs = $yr . ".$per";
	if ($pd == 4) {
	    $endobs = do_qtrs($endobs);
	}
	if ($pd == 12) {
	    $endobs = do_months($endobs);
	}
    }
    return $endobs;
}

sub do_qtrs {
    my ($obs) = @_;
    my ($yr, $per) = split('\.', $obs);
    if ($per =~ s/^0//) {
	$obs = $yr . ".$per";
    }
    return $obs;
}

sub do_months {
    my ($obs) = @_;
    my ($yr, $per) = split('\.', $obs);
    if ($per < 10 && $per =~ /^[1-9]/) {
	$per = "0" . $per;
	$obs = $yr . ".$per";
    }
    return $obs;
}

# main program starts

$id > 16 && die "Chapter $id is out of range!\n";

if (($id < 10) && !($id =~ /^0/)) {
    $id = "0$id";
}

$idxtitle = "# NBER Macrohistory Data -- $chapters[$id - 1]";

unlink "nberhist$id.bin";
unlink "nberhist$id.idx";
open (BIN, ">>nberhist$id.bin") || die "Can't open nberhist$id.bin";
open (IDX, ">>nberhist$id.idx") || die "Can't open nberhist$id.idx";

print IDX "$idxtitle\n";

foreach $dbfile (<www.nber.org/$id/*.db>) {

    open (TSP, "$dbfile") || die "Can't open tsp source $dbfile";

    $sername = $dbfile;
    $sername =~ s/\.db//;
    @fields = split('/', $sername);
    $sername = $fields[-1];
    $descrip = "";
    $sa = "";
    $pd = "";
    $stobs = "";
    $endobs = "";
    $area = "";
    $units = "";

    # deal with the comments first
    while ($line = <TSP>) {
        if ($line =~ /^\"/) {
	    $line =~ s/^\"c//;
	    $line =~ s/^\"//;
	    $line =~ s/^\s+//;
	    $line =~ s/\s+$//;
	    $line =~ s/\"$//;
	    $line =~ s/\s+$//;
	    # print "$line\n";
	} else {
	    last;
	}
	if ($line =~ s/AREA COVERED:\s+//) {
	    $area = $line;
	}
	if ($line =~ /SEASONAL ADJUSTMENT:  NONE/) {
	    $sa = "NSA";
	}
	if ($line =~ /SEASONALLY ADJUSTED/) {
	    $sa = "SA";
	}
	if ($line =~ s/UNITS:\s+//) {
	    $units = $line;
	}
	if ($line =~ /-------/) {
	    $descrip = "$linebak2 $linebak";
	    $descrip =~ s/, SEASONALLY ADJUSTED//;
	    $descrip = top_n_tail($descrip);
	}
	$linebak2 = $linebak;
	$linebak = $line;
    }
    if ($descrip =~ /^CHECKED/ || ! $descrip =~ /^\w/ ||
	$descrip =~ /^\./) {
      print("$sername: missing definition, skipping...\n");
      close (TSP);
      next;
    }

    # next deal with sample information
    $pd = $line;
    chomp($pd);
    $pd = -1 * $pd;
    $stobs = <TSP>;
    $stobs =~ s/\s+$//;
    $endobs = <TSP>;
    $endobs =~ s/\s+$//;
    if ($pd == 1) {
        $pdstr = "A";
	$stobs =~ s/\.//;
	$endobs =~ s/\.//;
    }
    elsif ($pd == 4) { 
        $pdstr = "Q"; 
	$stobs = do_qtrs($stobs);
	$endobs = do_qtrs($endobs);
    }
    elsif ($pd == 12) { 
        $pdstr = "M"; 
	$stobs = do_months($stobs);
	$endobs = do_months($endobs);
    }
    else {
      print "$sername: unrecognized frequency: $pd: skipping...\n";
      close (TSP);
      next;
    }
    $nobs = 0;
    $n_2 = 0;
    @obss = ();
    while ($obs = <TSP>) {
        $nobs++;
	chomp($obs);
	do_missval();
	# print "$obs\n";
	push(@obss, $obs);
	# print BIN pack("f", $obs);
    }
    close (TSP);

    $n_2 = dateton($endobs, $pd, $stobs) + 1;
    $ndiff = $nobs - $n_2;
    if ($ndiff < 0) {
	print "$sername: n counted ($nobs) < n computed ($n_2): skipping\n";
    }
    elsif ($ndiff > 0) {
	print "Warning: $sername: n (counted) = $nobs; n (computed) = $n_2\n";
	$i = $nobs - 1;
	while (is_missing($obss[$i])) {
	    if ($ndiff > 0) {
		$ndiff -= 1;
	    }
	    if ($ndiff == 0) {
		last;
	    }
	    $i -= 1;
	}
	if ($ndiff == 0) {
	    print "But this is due to missing values at end, OK\n";
	    $endobs = decrement($endobs, $pd, $nobs - $n_2);
	    $nobs = $n_2;
	} else {
	    print "Can't fix this, skipping...\n";
	    next;
	}
    }
    if ($ndiff == 0) {
	for ($i = 0; $i < $nobs; $i++) {
	    print BIN pack("f", $obss[$i]);
	}
	# print to index file
	print IDX "$sername $descrip, $area, $units $sa\n";
	print IDX "$pdstr  $stobs - $endobs  n = $nobs\n";
    }

} # end of one dbfile

close (BIN);
close (IDX);
