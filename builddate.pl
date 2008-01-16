#!/usr/bin/perl

my $update = 1;
my $test, $line;
my @fields;

my $datestr =`date +%Y-%m-%d`;
chomp ($datestr);

if (open (BUILD_H, "< build.h")) {
    while (<BUILD_H>) {
    	chomp();
	@fields = split(/ +/);
	$test = $fields[4];
	$test =~ s/\\n\"//;
	if ($test eq $datestr) { 
	   print "No need to update build.h\n";
	   $update = 0; 
	}
    }
    close BUILD_H;
}

if ($update) {
    print "Updating build.h\n";
    open (BUILD_H, "> build.h");
    print BUILD_H "#define BUILD_DATE \"build date ", $datestr, "\\n\"\n";
    close BUILD_H;
}
