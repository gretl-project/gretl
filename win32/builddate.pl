#!/usr/bin/perl

my $update = 1;
my $test, $line;
my @fields;

my $datestr =`date +%m/%d/%Y`;
chomp ($datestr);

if (open (BUILD_H, "< build.h")) {
    while (<BUILD_H>) {
	@fields = split(/ +/);
	$test = $fields[4];
	$test =~ s/\\n\"//;
	if ($test eq $datestr) { $update = 0; }
    }
    close BUILD_H;
}

if ($update) {
    open (BUILD_H, "> build.h");
    print BUILD_H "#define BUILD_DATE \"build date ", $datestr, "\\n\"\n";
    close BUILD_H;
}
