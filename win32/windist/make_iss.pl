#!/usr/bin/perl

# convert MANIFEST to .iss script for making Win32 installer

use strict;

my $zname;
my $ldash;
my $rdash;
my $verstr;
my $line; 
my $subdir;
my @pathbits;
my $i;

sub usage
{
    die <<"EndUsage";
usage: make_iss gretl_version_number

EndUsage
}


if (@ARGV == 0) { &usage; }
$verstr = $ARGV[0];

print STDERR "Making installer script for gretl version $verstr...\n";

print "; -- gretl.iss --\n";
print "[Setup]\n";
print "AppName=gretl\n";
print "AppVerName=gretl version $verstr\n";
print "AppCopyright=Copyright C 1999-2002 Allin Cottrell\n";
print "DefaultDirName=c:\\userdata\\gretl\n";
print "DefaultGroupName=gretl\n";
print "UninstallDisplayIcon={app}\\gretlw32.exe\n";
print "ChangesAssociations=yes\n";
print "MinVersion=4,3.51\n";
# print "UsePreviousAppDir=no\n";

print "\n[Dirs]\n";
print "Name: \"{app}\\user\"\n";
print "Name: \"{app}\\data\"\n";
print "Name: \"{app}\\data\\greene\"\n";
print "Name: \"{app}\\scripts\"\n";
print "Name: \"{app}\\db\"\n";
print "Name: \"{app}\\plugins\"\n";
print "Name: \"{app}\\nist\"\n";
print "Name: \"{app}\\locale\"\n";
print "Name: \"{app}\\locale\\es\"\n";
print "Name: \"{app}\\locale\\es\\LC_MESSAGES\"\n";
print "Name: \"{app}\\..\\gnuplot\"\n";
print "Name: \"{app}\\..\\gnuplot\\demo\"\n";

# extra GTK module dirs for gtk-2.0
print "Name: \"{app}\\lib\"\n";
print "Name: \"{app}\\lib\\gtk-2.0\"\n";
print "Name: \"{app}\\lib\\gtk-2.0\"\n";
print "Name: \"{app}\\lib\\gtk-2.0\\2.2.0\"\n";
print "Name: \"{app}\\lib\\gtk-2.0\\2.2.0\\loaders\"\n";
print "Name: \"{app}\\lib\\pango\"\n";
print "Name: \"{app}\\lib\\pango\\1.0.0\"\n";
print "Name: \"{app}\\lib\\pango\\1.0.0\\modules\"\n";

# GTK message catalog
print "Name: \"{app}\\lib\\locale\"\n";
print "Name: \"{app}\\lib\\locale\\es\"\n";
print "Name: \"{app}\\lib\\locale\\es\\LC_MESSAGES\"\n";
print "Name: \"{app}\\lib\\locale\\fr\"\n";
print "Name: \"{app}\\lib\\locale\\fr\\LC_MESSAGES\"\n";

# module catalogs
print "Name: \"{app}\\etc\"\n";
print "Name: \"{app}\\etc\\gtk-2.0\"\n";
print "Name: \"{app}\\etc\\pango\"\n";

print "\n[Files]\n";

while ($line = <STDIN>) {
   $line =~ s+/+\\+g;
   chomp($line);
   @pathbits = split(/\\/, $line);
   print "Source: \"$line\"; "; 
   if ($line =~ /gnuplot/) {
      print "Destdir: \"{app}\\..";         
      for ($i = 0; $i < @pathbits - 1; $i++) {
         print "\\$pathbits[$i]";
      }
   } else {
      print "Destdir: \"{app}";
      for ($i = 1; $i < @pathbits - 1; $i++) {
         print "\\$pathbits[$i]";
      }
   }     
   print "\"\n";
}

print "\n[Icons]\n";
print "Name: \"{group}\\gretl\"; Filename: \"{app}\\gretlw32.exe\"\n";
print "Name: \"{group}\\Gretl Web Site\"; Filename: \"{app}\\gretl_website.url\"\n";
print "Name: \"{group}\\gretl updater\"; Filename: \"{app}\\gretl_updater.exe\"\n";

print "\n[Registry]\n";
print "; Start \"gretl\" keys under HKEY_CURRENT_USER and HKEY_CLASSES_ROOT.\n"; 
print "Root: HKCR; Subkey: \"Software\\gretl\"; Flags: uninsdeletekey\n"; 
print "Root: HKCU; Subkey: \"Software\\gretl\"; Flags: uninsdeletekey\n";
print "Root: HKCR; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"gretldir\"; ValueData: \"{app}\"\n"; 
print "Root: HKCR; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"gnuplot\"; ValueData: \"{app}\\..\\gnuplot\\wgnuplot.exe\"\n"; 
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"userdir\"; ValueData: \"{app}\\user\\\"\n";   
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"expert\"; ValueData: \"false\"\n";    
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"binbase\"; ValueData: \"{app}\\db\\\"\n";   
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"ratsbase\"; ValueData: \"f:\\\"\n";   
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"dbhost_ip\"; ValueData: \"152.17.150.2\"\n";   
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"fontspec\"; ValueData: \"Courier New 10\"\n";
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"Rcommand\"; ValueData: \"RGui.exe\"\n";
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"updater\"; ValueData: \"false\"\n";
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"calculator\"; ValueData: \"calc.exe\"\n";   
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"editor\"; ValueData: \"winword.exe\"\n";   
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"toolbar\"; ValueData: \"true\"\n";   

# Establish file associations

# first for data files...
print "Root: HKCR; Subkey: \".gdt\"; ValueType: string; ValueName: ";
print "\"\"; ValueData: \"GretlDataFile\"; Flags: uninsdeletevalue\n"; 
print "Root: HKCR; Subkey: \".gdt\"; ValueType: string; ";
print "ValueName: \"Content Type\"; ValueData: \"application/x-gretldata\"\n";
print "Root: HKCR; Subkey: \"GretlDataFile\"; ValueType: string; ValueName: ";
print "\"\"; ValueData: \"Gretl data file\"; Flags: uninsdeletekey\n";
print "Root: HKCR; Subkey: \"GretlDataFile\\DefaultIcon\"; ValueType: string; ";
print "ValueName: \"\"; ValueData: \"{app}\\gretlw32.exe,1\"\n"; 
print "Root: HKCR; Subkey: \"GretlDataFile\\shell\\open\\command\"; ";
print "ValueType: string; ValueName: \"\"; ValueData: ";
print "\"\"\"{app}\\gretlw32.exe\"\" \"\"%1\"\"\"\n";

# then for session files
print "Root: HKCR; Subkey: \".gretl\"; ValueType: string; ValueName: ";
print "\"\"; ValueData: \"GretlSessionFile\"; Flags: uninsdeletevalue\n"; 
print "Root: HKCR; Subkey: \".gretl\"; ValueType: string; ";
print "ValueName: \"Content Type\"; ValueData: \"application/x-gretlsession\"\n"; 
print "Root: HKCR; Subkey: \"GretlSessionFile\"; ValueType: string; ValueName: ";
print "\"\"; ValueData: \"Gretl session file\"; Flags: uninsdeletekey\n";
print "Root: HKCR; Subkey: \"GretlSessionFile\\DefaultIcon\"; ValueType: string; ";
print "ValueName: \"\"; ValueData: \"{app}\\gretlw32.exe,2\"\n"; 
print "Root: HKCR; Subkey: \"GretlSessionFile\\shell\\open\\command\"; ";
print "ValueType: string; ValueName: \"\"; ValueData: ";
print "\"\"\"{app}\\gretlw32.exe\"\" -r \"\"%1\"\"\"\n";

# and for simple script files
print "Root: HKCR; Subkey: \".inp\"; ValueType: string; ValueName: ";
print "\"\"; ValueData: \"GretlScriptFile\"; Flags: uninsdeletevalue\n"; 
print "Root: HKCR; Subkey: \".inp\"; ValueType: string; ";
print "ValueName: \"Content Type\"; ValueData: \"application/x-gretlscript\"\n"; 
print "Root: HKCR; Subkey: \"GretlScriptFile\"; ValueType: string; ValueName: ";
print "\"\"; ValueData: \"Gretl script file\"; Flags: uninsdeletekey\n";
print "Root: HKCR; Subkey: \"GretlScriptFile\\DefaultIcon\"; ValueType: string; ";
print "ValueName: \"\"; ValueData: \"{app}\\gretlw32.exe,2\"\n"; 
print "Root: HKCR; Subkey: \"GretlScriptFile\\shell\\open\\command\"; ";
print "ValueType: string; ValueName: \"\"; ValueData: ";
print "\"\"\"{app}\\gretlw32.exe\"\" -r \"\"%1\"\"\"\n";
