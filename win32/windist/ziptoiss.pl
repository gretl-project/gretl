#!/usr/bin/perl

# convert zipfile listing to .iss script for making Windows
# installer

use strict;

my $zname;
my $ldash;
my $rdash;
my $verstr;
my $line; 
my $subdir;
my @fields;
my @pathbits;
my $i;

$zname = `/bin/ls -1 *.zip`;
chomp($zname);
$ldash = index($zname, "-") + 1;
$rdash = rindex($zname, "-");
$verstr = substr($zname, $ldash, $rdash - $ldash);
print STDERR "processing $zname (gretl version $verstr)...\n";

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

print "[Dirs]\n";
print "Name: \"{app}\\user\"\n";
print "Name: \"{app}\\data\"\n";
print "Name: \"{app}\\scripts\"\n";
print "Name: \"{app}\\db\"\n";
print "Name: \"{app}\\locale\"\n";
print "Name: \"{app}\\locale\\es\"\n";
print "Name: \"{app}\\locale\\es\\LC_MESSAGES\"\n\n";

print "[Files]\n";
while ($line = <>) {
   $line =~ s+/+\\+g;
   chomp($line);
   if ($line =~ /\\$/ || $line =~ /----/ || 
       $line =~ /Length/ || $line =~ /Archive/ ||
       $line =~ / files/) { 
      next; 
   }
   @fields = split(/ +/, $line);
   @pathbits = split(/\\/, $fields[4]);
   print "Source: \"c:\\userdata\\$fields[4]\"; "; 
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

# print "\n[Run]\n";
# print "Filename: \"{app}\\makecfg.exe\"; Parameters: \"\"\"{app}\"\"\"\n";

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
print "\"fontspec\"; ValueData: \"-unknown-Courier ";
print "New-normal-r-normal-*-*-120-*-*-m-*-iso8859-1\"\n";
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
