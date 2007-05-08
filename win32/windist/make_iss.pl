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
print "\n[Setup]\n";
print "AppName=gretl\n";
print "AppVerName=gretl version $verstr\n";
print "AppPublisherURL=http://gretl.sourceforge.net/\n";
print "AppSupportURL=http://gretl.sourceforge.net/\n";
print "AppCopyright=Copyright C 1999-2006 Allin Cottrell\n";
print "DefaultDirName={pf}\\gretl\n";
print "DefaultGroupName=gretl\n";
print "UninstallDisplayIcon={app}\\gretlw32.exe\n";
print "ChangesAssociations=yes\n";
print "DirExistsWarning=no\n";

print "\n[INI]\n";
print "Filename: \"{win}\\WGNUPLOT.INI\"; Section: \"WGNUPLOT\"; Flags: createkeyifdoesntexist\n";
print "Filename: \"{win}\\WGNUPLOT.INI\"; Section: \"WGNUPLOT\"; Key: \"TextFont\"; String: \"Courier New, 9\"\n";

print "\n[InstallDelete]\n";
print "Type: files; Name: \"{app}\\*.dll\"\n";

print "\n[Dirs]\n";
print "Name: {code:GetDataDir}; Flags: uninsneveruninstall\n";

print "\n[Files]\n";

while ($line = <STDIN>) {
    if ($line =~ /VERSION/ || $line =~ /DATE/) {
        next;
    }   
    @pathbits = split(/ +/, $line); 
    $line = $pathbits[@pathbits - 1];
    $line =~ s+/+\\+g;
    chomp($line);
    @pathbits = split(/\\/, $line);
    print "Source: \"$line\"; "; 
    if ($line =~ /\.ini$/) {
	print "Destdir: \"{win}";
    } elsif ($line =~ /README/) {
        print "Destdir: \"{code:GetDataDir}";
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
print "Name: \"{group}\\uninstall gretl\"; Filename: \"{app}\\unins000.exe\"\n";
print "Name: \"{userdesktop}\\gretl\"; Filename: \"{app}\\gretlw32.exe\"; WorkingDir: \"{app}\"\n";

print "\n[Registry]\n";
print "; Start \"gretl\" keys under HKEY_CURRENT_USER and HKEY_CLASSES_ROOT.\n"; 
print "Root: HKCR; Subkey: \"Software\\gretl\"; Flags: uninsdeletekey\n"; 
print "Root: HKCU; Subkey: \"Software\\gretl\"; Flags: uninsdeletekey\n";
print "Root: HKLM; Subkey: \"Software\\gretl\"; Flags: uninsdeletekey\n";
print "Root: HKCR; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"gretldir\"; ValueData: \"{app}\"\n"; 
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"userdir\"; ValueData: \"{code:GetDataDir}\"\n"; 
print "Root: HKLM; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"gnuplot\"; ValueData: \"{app}\\wgnuplot.exe\"\n";
print "Root: HKLM; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"Rcommand\"; ValueData: \"RGui.exe\"\n";
print "Root: HKLM; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"viewdvi\"; ValueData: \"windvi.exe\"\n";
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"expert\"; ValueData: \"false\"\n";    
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"binbase\"; ValueData: \"{app}\\db\\\"\n";   
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"ratsbase\"; ValueData: \"f:\\\"\n";   
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"dbhost\"; ValueData: \"ricardo.ecn.wfu.edu\"\n"; 
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"dbproxy\"; ValueData: \"\"\n"; 
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"useproxy\"; ValueData: \"false\"\n"; 
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"updater\"; ValueData: \"false\"\n";
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"calculator\"; ValueData: \"calc.exe\"\n";   
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"toolbar\"; ValueData: \"true\"\n"; 
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"Fixed_font\"; ValueData: \"Courier New 10\"\n";   
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"App_font\"; ValueData: \"\"\n";   
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"Png_font\"; ValueData: \"verdana 8\"\n"; 
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"Gp_colors\"; ValueData: \"\"\n";
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"DataPage\"; ValueData: \"Ramanathan\"\n";
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"ScriptPage\"; ValueData: \"Ramanathan\"\n";
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"manpref\"; ValueData: \"0\"\n";
print "Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: ";
print "\"useqr\"; ValueData: \"false\"\n";
 
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

print "\n[Code]\n";
print "var\n";
print "  DataDirPage: TInputDirWizardPage;\n";
print "\n";  
print "procedure InitializeWizard;\n";
print "begin\n";
print "  { Create the pages }\n";
print "   DataDirPage := CreateInputDirPage(wpSelectDir,\n";
print "                  'Select Personal Data Directory', \n";
print "                  'Where should personal data files be kept?',\n";
print "                  'Select the folder for storing personal data files, then click Next.',\n";
print "                   False, '');\n";
print "  DataDirPage.Add('');\n";
print "  DataDirPage.Values[0] := GetPreviousData('DataDir', '');\n";
print "  if DataDirPage.Values[0] = '' then\n";
print "     DataDirPage.Values[0] := ExpandConstant('{userappdata}\\gretl');\n";
print "end;\n";
print "\n";
print "procedure RegisterPreviousData(PreviousDataKey: Integer);\n";
print "begin\n";
print "  { Store the settings so we can restore them next time }\n";
print "  SetPreviousData(PreviousDataKey, 'DataDir', DataDirPage.Values[0]);\n";
print "end;\n";
print "\n";
print "function UpdateReadyMemo(Space, NewLine, MemoUserInfoInfo, MemoDirInfo, MemoTypeInfo,\n";
print "  MemoComponentsInfo, MemoGroupInfo, MemoTasksInfo: String): String;\n";
print "var\n";
print "  S: String;\n";
print "begin\n";
print "  { Fill the 'Ready Memo' with the normal settings and the custom settings }\n";
print "  S := '';\n";
print "  S := S + MemoDirInfo + NewLine;\n";
print "  S := S + Space + DataDirPage.Values[0] + ' (personal data files)' + NewLine;\n";
print "  Result := S;\n";
print "end;\n\n";
print "function GetDataDir(Param: String): String;\n";
print "begin\n";
print "  { Return the selected DataDir }\n";
print "  Result := DataDirPage.Values[0];\n";
print "end;\n";

