#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* filter: convert MANIFEST file to .iss script for making Win32 installer */

void usage (const char *s)
{
    fprintf(stderr, "usage: %s < MANIFEST\n", s);
    exit(EXIT_FAILURE);
}

void define_program_icons (void)
{
    printf("\n[Icons]\n");
    printf("Name: \"{group}\\gretl\"; Filename: \"{app}\\gretlw32.exe\"\n");
    printf("Name: \"{group}\\Gretl Web Site\"; Filename: \"{app}\\gretl_website.url\"\n");
    printf("Name: \"{group}\\gretl updater\"; Filename: \"{app}\\gretl_updater.exe\"\n");
    printf("Name: \"{group}\\uninstall gretl\"; Filename: \"{app}\\unins000.exe\"\n");
    printf("Name: \"{userdesktop}\\gretl\"; Filename: \"{app}\\gretlw32.exe\"; WorkingDir: \"{app}\"\n");
}

void LM_entry (const char *valname, const char *valdata)
{
    printf("Root: HKLM; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: "
	   "\"%s\"; ValueData: \"%s\"\n", valname, valdata);
}

void CU_entry (const char *valname, const char *valdata)
{
    printf("Root: HKCU; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: "
	   "\"%s\"; ValueData: \"%s\"\n", valname, valdata);
}

void reg_suffix (const char *sfx, const char *name, const char *mime,
		 const char *descrip, int iconnum)
{
    /* map from suffix to filetype name */
    printf("Root: HKCR; Subkey: \"%s\"; ValueType: string; ValueName: "
	   "\"\"; ValueData: \"%s\"; Flags: uninsdeletevalue\n",
	   sfx, name);
    /* map from suffix to mime type */
    printf("Root: HKCR; Subkey: \"%s\"; ValueType: string; "
	   "ValueName: \"Content Type\"; ValueData: \"%s\"\n",
	   sfx, mime);
    /* map from filetype name to description */
    printf("Root: HKCR; Subkey: \"%s\"; ValueType: string; ValueName: "
	   "\"\"; ValueData: \"%s\"; Flags: uninsdeletekey\n",
	   name, descrip);
    /* map from filetype name to icon */
    printf("Root: HKCR; Subkey: \"%s\\DefaultIcon\"; ValueType: string; "
	   "ValueName: \"\"; ValueData: \"{app}\\gretlw32.exe,%d\"\n",
	   name, iconnum);
    /* map from filetype name to 'open' action */
    if (iconnum == 1) {
	printf("Root: HKCR; Subkey: \"%s\\shell\\open\\command\"; "
	       "ValueType: string; ValueName: \"\"; ValueData: "
	       "\"\"\"{app}\\gretlw32.exe\"\" \"\"%%1\"\"\"\n", name);
    } else {
	printf("Root: HKCR; Subkey: \"%s\\shell\\open\\command\"; "
	       "ValueType: string; ValueName: \"\"; ValueData: "
	       "\"\"\"{app}\\gretlw32.exe\"\" -r \"\"%%1\"\"\"\n", name);
    }
}

void set_registry_entries (void)
{
    printf("\n[Registry]\n");

    /* base paths */
    printf("; Start \"gretl\" registry keys.\n");
    printf("Root: HKCR; Subkey: \"Software\\gretl\"; Flags: uninsdeletekey\n");
    printf("Root: HKCU; Subkey: \"Software\\gretl\"; Flags: uninsdeletekey\n");
    printf("Root: HKLM; Subkey: \"Software\\gretl\"; Flags: uninsdeletekey\n");

    /* specific entries, Local Machine */
    LM_entry("gretldir", "{app}");
    LM_entry("Rcommand", "RGui.exe");
    LM_entry("viewdvi", "windvi.exe");

    /* specific entries, Current User */
    CU_entry("binbase",    "{app}\\db\\");
    CU_entry("ratsbase",   "f:\\");
    CU_entry("dbhost",     "ricardo.ecn.wfu.edu");
    CU_entry("dbproxy",    "");
    CU_entry("useproxy",   "false");
    CU_entry("updater",    "false");
    CU_entry("calculator", "calc.exe");
    CU_entry("Fixed_font", "Courier New 10");
    CU_entry("App_font",   "");
    CU_entry("Png_font",   "verdana 8");
    CU_entry("Gp_colors",  "");
    CU_entry("DataPage",   "Gretl");
    CU_entry("ScriptPage", "Gretl");
    CU_entry("manpref",    "0");

    /* Establish file associations */

    reg_suffix(".gdt", "GretlDataFile", "application/x-gretldata", 
	       "Gretl data file", 1);

    reg_suffix(".gretl", "GretlSessionFile", "application/x-gretlsession", 
	       "Gretl session file", 2);

    reg_suffix(".inp", "GretlScriptFile", "application/x-gretlscript", 
	       "Gretl script file", 2);
}

void preamble (const char *s)
{
    printf("; -- gretl.iss --\n");
    printf("\n[Setup]\n");
    printf("AppName=gretl\n");
    printf("AppVerName=gretl version %s\n", s);
    printf("AppVersion=%s\n", s);
    printf("AppPublisherURL=http://gretl.sourceforge.net/\n");
    printf("AppSupportURL=http://gretl.sourceforge.net/\n");
    printf("DefaultDirName={pf}\\gretl\n");
    printf("DefaultGroupName=gretl\n");
    printf("PrivilegesRequired=poweruser\n");
    printf("UninstallDisplayIcon={app}\\gretlw32.exe\n");
    printf("ChangesAssociations=yes\n");
    printf("DirExistsWarning=no\n");

    printf("\n[InstallDelete]\n");
    printf("Type: files; Name: \"{app}\\*.dll\"\n");
}

void tailstrip (char *s)
{
    int i, n = strlen(s);

    for (i=n-1; i>0; i--) {
	if (s[i] == '\n' || s[i] == '\r' || s[i] == ' ') {
	    s[i] = '\0';
	} else {
	    break;
	}
    }
}

void modpath (char *s)
{
    while (*s) {
	if (*s == '/') {
	    *s = '\\';
	} 
	s++;
    }
}

#define MAXBITS 8

int split_path (char *s, char **p)
{
    int i;

    for (i=0; i<MAXBITS; i++) {
	p[i] = strtok((i == 0)? s : NULL, "\\");
	if (p[i] == NULL) {
	    break;
	}
    }

    return i;
}

int main (int argc, char **argv)
{
    char *pathbits[MAXBITS];
    char line[512];
    char version[16];
    char *path;
    int i, n;

    if (fgets(line, sizeof line, stdin) == NULL) {
	usage(argv[0]);
    }

    if (strlen(line) < 11) {
	fprintf(stderr, "malformed MANIFEST: expected VERSION ...\n");
	exit(EXIT_FAILURE);
    }

    *version = '\0';
    strncat(version, line + 8, 15);
    tailstrip(version);

    fprintf(stderr, "Making installer script for gretl version %s...\n",
	    version);
    preamble(version);

    printf("\n[Files]\n");

    /* Read MANIFEST from stdin.  Format is size date time pathname,
       for example:

       36352 2009-02-22 14:12 gretl/libprob.dll

    */
    while (fgets(line, sizeof line, stdin)) {
	if (strstr(line, "VERSION") || strstr(line, "DATE")) {
	    continue;
	}
	tailstrip(line);
	/* path is last space-separated field */
	path = strrchr(line, ' ');
	if (path == NULL) {
 	    continue;
	}
	path++;
	/* substitute slash -> backslash in line */
	modpath(path);
	printf("Source: \"%s\"; ", path); 
	n = split_path(path, pathbits);
	printf("Destdir: \"{app}");
	for (i=1; i<n-1; i++) {
	    printf("\\%s", pathbits[i]);
	}
	printf("\"\n");
    }

    define_program_icons();
    set_registry_entries();

    return 0;
}


