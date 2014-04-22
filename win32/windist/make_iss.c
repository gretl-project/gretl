#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* filter: convert MANIFEST file to .iss script for making Windows installer */

void usage (const char *s)
{
    fprintf(stderr, "usage: %s < MANIFEST\n", s);
    exit(EXIT_FAILURE);
}

struct isl {
    const char *id;
    const char *fname;
};

void add_languages (void)
{
    struct isl languages[] = {
	{ "en", "Default" },
	/* { "eu", "Basque" }, out of date translation */
	{ "eu", "Spanish" },
	{ "gl", "Spanish" },
	{ "pt_BR", "BrazilianPortuguese" },
	{ "cs", "Czech" },
	{ "fr", "French" },
	{ "de", "German" },
	{ "it", "Italian" },
	{ "pl", "Polish" },
	{ "pt", "Portuguese" },
	{ "ru", "Russian" },
	{ "es", "Spanish" },
	{ "tr", "Turkish" },
	{ "ca", "Catalan" },
	{ "el", "Greek" },
	{ NULL, NULL }
    };
    int i;

    printf("\n[Languages]\n");

    printf("Name: \"%s\"; MessagesFile: \"compiler:%s.isl\"\n", 
	   languages[0].id, languages[0].fname);

    for (i=1; languages[i].id != NULL; i++) {
	printf("Name: \"%s\"; MessagesFile: \"compiler:Languages\\%s.isl\"\n", 
	       languages[i].id, languages[i].fname);
    }
}

void define_program_icons (void)
{
    printf("\n[Icons]\n");
    printf("Name: \"{group}\\gretl\"; Filename: \"{app}\\gretl.exe\"\n");
    printf("Name: \"{group}\\Gretl Web Site\"; Filename: \"{app}\\gretl_website.url\"\n");
    printf("Name: \"{group}\\uninstall gretl\"; Filename: \"{app}\\unins000.exe\"\n");
    printf("Name: \"{userdesktop}\\gretl\"; Filename: \"{app}\\gretl.exe\"; WorkingDir: \"{app}\"\n");
}

void LM_entry (const char *valname, const char *valdata)
{
    printf("Root: HKLM; Subkey: \"Software\\gretl\"; ValueType: string; ValueName: "
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
	   "ValueName: \"\"; ValueData: \"{app}\\gretl.exe,%d\"\n",
	   name, iconnum);
    /* map from filetype name to 'open' action */
    if (iconnum == 2) {
	/* script: use "-r" */
	printf("Root: HKCR; Subkey: \"%s\\shell\\open\\command\"; "
	       "ValueType: string; ValueName: \"\"; ValueData: "
	       "\"\"\"{app}\\gretl.exe\"\" -r \"\"%%1\"\"\"\n", name);
    } else {
	/* data or session file */
	printf("Root: HKCR; Subkey: \"%s\\shell\\open\\command\"; "
	       "ValueType: string; ValueName: \"\"; ValueData: "
	       "\"\"\"{app}\\gretl.exe\"\" \"\"%%1\"\"\"\n", name);
    }
}

void set_registry_entries (void)
{
    printf("\n[Registry]\n");

    /* base paths */
    printf("; Start \"gretl\" registry keys.\n");
    printf("Root: HKCR; Subkey: \"Software\\gretl\"; Flags: uninsdeletekey\n");
    printf("Root: HKLM; Subkey: \"Software\\gretl\"; Flags: uninsdeletekey\n");

    /* specific entries, Local Machine */
    LM_entry("gretldir", "{app}");

    /* Establish file associations */
    reg_suffix(".gdt", "GretlDataFile", "application/x-gretldata", 
	       "Gretl data file", 1);
    reg_suffix(".inp", "GretlScriptFile", "application/x-gretlscript", 
	       "Gretl script file", 2);
    reg_suffix(".gretl", "GretlSessionFile", "application/x-gretlsession", 
	       "Gretl session file", 3);
    reg_suffix(".gdtb", "GretlBindataFile", "application/x-gretlbindata", 
	       "Gretl binary data file", 4);
}

void preamble (const char *s, int x64)
{
    printf("; -- gretl.iss --\n");
    printf("\n[Setup]\n");
    printf("AppName=gretl\n");
    if (x64) {
	printf("AppVerName=gretl version %s (x86_64)\n", s);
    } else {
	printf("AppVerName=gretl version %s\n", s);
    }
    printf("AppVersion=%s\n", s);
    if (x64) {
	printf("ArchitecturesAllowed=x64\n");
	printf("ArchitecturesInstallIn64BitMode=x64\n");
    }
    printf("AppPublisher=The gretl team\n");
    printf("AppPublisherURL=http://gretl.sourceforge.net/\n");
    printf("AppSupportURL=http://gretl.sourceforge.net/\n");
    printf("DefaultDirName={pf}\\gretl\n");
    printf("DefaultGroupName=gretl\n");
    printf("PrivilegesRequired=poweruser\n");
    printf("UninstallDisplayIcon={app}\\gretl.exe\n");
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
    char arch[8];
    char *path;
    int x64 = 0;
    int i, n;

    if (fgets(line, sizeof line, stdin) == NULL) {
	usage(argv[0]);
    }

    if (strlen(line) < 11) {
	fprintf(stderr, "malformed MANIFEST: expected VERSION ...\n");
	exit(EXIT_FAILURE);
    }

    *version = *arch = '\0';

    n = sscanf(line + 8, "%15s %7s", version, arch);

    if (n < 1) {
	fprintf(stderr, "malformed MANIFEST: expected VERSION ...\n");
	exit(EXIT_FAILURE);
    } else if (n == 2) {
	if (!strcmp(arch, "x64")) {
	    x64 = 1;
	    fprintf(stderr, "Making installer script for gretl version %s (x64)...\n",
		    version);
	} else {
	    fprintf(stderr, "malformed MANIFEST: if arch is present, it must be \"x64\"\n");
	    exit(EXIT_FAILURE);
	}
    } else {
	fprintf(stderr, "Making installer script for gretl version %s...\n",
		version);
    }	

    preamble(version, x64);

    add_languages();

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


