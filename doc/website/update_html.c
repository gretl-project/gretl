#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define WEBDIR "/home/allin/stats/esl/website"
#define SRCDIR "/home/allin/src"
#define WEBSRC "/home/allin/src/doc/website"

#define MYLEN 96
#define HTMLLEN 1024
#define BUFFER_SIZE 4096

#define NLANGS 3

enum {
    EN,
    ES,
    IT
};

const char *lang_names[] = {
    "en_US",
    "es_ES",
    "it_IT"
};

struct lang_strings_t {
    char lang[8];
    char longdate[48];
    char shortdate[16];
};

struct lang_strings_t lang_strings[NLANGS];

void lang_strings_init (void)
{
    int i;

    for (i=0; i<NLANGS; i++) {
	strcpy(lang_strings[i].lang, lang_names[i]);
	lang_strings[i].longdate[0] = '\0';
	lang_strings[i].shortdate[0] = '\0';
    }
}

void tail_strip (char *s)
{
    int i, n = strlen(s) - 1;

    for (i=n; i>0; i--) {
	if (s[i] == '\n' || s[i] == ' ') {
	    s[i] = '\0';
	} else {
	    break;
	}
    }
}

int read_lang_info (const char *s, int i)
{
    char pat[16];
    int err = 0;

    if (!sscanf(s, "%15s", pat)) return 1;

    if (!strcmp(pat, "LONGDATE")) {
	strcpy(lang_strings[i].longdate, s + strlen(pat) + 1);
	tail_strip(lang_strings[i].longdate);
    } else if (!strcmp(pat, "SHORTDATE")) {
	strcpy(lang_strings[i].shortdate, s + strlen(pat) + 1);
	tail_strip(lang_strings[i].shortdate);
    } else {
	err = 1;
    }

    if (!err) {
	printf("got '%s' for '%s' (lang=%s)\n", 
	       lang_strings[i].longdate,
	       pat, lang_strings[i].lang);
    }

    return err;
}

int read_subst_file (void)
{
    FILE *fp;
    char fname[MYLEN];
    char line[MYLEN];
    int err = 0;

    sprintf(fname, "%s/subst", WEBDIR);
    fp = fopen(fname, "r");
    if (fp == NULL) return 1;

    while (!err && fgets(line, sizeof line, fp)) {
	if (!strncmp(line, "lang = ", 7)) {
	    char lang[8];
	    int i;
	    
	    sscanf(line + 7, "%7s", lang);
	    for (i=0; i<NLANGS; i++) {
		if (!strcmp(lang, lang_strings[i].lang)) {
		    err = read_lang_info(line + 13, i);
		}
		if (err) break;
	    }
	}
    }

    fclose(fp);

    return 0;
}

int print_subst (const char *verstr)
{
    FILE *fp;
    char fname[MYLEN];
    int i;

    sprintf(fname, "%s/subst", WEBDIR);
    fp = fopen(fname, "w");
    if (fp == NULL) return 1;

    fprintf(fp, "VERSION %s\n", verstr);

    for (i=0; i<NLANGS; i++) {
	fprintf(fp, "lang = %s LONGDATE %s", lang_strings[i].lang,
		lang_strings[i].longdate);
	fprintf(fp, "lang = %s SHORTDATE %s", lang_strings[i].lang,
		lang_strings[i].shortdate);
    }

    fclose(fp);

    return 0;
}

int make_subst_file (const char *verstr)
{
    FILE *fp;
    char syscmd[32];
    char line[48];
    const char *tmpfile = "tmp.txt";
    int i;

    for (i=0; i<NLANGS; i++) {
	char langbit[16];

	if (!strcmp(lang_strings[i].lang, "en_US")) {
	    *langbit = '\0';
	} else {
	    sprintf(langbit, "LANG=%s ", lang_strings[i].lang);
	}

	sprintf(syscmd, "%sdate > %s", langbit, tmpfile);
	system(syscmd);
	fp = fopen(tmpfile, "r");
	if (fp == NULL) return 1;
	if (fgets(line, sizeof line, fp) == NULL) {
	    fclose(fp);
	    return 1;
	} else {
	    strcpy(lang_strings[i].longdate, line);
	}
	fclose(fp);
	remove(tmpfile);

	sprintf(syscmd, "%sdate +\"%%b %%e, %%Y\" > %s", langbit, tmpfile);
	system(syscmd);
	fp = fopen(tmpfile, "r");
	if (fp == NULL) return 1;
	if (fgets(line, sizeof line, fp) == NULL) {
	    fclose(fp);
	    return 1;
	} else {
	    strcpy(lang_strings[i].shortdate, line);
	}
	fclose(fp);
	remove(tmpfile);
    }

    return print_subst(verstr);
}

int copyfile (const char *src, const char *dest, int append) 
{
    FILE *srcfd, *destfd;
    char buf[BUFFER_SIZE];
    size_t n;

    if (!strcmp(src, dest)) return 1;
   
    if ((srcfd = fopen(src, "rb")) == NULL) {
        fprintf(stderr, "Couldn't open %s\n", src);
        return 1; 
    }

    if (append) {
	destfd = fopen(dest, "ab");
    } else {
	destfd = fopen(dest, "wb");
    }

    if (destfd == NULL) {
        fprintf(stderr, "Couldn't write to %s\n", dest);
        fclose(srcfd);
        return 1;
    }

    while ((n = fread(buf, 1, sizeof buf, srcfd)) > 0) {
        fwrite(buf, 1, n, destfd);
    }

    fclose(srcfd);
    fclose(destfd);

    return 0;
}

char *get_html_version (void)
{
    FILE *fp;
    char fname[MYLEN];
    char line[MYLEN];
    static char verstr[8];
    int err = 1;

    sprintf(fname, "%s/subst", WEBDIR);
    fp = fopen(fname, "r");
    if (fp == NULL) return NULL;

    while (fgets(line, sizeof line, fp)) {
	if (!strncmp(line, "VERSION", 7)) {
	    if (sscanf(line, "%*s %7s", verstr)) {
		err = 0;
		break;
	    }
	}
    }

    fclose(fp);

    return (err)? NULL : verstr;
}

char *get_src_version (void)
{
    FILE *fp;
    char fname[MYLEN];
    char line[MYLEN];
    static char verstr[8];
    int err = 1;

    sprintf(fname, "%s/lib/src/version.h", SRCDIR);
    fp = fopen(fname, "r");
    if (fp == NULL) return NULL;

    while (fgets(line, sizeof line, fp)) {
	if (strstr(line, "version_string")) {
	    if (sscanf(line, "const char *version_string = \"%7[^\"]", 
		       verstr)) {
		err = 0;
		break;
	    }
	}
    }

    fclose(fp);

    return (err)? NULL : verstr;
}

int make_substitutions (char *line, int lang, char *verstr)
{
    char save[HTMLLEN];
    char *p, *sub;

    p = line;
    while ((p = strstr(p, "VERSION")) != NULL) {
	sub = verstr;
	strcpy(save, p + strlen("VERSION"));
	strcpy(p, sub);
	strcpy(p + strlen(sub), save);
    }

    p = line;
    while ((p = strstr(p, "LONGDATE")) != NULL) {
	sub = lang_strings[lang].longdate;
	strcpy(save, p + strlen("LONGDATE"));
	strcpy(p, sub);
	strcpy(p + strlen(sub), save);
    }    

    p = line;
    while ((p = strstr(p, "SHORTDATE")) != NULL) {
	sub = lang_strings[lang].shortdate;
	strcpy(save, p + strlen("SHORTDATE"));
	strcpy(p, sub);
	strcpy(p + strlen(sub), save);
    }
    
    return 0;
}

int copy_substitute (const char *fsrc, const char *ftarg,
		     char *verstr)
{
    FILE *fin, *fout;
    char line[HTMLLEN];
    int lang;

    if (strstr(fsrc, "espan") || strstr(fsrc, "_es")) {
	lang = ES;
    } else if (strstr(fsrc, "itali") || strstr(fsrc, "_it")) {
	lang = IT;
    } else {
	lang = EN;
    }

    fin = fopen(fsrc, "r");
    if (fin == NULL) {
	fprintf(stderr, "Couldn't open %s for reading\n", fsrc);
	return 1;
    }

    fout = fopen(ftarg, "w");
    if (fout == NULL) {
	fprintf(stderr, "Couldn't open %s for reading\n", ftarg);
	fclose(fin);
	return 1;
    }

    while (fgets(line, sizeof line, fin)) {
	make_substitutions(line, lang, verstr);
	fputs(line, fout);
    }

    fclose(fin);
    fclose(fout);

    return 0;
}

struct from_to {
    char *src;
    char *targ;
};

int process_templates (char *verstr)
{
    struct from_to templates[] = {
	{ "index_pat.html",          "index.html" },
	{ "win32_pat.html",          "win32/index.html" },
	{ "osx_pat.html",            "osx.html" },	
	{ "gretl_espanol_pat.html",  "gretl_espanol.html" },
	{ "win32_pat_es.html",       "win32/index_es.html", },
	{ "osx_pat_es.html",         "osx_es.html" },
	{ "gretl_italiano_pat.html", "gretl_italiano.html" },
	{ "win32_pat_it.html",       "win32/index_it.html" },
	{ "osx_pat_it.html",         "osx_it.html" },
	{ NULL, NULL }
    };
    struct from_to *ptr = templates;
    char fullsrc[MYLEN];
    char fulltarg[MYLEN];
    int err = 0;

    while (ptr->src && !err) {

	sprintf(fullsrc, "%s/template/%s", WEBSRC, ptr->src);
	sprintf(fulltarg, "%s/%s", WEBDIR, ptr->targ);

	fprintf(stderr, "substituting from\n %s to\n %s\n", 
		fullsrc, fulltarg);

	err = copy_substitute(fullsrc, fulltarg, verstr);

	if (err) {
	    fputs("Error in process_templates\n", stderr);
	} else {
	    ptr++;
	}
    }
    
    return err;
}

int make_html_changelog (void)
{
    char targ[MYLEN];
    char src[MYLEN];
    int err = 0;

    sprintf(targ, "%s/ChangeLog.html", WEBDIR);
    
    sprintf(src, "%s/template/changelog.top", WEBSRC);
    err = copyfile(src, targ, 0);
    if (err) return err;

    sprintf(src, "%s/ChangeLog", SRCDIR);
    err = copyfile(src, targ, 1);
    if (err) return err;

    sprintf(src, "%s/template/changelog.end", WEBSRC);
    err = copyfile(src, targ, 1);

    return err;
}

int main (void)
{
    char *html_version;
    char *src_version;
    int up_to_date;
    int err = 0;

    src_version = get_src_version();
    if (src_version == NULL) {
	fputs("Couldn't find source version\n", stderr);
	exit(EXIT_FAILURE);
    }

    html_version = get_html_version();
    if (html_version == NULL) {
	fputs("Couldn't find html version\n", stderr);
	exit(EXIT_FAILURE);
    }

    fprintf(stderr, "source version is %s\n", src_version);
    fprintf(stderr, "html version is %s\n", html_version);

    if (strcmp(src_version, html_version)) {
	up_to_date = 0;
	fputs("version needs updating\n", stderr);
    } else {
	up_to_date = 1;
	fputs("version is current\n", stderr);
    }

    lang_strings_init();

    if (up_to_date) {
	err = read_subst_file();
    } else {
	err = make_subst_file(src_version);
    }

    if (err) {
	fputs("Error doing stuff\n", stderr);
    } else {
	err = process_templates(src_version);
    }

    if (!err && !up_to_date) {
	err = make_html_changelog();
    }

    return err;
}
