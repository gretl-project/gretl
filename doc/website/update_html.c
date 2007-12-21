#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <locale.h>

#define JACK 0

#if JACK
#define WEBDIR "/home/jack/src/gretl/doc/website"
#define SRCDIR "/home/jack/src/gretl/"
#define WEBSRC "/home/jack/src/gretl/doc/website"
#else
#define WEBDIR "/home/cottrell/stats/esl/website"
#define SRCDIR "/home/cottrell/src"
#define WEBSRC "/home/cottrell/src/doc/website"
#endif

#define MYLEN 96
#define HTMLLEN 1024
#define BUFFER_SIZE 4096

#define NLANGS 4

enum {
    EN,
    ES,
    IT,
    PT
};

char substfile[FILENAME_MAX];

const char *lang_names[] = {
    "en_US",
    "es_ES",
    "it_IT",
    "pt_PT"
};

struct lang_strings_t {
    char lang[8];
    char longdate[48];
    char shortdate[16];
};

typedef struct _gretl_version gretl_version;
struct _gretl_version {
    int major;
    int minor;
    int rev;
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
	if (!strcmp(pat, "LONGDATE")) {
	    printf("got '%s' for '%s' (lang=%s)\n", 
		   lang_strings[i].longdate,
		   pat, lang_strings[i].lang);
	} else if (!strcmp(pat, "SHORTDATE")) {
	    printf("got '%s' for '%s' (lang=%s)\n", 
		   lang_strings[i].shortdate,
		   pat, lang_strings[i].lang);
	}	    
    }

    return err;
}

/* read info with which to perform substitutions, from the
   little file "subst" */

int read_subst_file (void)
{
    FILE *fp;
    char line[MYLEN];
    int err = 0;

    fp = fopen(substfile, "r");
    if (fp == NULL) {
	return 1;
    }

    while (!err && fgets(line, sizeof line, fp)) {
	if (!strncmp(line, "lang = ", 7)) {
	    char lang[8];
	    int i;
	    
	    sscanf(line + 7, "%7s", lang);
	    for (i=0; i<NLANGS; i++) {
		if (!strcmp(lang, lang_strings[i].lang)) {
		    err = read_lang_info(line + 13, i);
		}
		if (err) {
		    break;
		}
	    }
	}
    }

    fclose(fp);

    return err;
}

/* print out version and date info to the file "subst" */

int print_subst (const char *verstr)
{
    FILE *fp;
    int i;

    fp = fopen(substfile, "w");
    if (fp == NULL) {
	return 1;
    }

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

/* convert an input date of the form YYYY.MM.DD to an appropriately
   internationalized string like "Mar 14, 2005" */

int get_intl_progdate (const char *s, int i)
{
    int yr, mon, day;
    int err = 0;

    if (sscanf(s, "%d.%d.%d", &yr, &mon, &day) != 3) {
	err = 1;
    } else {
	char pdate[32];
	struct tm tm = {0};

	tm.tm_year = yr - 1900;
	tm.tm_mon = mon - 1;
	tm.tm_mday = day;

	setlocale(LC_TIME, lang_strings[i].lang);
	strftime(pdate, sizeof pdate, "%b %e, %Y\n", &tm);
	setlocale(LC_TIME, "C");

	strcpy(lang_strings[i].shortdate, pdate);
    }

    return err;
}

int syscmd_to_string (const char *syscmd, char *targ, const char *tmpfile)
{
    FILE *fp;
    char line[64];
    int err = 0;

    system(syscmd);

    /* open tmpfile and retrieve the result */

    fp = fopen(tmpfile, "r");

    if (fp == NULL) {
	err = 1;
    } else {
	if (fgets(line, sizeof line, fp) == NULL) {
	    err = 1;
	} else {
	    strcpy(targ, line);
	}

	fclose(fp);
	remove(tmpfile);
    }    

    return err;
}

/* get version and data info organized, and call for this info
   to be written to file */

int make_subst_file (const char *verstr, const char *progdate)
{
    char syscmd[32];
    const char *tmpfile = "tmp.txt";
    int i, err = 0;

    for (i=0; i<NLANGS; i++) {
	char langbit[16];

	if (!strcmp(lang_strings[i].lang, "en_US")) {
	    *langbit = '\0';
	} else {
	    sprintf(langbit, "LANG=%s ", lang_strings[i].lang);
	}

	sprintf(syscmd, "%sdate > %s", langbit, tmpfile);
	err = syscmd_to_string(syscmd, lang_strings[i].longdate, tmpfile);

	if (progdate != NULL) {
	    get_intl_progdate(progdate, i);
	} else {
	    sprintf(syscmd, "%sdate +\"%%b %%e, %%Y\" > %s", langbit, tmpfile);
	    err = syscmd_to_string(syscmd, lang_strings[i].shortdate, tmpfile);
	}
    }

    if (!err) {
	err = print_subst(verstr);
    }

    return err;
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

/* read the little "subst" file for a version number */

int get_version_from_subst_file (char *verstr)
{
    FILE *fp;
    char line[MYLEN];
    int err = 1;

    *verstr = '\0';

    fp = fopen(substfile, "r");
    if (fp == NULL) {
	return 1;
    }

    while (fgets(line, sizeof line, fp)) {
	if (!strncmp(line, "VERSION", 7)) {
	    if (sscanf(line, "%*s %7s", verstr)) {
		err = 0;
		break;
	    }
	}
    }

    fclose(fp);

    return err;
}

char *get_src_version (void)
{
    FILE *fp;
    char fname[MYLEN];
    char line[MYLEN];
    static char verstr[16] = {0};
    char *p;
    int err = 1;

    sprintf(fname, "%s/lib/src/version.h", SRCDIR);
    fp = fopen(fname, "r");
    if (fp == NULL) return NULL;

    while (fgets(line, sizeof line, fp)) {
	if (sscanf(line, "#define GRETL_VERSION \"%15[^\"]", verstr)) {
	    err = 0;
	    break;
	}
    }

    fclose(fp);

    p = strstr(verstr, "cvs");
    if (p != NULL) {
	*p = '\0';
    }

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
    } else if (strstr(fsrc, "portu") || strstr(fsrc, "_pt")) {
	lang = PT;
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
	{ "index_pat.html",           "index.html" },
	{ "win32_pat.html",           "win32/index.html" },
	{ "osx_pat.html",             "osx.html" },	
	{ "gretl_espanol_pat.html",   "gretl_espanol.html" },
	{ "win32_pat_es.html",        "win32/index_es.html", },
	{ "osx_pat_es.html",          "osx_es.html" },
	{ "gretl_italiano_pat.html",  "gretl_italiano.html" },
	{ "win32_pat_it.html",        "win32/index_it.html" },
	{ "osx_pat_it.html",          "osx_it.html" },
	{ "gretl_portugues_pat.html", "gretl_portugues.html" },
	{ "win32_pat_pt.html",        "win32/index_pt.html" },
	{ "osx_pat_pt.html",          "osx_pt.html" },
	{ NULL, NULL }
    };
    struct from_to *ptr = templates;
    char fullsrc[MYLEN];
    char fulltarg[MYLEN];
    int err = 0;

    while (ptr->src != NULL && !err) {

	/* take source from CVS-controlled directory; write 
	   output to web directory */

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

int version_history (char *fname, gretl_version *versions)
{
    int ret, M, m, r, n = 0;
    FILE *clog;
    char line[MYLEN];
    int ye,mo,da;

    clog = fopen(fname, "r");

    while (fgets(line, MYLEN, clog)) {
	ret = sscanf(line, "%d/%d/%d version %d.%d.%d\n", 
		     &da, &mo, &ye, &M, &m, &r);
	if (ret == 6) {
	    versions[n].major = M;
	    versions[n].minor = m;
	    versions[n].rev = r;
	    n++;
	}
    }

    fclose(clog);

    return n;
}

void bug_print_line (char *s, FILE *fp)
{
    int n, bugnum;

    while (*s) {
	n = strspn(s, "0123456789");
	if (n == 7) {
	    sscanf(s, "%d", &bugnum);
	    fprintf(fp, "<a href=\"http://sourceforge.net/tracker/index.php?"
		    "func=detail&aid=%d&group_id=36234&atid=416803\">%d</a>", 
		    bugnum, bugnum);
	    s += n;
	} else {
	    fputc(*s, fp);
	    s++;
	}
    }
}

void write_changelog (char *src, char *targ, gretl_version *gv, int nv)
{
    int ret, M, m, r, i;
    FILE *clog;
    FILE *clogh;
    char line[MYLEN];
    int ye, mo, da;

    clogh = fopen(targ, "a");
    clog = fopen(src, "r");

    fputs("<table>\n", clogh);
    fputs("<colgroup span=\"2\">\n", clogh);
    fputs("<col width=\"150\"></col>\n", clogh);
    fputs("</colgroup>\n", clogh);

    fputs("<a name=\"top\">\n", clogh);
    fputs("<tr>\n<td valign=\"top\">\n<pre>", clogh);

    for (i=0; i<nv; i++) {
	M = gv[i].major;
	m = gv[i].minor;
	r = gv[i].rev;

	fprintf(clogh, "<a href=\"#v%d-%d-%d\">", M, m, r);
	fprintf(clogh, "Version %d.%d.%d</a>\n", M, m, r);
    }

    fputs("</pre>\n</td>\n<td valign=\"top\">\n<pre>\n", clogh);

    while (fgets(line, MYLEN, clog)) {
	ret = sscanf(line, "%d/%d/%d version %d.%d.%d\n", 
		     &da, &mo, &ye, &M, &m, &r);
	if (ret == 6) {
	    fprintf(clogh, "<a name=\"v%d-%d-%d\">", M, m, r); 
	    fprintf(clogh, " %d/%d/%d Version %d.%d.%d</a>", 
		    da, mo, ye, M, m, r);
	    fputs(" [<a href=\"#top\">Back to top</a>]\n", clogh);
	} else {
	    bug_print_line(line, clogh);
	}
    }

    fputs("</pre>\n</td>\n</tr>\n", clogh);

    fclose(clog);
    fclose(clogh);
}

int make_html_changelog (void)
{
    gretl_version gv[1000];
    char targ[MYLEN];
    char src[MYLEN];
    int nv, err = 0;

    sprintf(targ, "%s/ChangeLog.html", WEBDIR);

    sprintf(src, "%s/template/changelog.top", WEBSRC);
    err = copyfile(src, targ, 0);
    if (err) return err;

    sprintf(src, "%s/ChangeLog", SRCDIR);

    nv = version_history(src, gv);

    write_changelog(src, targ, gv, nv);

    sprintf(src, "%s/template/changelog.end", WEBSRC);
    err = copyfile(src, targ, 1);

    return err;
}

int main (int argc, char **argv)
{
    char *progdate;
    char *src_version;
    char html_version[8];
    int up_to_date;
    int err = 0;

    if (argc < 2) {
	fprintf(stderr, "%s: Please give a gretl version number,\n"
		"or say \"auto\" to read version info from the source tree.\n",
		argv[0]);
	fputs("You can specify a program date, YYYY.MM.DD, as a second arg\n", stderr);
	exit(EXIT_FAILURE);
    }

    if (strcmp(argv[1], "auto")) {
	src_version = argv[1];
    } else {
	src_version = get_src_version();
    }

    if (argc == 3) {
	progdate = argv[2];
    } else {
	progdate = NULL;
    }

    if (src_version == NULL) {
	fputs("Couldn't find source version\n", stderr);
	exit(EXIT_FAILURE);
    }

    sprintf(substfile, "%s/subst", WEBDIR);

    get_version_from_subst_file(html_version);
    if (*html_version == '\0') {
	fprintf(stderr, "Couldn't get html version from %s\n", substfile);
	strcpy(html_version, "unknown");
    }

    fprintf(stderr, "source version is %s\n", src_version);
    fprintf(stderr, "html version is %s\n", html_version);

    if (strcmp(src_version, html_version)) {
	up_to_date = 0;
	fputs("version number needs updating\n", stderr);
    } else {
	up_to_date = 1;
	fputs("version number is current\n", stderr);
    }

    lang_strings_init();

    if (up_to_date) {
	err = read_subst_file();
    } else {
	err = make_subst_file(src_version, progdate);
    }

    if (err) {
	fputs("Error doing stuff\n", stderr);
    } else {
	err = process_templates(src_version);
    }

    if (!err) {
	err = make_html_changelog();
    }

    return err;
}
