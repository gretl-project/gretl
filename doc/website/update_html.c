#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <locale.h>
#include <ctype.h>

#define MYLEN 512
#define HTMLLEN 1024
#define BUFFER_SIZE 4096

#define NLANGS 6

static char WEBDIR[256];
static char SRCDIR[256];
static char WEBSRC[256];

enum {
    EN,
    ES,
    IT,
    PT,
    TR,
    RU
};

const char *lang_names[] = {
    "en_US",
    "es_ES",
    "it_IT",
    "pt_PT",
    "tr_TR",
    "ru_RU"
};

struct lang_strings_t {
    char lang[16];
    char longdate[48];
    char shortdate[32];
};

typedef struct _gretl_version gretl_version;

struct _gretl_version {
    int major;
    int minor;
    int rev;
    char letter;
    char pre;
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

int read_subst_file_full (const char *fname)
{
    FILE *fp;
    char line[MYLEN];
    int err = 0;

    fp = fopen(fname, "r");
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

int print_subst (const char *fname, const char *verstr)
{
    FILE *fp;
    int i;

    fp = fopen(fname, "w");
    if (fp == NULL) {
        fprintf(stderr, "Couldn't open '%s'\n", fname);
        return 1;
    }

    fprintf(stderr, "opened '%s' for writing\n", fname);

    fprintf(fp, "VERSION %s\n", verstr);

    for (i=0; i<NLANGS; i++) {
        fprintf(fp, "lang = %s LONGDATE %s\n", lang_strings[i].lang,
                lang_strings[i].longdate);
        fprintf(fp, "lang = %s SHORTDATE %s\n", lang_strings[i].lang,
                lang_strings[i].shortdate);
    }

    fclose(fp);

    return 0;
}

/* convert an input date of the form YYYY-MM-DD to an appropriately
   internationalized string like "Mar 14, 2005" */

int get_intl_progdate (struct tm *tm, int i)
{
    char locale[32];
    char pdate[32];

    sprintf(locale, "%s.UTF-8", lang_strings[i].lang);
    setlocale(LC_TIME, locale);

    if (!strncmp(lang_strings[i].lang, "ru", 2)) {
        strftime(pdate, sizeof pdate, "%b %e, %Y г.", tm);
    } else {
        strftime(pdate, sizeof pdate, "%b %e, %Y", tm);
    }

    strcpy(lang_strings[i].shortdate, pdate);
    setlocale(LC_TIME, "C");

    return 0;
}

int make_subst_file (const char *fname,
                     const char *release_version,
                     const char *release_date)
{
    struct tm tm = {0};
    int y, m, d;
    int i, err = 0;

    if (sscanf(release_date, "%d-%d-%d", &y, &m, &d) != 3) {
        return -1;
    }

    tm.tm_year = y - 1900;
    tm.tm_mon  = m - 1;
    tm.tm_mday = d;

    for (i=0; i<NLANGS; i++) {
        strcpy(lang_strings[i].longdate, release_date);
        err = get_intl_progdate(&tm, i);
        fprintf(stderr, "lang %d: '%s', '%s'\n", i, lang_strings[i].longdate,
                lang_strings[i].shortdate);
    }

    if (!err) {
        err = print_subst(fname, release_version);
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

/* read the little "subst" file for version number and dates:
   returns non-zero on error.
*/

int get_info_from_subst_file (const char *fname, char *verstr,
                              char *progdate)
{
    FILE *fp;
    char line[MYLEN];
    int y, m, d;
    int n_read = 0;

    *verstr = '\0';
    *progdate = '\0';

    fp = fopen(fname, "r");
    if (fp == NULL) {
        return -1;
    }

    while (fgets(line, sizeof line, fp) && n_read < 2) {
        if (!strncmp(line, "VERSION", 7)) {
            if (sscanf(line, "%*s %7s", verstr)) {
                n_read++;
            }
        } else if (strstr(line, "en_US LONGDATE")) {
            if (sscanf(line + 22, "%d-%d-%d", &y, &m, &d) == 3) {
                sprintf(progdate, "%d-%02d-%02d", y, m, d);
                n_read++;
            }
        }
    }

    fclose(fp);

    return n_read == 2 ? 0 : -1;
}

int get_version_from_log (char *relver, char *reldate)
{
    char fname[MYLEN];
    char line[MYLEN];
    FILE *fp;
    int found = 0;

    sprintf(fname, "%s/ChangeLog", SRCDIR);
    fp = fopen(fname, "r");
    if (fp == NULL) {
        return 1;
    }

    while (fgets(line, sizeof line, fp)) {
        if (isdigit(line[0]) && isdigit(line[1])) {
            if (sscanf(line, "%s version %s", reldate, relver) == 2) {
                found = 1;
                break;
            }
        }
    }

    fclose(fp);

    return found ? 0 : 1;
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
    } else if (strstr(fsrc, "turk") || strstr(fsrc, "_tr")) {
        lang = TR;
    } else if (strstr(fsrc, "russi") || strstr(fsrc, "_ru")) {
        lang = RU;
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
        { "gretl_espanol_pat.html",   "es.html" },
        { "win32_pat_es.html",        "win32/index_es.html", },
        { "gretl_italiano_pat.html",  "it.html" },
        { "win32_pat_it.html",        "win32/index_it.html" },
        { "osx_pat_it.html",          "osx_it.html" },
        { "gretl_portugues_pat.html", "pt.html" },
        { "win32_pat_pt.html",        "win32/index_pt.html" },
        { "gretl_turkish_pat.html",   "tr.html" },
        { "win32_pat_tr.html",        "win32/index_tr.html" },
        { "gretl_russian_pat.html",   "ru.html" },
        { NULL, NULL }
    };
    struct from_to *ptr = templates;
    char fullsrc[MYLEN];
    char fulltarg[MYLEN];
    int err = 0;

    while (ptr->src != NULL && !err) {

        /* take source from git-controlled directory; write
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
    int year, mon, day;

    clog = fopen(fname, "r");
    if (clog == NULL) {
        fprintf(stderr, "Couldn't open '%s' for reading\n", fname);
        return 0;
    }

    while (fgets(line, MYLEN, clog)) {
        ret = sscanf(line, "%d-%d-%d version %d.%d.%d\n",
                     &year, &mon, &day, &M, &m, &r);
        if (ret == 6) {
            versions[n].major = M;
            versions[n].minor = m;
            versions[n].rev = r;
            versions[n].letter = 0;
            versions[n].pre = 0;
            n++;
        } else if (ret == 5 && year == 2002 && mon >= 9 && (mon < 11 || day <= 15)) {
            /* version 1.0 and its pre-releases */
            versions[n].major = 1;
            versions[n].minor = 0;
            versions[n].rev = 0;
            versions[n].letter = 0;
            if (mon == 11 && day == 15) {
                versions[n].pre = 0;
            } else {
                versions[n].pre = line[25];
            }
            n++;
        } else {
            char c;

            ret = sscanf(line, "%d-%d-%d version %d%c\n",
                         &year, &mon, &day, &M, &c);
            if (ret == 5 && M >= 2015) {
                versions[n].major = M;
                versions[n].minor = 0;
                versions[n].rev = 0;
                versions[n].letter = c;
                versions[n].pre = 0;
                n++;
            }
        }
    }

    fclose(clog);

    return n;
}

void bug_print_line (char *s, FILE *fp)
{
    int bnum = 0, bnum2 = 0;

    if (!strncmp(s, "version", 7) && strstr(s, ", in progress")) {
        fputs("<b>", fp);
        while (*s && *s != '\n') {
            fputc(*s, fp);
            s++;
        }
        fputs("</b>\n", fp);
    } else if (sscanf(s, "- Fix bug %d\n", &bnum) == 1 ||
        sscanf(s, "- Fix bugs %d, %d\n", &bnum, &bnum2) == 2) {
        int n;

        while (*s) {
            n = strspn(s, "0123456789");
            if (n == 7) {
                /* old-style SF bug reference */
                sscanf(s, "%d", &bnum);
                fprintf(fp, "<a href=\"http://sourceforge.net/tracker/index.php?"
                        "func=detail&aid=%d&group_id=36234&atid=416803\">%d</a>",
                        bnum, bnum);
                s += n;
            } else if (n > 0) {
                /* new-style reference */
                sscanf(s, "%d", &bnum);
                fprintf(fp, "<a href=\"http://sourceforge.net/p/gretl/bugs/%d/\">%d</a>",
                        bnum, bnum);
                s += n;
            } else {
                fputc(*s, fp);
                s++;
            }
        }
    } else {
        fputs(s, fp);
    }
}

void write_changelog (char *src, char *targ, gretl_version *gv, int nv)
{
    int ns, M, m, r, p, i;
    FILE *fin;
    FILE *fout;
    char c, line[MYLEN];
    int year, mon, day;

    fout = fopen(targ, "a");
    fin = fopen(src, "r");

    fputs("<table cellspacing=\"10\">\n", fout);

    fputs("<a name=\"top\">\n", fout);
    fputs("<tr>\n<td valign=\"top\">\n<pre>\n&nbsp;\n", fout);

    for (i=0; i<nv; i++) {
        M = gv[i].major;
        m = gv[i].minor;
        r = gv[i].rev;
        c = gv[i].letter;
        p = gv[i].pre;

        if (c != 0) {
            fprintf(fout, "<a href=\"#v%d%c\">", M, c);
        } else if (p == 0) {
            fprintf(fout, "<a href=\"#v%d-%d-%d\">", M, m, r);
        } else {
            fprintf(fout, "<a href=\"#v1pre%c\">", p);
        }

        if (i == 0) {
            if (c != 0) {
                fprintf(fout, "Version %d%c</a>&nbsp&nbsp;\n", M, c);
            } else {
                fprintf(fout, "Version %d.%d.%d</a>&nbsp&nbsp;\n", M, m, r);
            }
        } else {
            if (c != 0) {
                fprintf(fout, "Version %d%c</a>\n", M, c);
            } else if (p == 0) {
                fprintf(fout, "Version %d.%d.%d</a>\n", M, m, r);
            } else {
                fprintf(fout, "Version 1.0pre%c</a>\n", p);
            }
        }
    }

    fputs("</pre>\n</td>\n<td valign=\"top\">\n<pre>\n", fout);

    while (fgets(line, MYLEN, fin)) {
        int doparen = 1;
        int err = 0;

        ns = sscanf(line, "%d-%d-%d version %d.%d.%d\n",
                    &year, &mon, &day, &M, &m, &r);
        if (ns == 6) {
            fprintf(fout, "<b><a name=\"v%d-%d-%d\">", M, m, r);
            fprintf(fout, " %d-%02d-%02d Version %d.%d.%d</a></b>",
                    year, mon, day, M, m, r);
        } else if (ns == 5 && year == 2002 && mon >= 9 && (mon < 11 || day <= 15)) {
            /* 1.0 and its pre-releases */
            if (mon == 11 && day == 15) {
                fprintf(fout, "<b><a name=\"v1-0-0\">");
                fprintf(fout, " %d-%02d-%02d Version 1.0</a></b>",
                        year, mon, day);
            } else {
                fprintf(fout, "<b><a name=\"v1pre%c\">", line[25]);
                fprintf(fout, " %d-%02d-%02d Version 1.0pre%c</a></b>",
                        year, mon, day, line[25]);
            }
            doparen = 0;
        } else {
            ns = sscanf(line, "%d-%d-%d version %d%c\n",
                        &year, &mon, &day, &M, &c);
            if (ns == 5) {
                if (M == 0 && c == '.') {
                    /* very old versions */
                    char tail[8];

                    sscanf(line, "%d-%d-%d version %s",
                           &year, &mon, &day, tail);
                    fprintf(fout, "<b><a name=\"v%s\">", tail);
                    fprintf(fout, " %d-%02d-%02d Version %s</a></b>",
                            year, mon, day, tail);
                } else {
                    fprintf(fout, "<b><a name=\"v%d%c\">", M, c);
                    fprintf(fout, " %d-%02d-%02d Version %d%c</a></b>",
                            year, mon, day, M, c);
                }
            } else {
                bug_print_line(line, fout);
                err = 1;
            }
        }
        if (!err) {
            if (doparen) {
                char *p = strchr(line, '(');

                if (p != NULL) {
                    fprintf(fout, " %.7s", p);
                }
            }
            fputs(" [<a href=\"#top\">Back to top</a>]\n", fout);
        }
    }

    fputs("</pre>\n</td>\n</tr>\n", fout);

    fclose(fin);
    fclose(fout);
}

enum {
    CHANGE_LOG,
    BACKWARD_LOG
};

int make_html_log (int type)
{
    gretl_version gv[1000];
    char targ[MYLEN];
    char src[MYLEN];
    int nv, err = 0;

    if (type == CHANGE_LOG) {
        sprintf(targ, "%s/ChangeLog.html", WEBDIR);
        sprintf(src, "%s/template/changelog.top", WEBSRC);
    } else {
        sprintf(targ, "%s/Backward.html", WEBDIR);
        sprintf(src, "%s/template/backward.top", WEBSRC);
    }

    err = copyfile(src, targ, 0);
    if (err) return err;

    if (type == CHANGE_LOG) {
        sprintf(src, "%s/ChangeLog", SRCDIR);
    } else {
        sprintf(src, "%s/CompatLog", SRCDIR);
    }

    nv = version_history(src, gv);

    write_changelog(src, targ, gv, nv);

    sprintf(src, "%s/template/changelog.end", WEBSRC);
    err = copyfile(src, targ, 1);

    return err;
}

int make_gretldata_dtd_page (void)
{
    char line[512];
    char targ[MYLEN];
    char src[MYLEN];
    FILE *fp, *fq;
    int i;

    sprintf(src, "%s/share/data/gretldata.dtd", SRCDIR);
    sprintf(targ, "%s/gretldata.dtd.html", WEBDIR);

    fp = fopen(src, "r");
    if (fp == NULL) {
        return 1;
    }

    fq = fopen(targ, "w");
    if (fq == NULL) {
        fclose(fp);
        return 1;
    }

    fputs("<html>\n<pre>\n", fq);
    while (fgets(line, sizeof line, fp)) {
        for (i=0; line[i]; i++) {
            if (line[i] == '<') {
                fputs("&lt;", fq);
            } else if (line[i] == '>') {
                fputs("&gt;", fq);
            } else {
                fputc(line[i], fq);
            }
        }
    }
    fputs("</pre>\n</html>", fq);

    fclose(fp);
    fclose(fq);

    return 0;
}

int set_working_directories (void)
{
    char *home = getenv("HOME");

    if (home == NULL) {
        return 1;
    }

    if (!strcmp(home, "/home/cottrell")) {
        strcpy(WEBDIR, "/home/cottrell/stats/esl/website");
        strcpy(SRCDIR, "/home/cottrell/src");
        strcpy(WEBSRC, "/home/cottrell/src/doc/website");
    } else if (!strcmp(home, "/home/allin")) {
        strcpy(WEBDIR, "/home/allin/stats/esl/website");
        strcpy(SRCDIR, "/home/allin/src");
        strcpy(WEBSRC, "/home/allin/src/doc/website");
    } else {
        sprintf(WEBDIR, "%s/src/gretl/doc/website", home);
        sprintf(SRCDIR, "%s/src/gretl", home);
        sprintf(WEBSRC, "%s/src/gretl/doc/website", home);
    }

    return 0;
}

int main (int argc, char **argv)
{
    char substfile[FILENAME_MAX];
    char release_version[8] = {0};
    char release_date[12] = {0};
    int got_subst = 0;
    int err = 0;

    if (argc == 2 && !strcmp(argv[1], "--help")) {
        fprintf(stderr, "%s: You can say \"auto\" to read version info from the ChangeLog.\n",
                argv[0]);
        fputs("* Or You can specify a version string and date (YYYY-MM-DD).\n", stderr);
        fputs("* With no args, we use the existing subst file, if present.\n", stderr);
        exit(EXIT_SUCCESS);
    } else {
        fprintf(stderr, "%s:\n", argv[0]);
    }

    err = set_working_directories();
    if (err) {
        fputs("Couldn't determine working directories (no $HOME)!\n", stderr);
        exit(EXIT_FAILURE);
    }

    /* we'll want this for reading or writing */
    sprintf(substfile, "%s/subst", WEBDIR);

    if (argc == 2 && !strcmp(argv[1], "auto")) {
        err = get_version_from_log(release_version, release_date);
        if (err) {
            fputs("Couldn't get version + date from ChangeLog\n", stderr);
            exit(EXIT_FAILURE);
        } else {
            fprintf(stderr, "Got version %s, date %s from ChangeLog\n",
                    release_version, release_date);
        }
    } else if (argc == 3) {
        strncat(release_version, argv[1], 5);
        strncat(release_date, argv[2], 10);
        fprintf(stderr, "Got version %s, date %s from command line\n",
                release_version, release_date);
    } else {
        fprintf(stderr, "Trying substfile '%s'\n", substfile);
        err = get_info_from_subst_file(substfile, release_version,
                                       release_date);
        if (err) {
            fprintf(stderr, "Couldn't get release + date from %s: try 'auto'\n",
                    substfile);
            exit(EXIT_FAILURE);
        } else {
            fprintf(stderr, "Got version %s, date %s from subst file\n",
                    release_version, release_date);
            got_subst = 1;
        }
    }

    lang_strings_init();

    if (got_subst) {
        err = read_subst_file_full(substfile);
    } else {
        make_subst_file(substfile, release_version, release_date);
    }

    if (err) {
        fputs("Stopping on error!\n", stderr);
    } else {
        err = process_templates(release_version);
    }
    if (!err) {
        err = make_html_log(CHANGE_LOG);
    }
    if (!err) {
        err = make_html_log(BACKWARD_LOG);
    }
    if (!err) {
        err = make_gretldata_dtd_page();
    }

    return err;
}
