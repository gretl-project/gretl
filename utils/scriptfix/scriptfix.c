/* Checker-and-fixer utility to process a gretl gfn file (or a plain
   inp file though that's not much tested), looking for occurrences
   of a single '=' where it appears that "==" is intended.

   Given an argument of foo.gfn, if there are any errors found,
   the program

   * writes foo.gfn.fixed (it doesn't touch the original input); and
   * writes to stdout an account of the changes made, with their
     associated line numbers.

   If no errors are found, that's reported on stdout and the "fixed"
   file is deleted on exit.

   If the --verbose flag is appended after the filename argument,
   then you get a verbose account of the checking process on stderr.

   No special libraries needed; compile with (e.g.)

   gcc -O2 -Wall -o fixeq fixeq.c

   Allin Cottrell, 2017-02-12
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define LLEN 8192 /* should be generous for one line? */

int verbose = 0;

void verbose_print (const char *s)
{
    if (verbose) {
	fputs(s, stderr);
    }
}

int is_operator (char c)
{
    const char *ops = "!.+-%^*/<>";
    int i;

    for (i=0; ops[i]; i++) {
	if (c == ops[i]) {
	    return 1;
	}
    }

    return 0;
}

int maybe_fix_eq (char *line, char *s, int xml)
{
    int offset = s - line;
    int err = 1;

    if (offset > 0 && is_operator(*(s-1))) {
	/* '!=', '.=' or similar, OK; or '>=', '<=' in
	   plain inp file
	*/
	err = 0;
    } else if (xml && offset > 4) {
	s -= 4;
	if (!strncmp(s, "&lt;", 4) || !strncmp(s, "&gt;", 4)) {
	    /* '<=' or '>=', OK */
	    err = 0;
	} else {
	    s += 4;
	}
    }

    if (err) {
	char fixup[LLEN];

	if (verbose) {
	    fprintf(stderr, "*** dodgy line? %s", line);
	}
	*fixup = '\0';
	strncat(fixup, s, strlen(s));
	s[0] = '=';
	s[1] = '\0';
	strncat(s+1, fixup, strlen(fixup));
    }

    return err;
}

int is_option (const char *s)
{
    if (s[0] == '-' && s[1] == '-' && isalpha(s[2])) {
	const char *optchars = "abcdefghijklmnopqrstuvwxyz-";

	s += 2;
	s += strspn(s, optchars);
	if (*s == '=' || isspace(*s)) {
	    return 1;
	}
    }

    return 0;
}

void get_pkg_info (char *targ, char *src)
{
    int n = strcspn(src, "<");

    *targ = '\0';
    strncat(targ, src, n);
}

char *xml_unescape (const char *s)
{
    char *ret = calloc(strlen(s) + 1, 1);
    int i = 0;

    while (*s) {
	if (*s == '&') {
	    s++;
	    if (!strncmp(s, "amp;", 4)) {
		ret[i++] = '&';
		s += 3;
	    } else if (!strncmp(s, "quot;", 5)) {
		ret[i++] = '"';
		s += 4;
	    } else if (!strncmp(s, "lt;", 3)) {
		ret[i++] = '<';
		s += 2;
	    } else if (!strncmp(s, "gt;", 3)) {
		ret[i++] = '>';
		s += 2;
	    } else {
		fprintf(stderr, "?? xml_unescape: %s", s - 1);
	    }
	} else {
	    ret[i++] = *s;
	}
	s++;
    }

    return ret;
}

void report_fix (const char *orig, const char *line, int i,
		 const char *funcname, int xml)
{
    if (*funcname != '\0') {
	printf("line %d, in function %s\n", i, funcname);
    } else {
	printf("line %d\n", i);
    }
    
    if (xml && strchr(orig, '&') != NULL) {
	char *tmp;

	tmp = xml_unescape(orig);
	printf(" original: %s", tmp);
	free(tmp);
	tmp = xml_unescape(line);
	printf("  revised: %s", tmp);
	free(tmp);
    } else {
	printf(" original: %s", orig);
	printf("  revised: %s", line);
    }
}

char *get_fixname (const char *fname)
{
    char *ret = malloc(strlen(fname) + 7);

    sprintf(ret, "%s.fixed", fname);
    return ret;
}

/* "type = name(args)" */

void set_funcname (char *targ, const char *src)
{
    /* skip return type */
    src += strspn(src, " \t");
    src += strcspn(src, " \t");
    src += strspn(src, " \t");
    /* and scan name */
    sscanf(src, "%31[^ \t(]", targ);
    if (verbose) {
	fprintf(stderr, "entering function %s\n", targ);
    }
}

/* name="fname" */

void xml_set_funcname (char *targ, const char *src)
{
    /* skip 'name="' */
    src += 6;
    /* and scan name */
    sscanf(src, "%31[^\"]", targ);
    if (verbose) {
	fprintf(stderr, "entering function %s\n", targ);
    }
}

void clear_funcname (char *name)
{
    if (verbose) {
	fprintf(stderr, "exiting function %s\n", name);
    }
    *name = '\0';
}

int has_suffix (const char *str, const char *sfx)
{
    const char *p;
    int ret = 0;

    if (str != NULL && sfx != NULL) {
	p = strrchr(str, *sfx);
	if (p != NULL && strlen(p) == strlen(sfx)) {
	    ret = 1;
	    while (*p) {
		if (*p != *sfx && *p != toupper(*sfx)) {
		    ret = 0;
		    break;
		}
		p++;
		sfx++;
	    }
	}
    }

    return ret;
}

const char *pkgname (const char *s)
{
    const char *p = strrchr(s, '/');

    if (p != NULL) {
	s = p + 1;
    }

    return s;
}

int is_quote (const char *s, int xml)
{
    if (xml) {
	return strncmp(s, "&quot;", 6) == 0;
    } else {
	return *s == '"' && (*s-1) != '\\';
    }
}

int main (int argc, char **argv)
{
    char *fname, *fixname = NULL;
    char line[LLEN], orig[LLEN];
    char funcname[32];
    char pkgver[16];
    char pkgdate[16];
    char *s, *p;
    int xml = 1;
    int in_code = 0;
    int in_comment = 0;
    int in_foreign = 0;
    int quoted = 0;
    int info_done = 0;
    int total_fixes = 0;
    int i = 0;
    FILE *fp, *fq;

    if (argc < 2) {
	fprintf(stderr, "%s: we need the name of a gfn or inp file\n", argv[0]);
	exit(EXIT_FAILURE);
    }

    fname = argv[1];
    if (has_suffix(fname, ".inp")) {
	xml = 0;
	in_code = 1;
    }

    if (argc == 3 && !strcmp(argv[2], "--verbose")) {
	verbose = 1;
    }

    fp = fopen(fname, "rb");
    if (fp == NULL) {
	fprintf(stderr, "%s: couldn't open %s\n", argv[0], fname);
	exit(EXIT_FAILURE);
    }

    fixname = get_fixname(fname);
    fq = fopen(fixname, "w");
    if (fq == NULL) {
	fprintf(stderr, "%s: couldn't write to %s\n", argv[0], fixname);
    }

    if (!xml) {
	printf("*** %s ***\n", pkgname(fname));
    }

    *pkgver = *pkgdate = *funcname = '\0';

    while (fgets(line, sizeof line, fp)) {
	int sample_start = 0;
	int code_start = 0;
	int fixups = 0;
	int rhs = 0;

	strcpy(orig, line); /* for reference */

	i++;

	if (xml) {
	    if (!strncmp(line, "<version>", 9)) {
		get_pkg_info(pkgver, line + 9);
	    } else if (!strncmp(line, "<date>", 6)) {
		get_pkg_info(pkgdate, line + 6);
	    } else if (!strncmp(line, "<code>", 6)) {
		verbose_print("entering code\n");
		code_start = in_code = 1;
	    } else if (!strncmp(line, "</code>", 7)) {
		verbose_print("exiting code\n");
		in_code = 0;
	    } else if (!strncmp(line, "<sample-script>", 15)) {
		verbose_print("entering sample script\n");
		sample_start = in_code = 1;
	    } else if (!strncmp(line, "</sample-script>", 16)) {
		verbose_print("exiting sample script\n");
		in_code = 0;
	    } else if (!strncmp(line, "<gretl-function ", 16)) {
		xml_set_funcname(funcname, line + 16);
	    } else if (!strncmp(line, "</gretl-function>", 17)) {
		clear_funcname(funcname);
	    }
	    if (!info_done && *pkgver && *pkgdate) {
		printf("*** %s, version %s, %s ***\n", pkgname(fname),
		       pkgver, pkgdate);
		info_done = 1;
	    }
	}

	if (in_code) {
	    s = line;
	    if (code_start) {
		s += 6;
	    } else if (sample_start) {
		s += 15;
	    }
	    if (verbose) {
		fprintf(stderr, "> %s", s);
	    }
	    s += strspn(s, " \t");
	    if (!strncmp(s, "catch ", 6)) {
		s += 6;
		s += strspn(s, " \t");
	    }
	    p = s;
	    while (*s) {
		if (!in_comment && !quoted && !in_foreign) {
		    if (!rhs && !strncmp(p, "foreign ", 8)) {
			verbose_print(" entering foreign\n");
			in_foreign = 1;
			break;
		    } else if (!xml && !rhs && !strncmp(p, "function ", 9)) {
			set_funcname(funcname, p + 9);
			break;
		    } else if (!xml && !rhs && !strncmp(p, "end function", 12)) {
			clear_funcname(funcname);
			*funcname = '\0';
			break;
		    } else if (!rhs && !strncmp(p, "if ", 3)) {
			verbose_print(" entering 'if'\n");
			rhs = 1;
			s += 2;
		    } else if (!rhs && !strncmp(p, "elif ", 5)) {
			verbose_print(" entering 'elif'\n");
			rhs = 1;
			s += 4;
		    } else if (!rhs && !strncmp(p, "smpl ", 5)) {
			verbose_print(" entering 'smpl'\n");
			rhs = 1;
			s += 4;
		    } else if (!strncmp(s, "/*", 2)) {
			verbose_print(" entering comment\n");
			in_comment = 1;
			s++;
		    } else if (is_quote(s, xml)) {
			verbose_print(" entering quote\n");
			quoted = 1;
			if (xml) s += 5; /* &quot; */
		    } else if (*s == '#') {
			verbose_print(" end-of-line comment\n");
			break;
		    } else if (is_option(s)) {
			verbose_print(" reached option flags\n");
			break;
		    } else if (*s == '=') {
			if (*(s+1) == '=') {
			    s++; /* OK */
			} else if (rhs) {
			    if (maybe_fix_eq(line, s, xml)) {
				fixups++;
				s++;
			    }
			}
			rhs = 1;
		    }
		} else if (quoted) {
		    if (is_quote(s, xml)) {
			verbose_print(" exiting quote\n");
			quoted = 0;
			if (xml) s += 5; /* &quot; */
		    }
		} else if (in_comment) {
		    if (!strncmp(s, "*/", 2)) {
			verbose_print(" exiting comment\n");
			in_comment = 0;
			s++;
		    }
		} else if (in_foreign) {
		    if (!strncmp(p, "end foreign", 11)) {
			verbose_print(" exiting foreign\n");
			in_foreign = 0;
			break;
		    }
		}
		s++;
	    }
	    if (fixups > 0) {
		report_fix(orig, line, i, funcname, xml);
	    }
	}
	total_fixes += fixups;
	fputs(line, fq);
    }

    fclose(fp);
    fclose(fq);

    if (total_fixes == 0) {
	printf("%s: no fixes required\n", fname);
	remove(fixname);
    } else {
	printf("Wrote %s (total fixes: %d)\n", fixname, total_fixes);
    }

    free(fixname);

    return 0;
}
