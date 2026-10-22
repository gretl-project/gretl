/* NLS setup code for macOS in common between gretl and gretl_edit: we
   put it here so that it stays in sync across the two executables.
*/

#include <CoreFoundation/CoreFoundation.h>

/* Let's see what we get from setlocale() */

static void macos_check_locale (void)
{
    CFLocaleRef cfloc = CFLocaleCopyCurrent();
    CFStringRef cfprop;
    const char *s;

    cfprop = (CFStringRef) CFLocaleGetValue(cfloc, kCFLocaleIdentifier);
    s = CFStringGetCStringPtr(cfprop, kCFStringEncodingASCII);
    if (s != NULL) {
	fprintf(stderr, "macos_check_locale: CF gave ID '%s'\n", s);
    }

    CFRelease(cfloc);
}

static void real_nls_init (void)
{
    char *gretlhome = getenv("GRETL_HOME");
    char localedir[MAXSTR];
    char *p;

    if (gretlhome == NULL) {
	return;
    }

    g_strlcpy(localedir, gretlhome, sizeof localedir);
    p = strstr(localedir, "share/gretl");
    if (p != NULL) {
	strcpy(p, "share/locale");
    }

    p = setlocale(LC_ALL, "");
    macos_check_locale();
    if (p != NULL) {
	fprintf(stderr, "NLS init: setlocale() gave '%s'\n", p);
        /* for the benefit of gettext */
        gretl_setenv("LANGUAGE", p);
    }
    bindtextdomain(PACKAGE, localedir);
    textdomain(PACKAGE);
    bind_textdomain_codeset(PACKAGE, "UTF-8");
}
