/* NLS setup code for macOS in common between gretl and gretl_edit: we
   put it here so that it stays in sync across the two executables.
   Added October 2026.
*/

#define LOCALE_DEBUG 1

#ifdef LOCALE_DEBUG

#include <CoreFoundation/CoreFoundation.h>

/* Consult macOS-native language listing */

static void show_preferred_langs (void)
{
    CFArrayRef langs = CFLocaleCopyPreferredLanguages();

    if (langs) {
	char buf[256];
	CFStringRef lstr;
        CFIndex nl = CFArrayGetCount(langs);
	CFIndex i;

        for (i=0; i<nl; i++) {
            lstr = (CFStringRef) CFArrayGetValueAtIndex(langs, i);
            if (CFStringGetCString(lstr, buf, sizeof(buf), kCFStringEncodingUTF8)) {
                fprintf(stderr, "Language %ld: '%s'\n", (long) i, buf);
            }
        }
        CFRelease(langs);
    }
}

/* See what we get from CFLocale */

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

#endif /* LOCALE_DEBUG */

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

#if LOCALE_DEBUG
    macos_check_locale();
    show_preferred_langs();
#endif

    if (p != NULL) {
	fprintf(stderr, "NLS init: setlocale() gave '%s'\n", p);
        /* for the benefit of gettext */
        gretl_setenv("LANGUAGE", p);
    }
    bindtextdomain(PACKAGE, localedir);
    textdomain(PACKAGE);
    bind_textdomain_codeset(PACKAGE, "UTF-8");
}
