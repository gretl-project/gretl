/*
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* osx_open.c for gretl -- open various sorts of files */

#include "gretl.h"
#include "osx_open.h"

#ifdef USE_CARBON

# include <Carbon/Carbon.h>

/* deprecated in macOS >= 10.10, removed in macOS 11 */

int osx_open_file (const char *path)
{
    FSRef ref;
    int err;

    err = FSPathMakeRef((const UInt8 *) path, &ref, NULL);
    if (!err) {
        err = LSOpenFSRef(&ref, NULL);
    }

    return err;
}

int osx_open_pdf (const char *path, const char *dest)
{
    FSRef ref;
    int done = 0;
    int err;

    err = FSPathMakeRef((const UInt8 *) path, &ref, NULL);

    if (!err) {
        guint8 exe[PATH_MAX] = {0};
        FSRef appref;

        err = LSGetApplicationForItem(&ref, kLSRolesViewer | kLSRolesEditor,
                                      &appref, NULL);

        if (!err) {
            FSRefMakePath(&appref, exe, PATH_MAX);
        }

        if (!err && strstr((const char *) exe, "Adobe") != NULL) {
            /* Adobe Acrobat or Acrobat Reader: try passing an
               option to open at the specified chapter.
            */
            LSLaunchFSRefSpec rspec;
            AEDesc desc;
            gchar *opt;
            int lserr;

            opt = g_strdup_printf("nameddest=%s", dest);
            AECreateDesc(typeChar, opt, strlen(opt), &desc);

            rspec.appRef = &appref;
            rspec.numDocs = 1;
            rspec.itemRefs = &ref;
            rspec.passThruParams = &desc;
            rspec.launchFlags = kLSLaunchAsync;
            rspec.asyncRefCon = NULL;

            lserr = LSOpenFromRefSpec(&rspec, NULL);
            if (lserr) {
                fprintf(stderr, "LSOpenFromRefSpec, err = %d\n", lserr);
            } else {
                done = 1;
            }

            AEDisposeDesc(&desc);
            g_free(opt);
        } else if (!err && strstr((const char *) exe, "Preview") != NULL) {
            /* The default Apple Preview.app: there's no option
               as per Adobe, but we can at least try to open the
               Table-of-Contents pane (Option-Control-3).
            */
            err = LSOpenFSRef(&ref, NULL);
            if (!err) {
                FILE *fp = popen("/usr/bin/osascript", "w");

                if (fp != NULL) {
                    /* try to get the table of contents shown */
                    fputs("activate application \"Preview\"\n", fp);
                    fputs("tell application \"System Events\"\n", fp);
                    fputs(" keystroke \"3\" using {option down, command down}\n", fp);
                    fputs("end tell\n", fp);
                    pclose(fp);
                }
                done = 1;
            }
        }
    }

    if (!err && !done) {
        err = LSOpenFSRef(&ref, NULL);
    }

    return err;
}

#else /* macOS >= 10.10, no Carbon */

# include <CoreFoundation/CoreFoundation.h>
# include <CoreServices/CoreServices.h>

int osx_open_file (const char *path)
{
    CFURLRef ref;
    int err = 0;

    ref = CFURLCreateFromFileSystemRepresentation(NULL,
                                                  (const UInt8 *) path,
                                                  strlen(path),
                                                  false);
    if (ref != NULL) {
        err = LSOpenCFURLRef(ref, NULL);
        CFRelease(ref);
    } else {
        err = 1;
    }

    return err;
}

int osx_open_pdf (const char *path, const char *dest)
{
    CFURLRef ref;
    CFURLRef appref;
    int viewer = 0;
    int done = 0;
    int err = 0;

    ref = CFURLCreateFromFileSystemRepresentation(NULL,
						  (const UInt8 *) path,
						  strlen(path),
						  false);
    if (ref == NULL) {
        return 1;
    }

    appref = LSCopyDefaultApplicationURLForURL(ref, kLSRolesAll, NULL);

    if (appref == NULL) {
        err = 1;
    } else {
        CFStringRef exe = CFURLGetString(appref);
        const char *s[] = {"Adobe", "Preview"};
        CFStringRef v[2];
        CFRange cfr;
        int i;

        v[0] = CFStringCreateWithCString(NULL, s[0], kCFStringEncodingASCII);
        v[1] = CFStringCreateWithCString(NULL, s[1], kCFStringEncodingASCII);

        for (i=0; i<2; i++) {
            cfr = CFStringFind(exe, v[i], 0);
            if (cfr.length > 0) {
                viewer = i+1;
            }
        }

        CFRelease(v[0]);
        CFRelease(v[1]);
    }

    if (!err && viewer == 1) {
        /* Adobe Acrobat or Acrobat Reader: try passing an
           option to open at the specified chapter.
        */
        const void *vals = {ref};
        CFArrayRef refs;
        LSLaunchURLSpec uspec;
        AEDesc desc;
        gchar *opt;
        int lserr;

        opt = g_strdup_printf("nameddest=%s", dest);
        AECreateDesc(typeChar, opt, strlen(opt), &desc);
        refs = CFArrayCreate(NULL, &vals, 1, &kCFTypeArrayCallBacks);

        uspec.appURL = appref;
        uspec.itemURLs = refs;
        uspec.passThruParams = &desc;
        uspec.launchFlags = kLSLaunchAsync;
        uspec.asyncRefCon = NULL;

        lserr = LSOpenFromURLSpec(&uspec, NULL);
        if (lserr) {
            fprintf(stderr, "LSOpenFromURLSpec, err = %d\n", lserr);
        } else {
            done = 1;
        }
        AEDisposeDesc(&desc);
        g_free(opt);
    } else if (!err && viewer == 2) {
        /* The default Apple Preview.app: there's no option
           as per Adobe, but we can at least try to open the
           Table-of-Contents pane (Option-Control-3).
        */
        err = LSOpenCFURLRef(ref, NULL);
        if (!err) {
            FILE *fp = popen("/usr/bin/osascript", "w");

            if (fp != NULL) {
                /* try to get the table of contents shown */
                fputs("activate application \"Preview\"\n", fp);
                fputs("tell application \"System Events\"\n", fp);
                fputs(" keystroke \"3\" using {option down, command down}\n", fp);
                fputs("end tell\n", fp);
                pclose(fp);
            }
            done = 1;
        }
    }

    if (!done) {
        /* try just opening the document */
        err = LSOpenCFURLRef(ref, NULL);
    }
    CFRelease(ref);

    return err;
}

# endif /* Carbon-free variant */

int osx_open_url (const char *url)
{
    CFStringRef s;
    CFURLRef u;
    int err;

    s = CFStringCreateWithBytes(NULL, (const UInt8 *) url, strlen(url),
                                kCFStringEncodingASCII,
                                0);
    if (s == NULL) {
        err = 1;
    } else {
        u = CFURLCreateWithString(NULL, s, NULL);
        if (u == NULL) {
            err = 1;
        } else {
            err = LSOpenCFURLRef(u, NULL);
            CFRelease(u);
        }
        CFRelease(s);
    }

    return err;
}
