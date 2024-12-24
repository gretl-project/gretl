#include <sys/param.h>
#include <mach-o/dyld.h>

/* This function provides an replacement for the old method
   of starting gretl on macOS via a shell script: instead, we
   get the environment set up correctly via C code.
*/

void osx_setup_paths (void)
{
    char userpath[MAXPATHLEN + 1];
    char execpath[MAXPATHLEN + 1];
    const char *gpshare = "6.0";
    uint32_t pathsz = sizeof execpath;
    gchar *tmp, *respath = NULL;
    char *c, *rhome;

    /* record the initial working directory */
    getcwd(userpath, sizeof userpath);
    setenv("GRETL_STARTDIR", userpath, 1);

    /* get the full path to the active executable */
    _NSGetExecutablePath(execpath, &pathsz);

    /* The gretl executable lives in Contents/MacOS under Gretl.app.
       We need to obtain the Contents/Resources path.
    */
    c = strstr(execpath, "/Contents/MacOS");
    *c = '\0';
    respath = g_strdup_printf("%s/Contents/Resources", execpath);

    tmp = g_strdup_printf("%s/share/gretl/", respath);
    setenv("GRETL_HOME", tmp, 1);
    g_free(tmp);

    setenv("GTK_PATH", respath, 1);
    setenv("GTK_DATA_PREFIX", respath, 1);
    setenv("GTK_EXE_PREFIX", respath, 1);
    setenv("G_FILENAME_ENCODING", "UTF-8", 1);

#if GTK_MAJOR_VERSION == 3
    tmp = g_strdup_printf("%s/share/glib-2.0/schemas", respath);
    setenv("GSETTINGS_SCHEMA_DIR", tmp, 1);
    g_free(tmp);
#else /* GTK2 */
    tmp = g_strdup_printf("%s/etc/gtk-2.0/gtk.immodules", respath);
    setenv("GTK_IM_MODULE_FILE", tmp, 1);
    g_free(tmp);
    tmp = g_strdup_printf("%s/etc/gtk-2.0/gtkrc", respath);
    setenv("GTK2_RC_FILES", tmp, 1);
    g_free(tmp);
    /* pango stuff */
    tmp = g_strdup_printf("%s/etc", respath);
    setenv("PANGO_SYSCONFDIR", tmp, 1);
    g_free(tmp);
    tmp = g_strdup_printf("%s/etc/pango/pangorc", respath);
    setenv("PANGO_RC_FILE", tmp, 1);
    g_free(tmp);
    tmp = g_strdup_printf("%s/lib", respath);
    setenv("PANGO_LIBDIR", tmp, 1);
    g_free(tmp);
#endif

    /* gnuplot variables */
    tmp = g_strdup_printf("%s/share/gnuplot/%s/gnuplot.gih", respath, gpshare);
    setenv("GNUHELP", tmp, 1);
    g_free(tmp);
    tmp = g_strdup_printf("%s/share/gnuplot/%s/PostScript", respath, gpshare);
    setenv("GNUPLOT_PS_DIR", tmp, 1);
    g_free(tmp);
    setenv("GNUTERM", "wxt", 1);

    /* R shared library: may require platform-specific symlink */
    rhome = getenv("R_HOME");
    if (rhome == NULL) {
        setenv("R_HOME", "/Library/Frameworks/R.framework/Resources", 1);
    }

    /* Add ourselves to PATH? */
    if (1) {
        char *path0 = getenv("PATH");

        if (path0 != NULL && *path0 != '\0') {
            tmp = g_strdup_printf("%s:%s/bin", path0, respath);
        } else {
            tmp = g_strdup_printf("%s/bin", respath);
        }
        setenv("PATH", tmp, 1);
        g_free(tmp);
    }

    g_free(respath);
}
