#include <sys/param.h>
#include <mach-o/dyld.h>

/* This function provides an replacement for the old method
   of starting gretl on macOS via a shell script: instead, we
   get the environment set up correctly via C code.
*/

void osx_setup_paths (void)
{
    char *c, execpath[MAXPATHLEN + 1];
    char respath[MAXPATHLEN + 1];
    uint32_t pathsz = sizeof execpath;
    char *rhome;
    gchar *tmp;

    getcwd(respath, sizeof respath);
    setenv("GRETL_STARTDIR", respath, 1);

    _NSGetExecutablePath(execpath, &pathsz);

    c = strrchr(execpath, '/');
    *c = '\0';

    strcat(execpath, "/../Resources");
    chdir(execpath);

    getcwd(respath, sizeof respath);

    tmp = g_strdup_printf("%s/share/gretl/", respath);
    setenv("GRETL_HOME", tmp, 1);
    g_free(tmp);

    setenv("GTK_PATH", respath, 1);
    setenv("GTK_DATA_PREFIX", respath, 1);
    setenv("GTK_EXE_PREFIX", respath, 1);
    setenv("G_FILENAME_ENCODING", "UTF-8", 1);

    tmp = g_strdup_printf("%s/etc/gtk-2.0/gtkrc", respath);
    setenv("GTK2_RC_FILES", tmp, 1);
    g_free(tmp);
    tmp = g_strdup_printf("%s/etc/gtk-2.0/gtk.immodules", respath);
    setenv("GTK_IM_MODULE_FILE", tmp, 1);
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

    /* gnuplot variables */
    tmp = g_strdup_printf("%s/share/gnuplot/5.2/gnuplot.gih", respath);
    setenv("GNUHELP", tmp, 1);
    g_free(tmp);
    tmp = g_strdup_printf("%s/share/gnuplot/5.2/PostScript", respath);
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
}
