#include <gnome.h>

static void
callback_1 (GtkWidget *button, gpointer data)
{
    GtkWidget *less;

    less = gnome_less_new();
    gnome_less_show_file(GNOME_LESS(less), "myfile");
}

static void
callback_2 (GtkWidget *button, gpointer data)
{
    gnome_url_show("http://www.ecn.wfu.edu/gretl/");
}

static void
callback_3 (GtkWidget *button, gpointer data)
{
    char * const gargv[] = {"gnomecal", NULL};

    gnome_execute_async(NULL, 1, gargv);
}

static void
foo_1 (void)
     /* See libgnome/gnome-config.h */
{
    int counter;
    char *text;
    gboolean def;

    counter = gnome_config_get_int_with_default("/gretl/section/counter=1",
						&def);
    if (def) g_print("Default used for counter!\n");
    text = gnome_config_get_string("/gretl/section/text=TEXT");

    g_free(text);
}

static void
foo_2 (void)
     /* See libgnome/gnome-config.h */
{
    char *text;
    int counter;

    /*after we have set text and counter to some values we can
      write them to our config location*/
    gnome_config_set_int("/gretl/section/counter", counter);
    gnome_config_set_string("/gretl/section/text", text);
    gnome_config_sync();
}

static void help_thing = {
    GNOMEUIINFO_HELP("app_name"),
    GNOMEUIINFO_MENU_ABOUT_ITEM(callback, data),
    GNOMEUIINFO_END
};


static void
my_gnome_error (const char *filename)
{
    FILE *fp;

    fp = fopen(filename, "r");
    if (!fp) {
	char *err = g_strdup_printf(_("Cannot open file %s"),
				    filename);
	gnome_app_error(GNOME_APP(app), err);
	g_free(err);

    }
}

static void
really_delete_handler(int reply, gpointer data)
{
    GtkWidget *some_widget;

    some_widget = GTK_WIDGET(data);

    if (reply == 0) {
	/* yes was selected */
	;
    }
}

static void
my_gnome_q (void)
{
    /* ask a yes no question, passing some_widget as the data
       argument to the handler */
    gnome_app_question(GNOME_APP(app),_("Really delete object?"),
		       really_delete_handler,
		       some widget);
}

static void about_thing (void)
{
    GtkWidget* dlg;
    char *authors[] = {
	"Allin Cottrell",
	NULL;
    };

    dlg = gnome_about_new("gretl", 
			  VERSION, 
			  "(c) 2000-2001 Allin Cottrell", 
			  authors, 
			  "Just some application, "
			  "blah blah blah",  /* Other comments. */
			  "gretl-logo.xpm"
			  );

    gnome_dialog_set_parent(GNOME_DIALOG(dlg), GTK_WINDOW(app));

    gtk_widget_show(dlg);
}

int main (void) {

    GtkWidget *app;

    gnome_init ("gretl", "0.1", argc, argv);

    app = gnome_app_new ("gretl", "Econometrics program");

    return 0;

}

