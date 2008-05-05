/*  gtksourcelanguagesmanager.c
 *
 *  Copyright (C) 2003 - Paolo Maggi <paolo.maggi@polito.it>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>

#include <libxml/xmlreader.h>

#include "gtksourcelanguagesmanager.h"
#include "gtksourcelanguage-private.h"
#include "gtksourcelanguage.h"
#include "gtksourcetag.h"

enum {
	PROP_0,
	PROP_LANG_SPECS_DIRS
};

struct _GtkSourceLanguagesManagerPrivate {

	GSList 		*available_languages;

	GSList		*language_specs_directories;
};

static GObjectClass *parent_class = NULL;

static void  gtk_source_languages_manager_class_init		(GtkSourceLanguagesManagerClass *klass);
static void  gtk_source_languages_manager_instance_init	(GtkSourceLanguagesManager *lm);
static void	 gtk_source_languages_manager_finalize	 	(GObject 		   *object);

static void	 slist_deep_free 				(GSList 		   *list);
static GSList 	*get_lang_files 				(GtkSourceLanguagesManager *lm);
static GSList	*build_file_listing 				(const gchar 		   *directory, 
					 			 GSList			   *filenames);

static void	 gtk_source_languages_manager_set_property 	(GObject 		   *object, 
					   			 guint 	 		    prop_id,
			    		   			 const GValue 		   *value, 
					   			 GParamSpec		   *pspec);
static void	 gtk_source_languages_manager_get_property 	(GObject 		   *object, 
					   			 guint 	 		    prop_id,
			    		   			 GValue 		   *value, 
					   			 GParamSpec		   *pspec);
static void	 gtk_source_languages_manager_set_specs_dirs	(GtkSourceLanguagesManager *lm,
								 const GSList 		   *dirs);

GType
gtk_source_languages_manager_get_type (void)
{
	static GType languages_manager_type = 0;

  	if (languages_manager_type == 0)
    	{
      		static const GTypeInfo our_info =
      		{
        		sizeof (GtkSourceLanguagesManagerClass),
        		NULL,		/* base_init */
        		NULL,		/* base_finalize */
        		(GClassInitFunc) gtk_source_languages_manager_class_init,
        		NULL,           /* class_finalize */
        		NULL,           /* class_data */
        		sizeof (GtkSourceLanguagesManager),
        		0,              /* n_preallocs */
        		(GInstanceInitFunc) gtk_source_languages_manager_instance_init
      		};

      		languages_manager_type = g_type_register_static (G_TYPE_OBJECT,
                					    "GtkSourceLanguagesManager",
							    &our_info,
							    0);
    	}

	return languages_manager_type;
}

static void
gtk_source_languages_manager_class_init (GtkSourceLanguagesManagerClass *klass)
{
	GObjectClass *object_class = G_OBJECT_CLASS (klass);

	parent_class		= g_type_class_peek_parent (klass);
	object_class->finalize	= gtk_source_languages_manager_finalize;

	object_class->set_property = gtk_source_languages_manager_set_property;
	object_class->get_property = gtk_source_languages_manager_get_property;
	
	g_object_class_install_property (object_class,
					 PROP_LANG_SPECS_DIRS,
					 g_param_spec_pointer ("lang_files_dirs",
						 	       "Language specification directories",
							       "List of directories where the "
								 "language specification files (.lang) "
								 "are located",
							       (G_PARAM_READWRITE | 
							        G_PARAM_CONSTRUCT_ONLY)));
}

static void
gtk_source_languages_manager_set_property (GObject 	*object, 
					   guint 	 prop_id,
			    		   const GValue *value, 
					   GParamSpec	*pspec)
{
	GtkSourceLanguagesManager *lm;

	g_return_if_fail (GTK_IS_SOURCE_LANGUAGES_MANAGER (object));

	lm = GTK_SOURCE_LANGUAGES_MANAGER (object);
    
	switch (prop_id) {	
	    case PROP_LANG_SPECS_DIRS:
		gtk_source_languages_manager_set_specs_dirs (
				lm, 
				g_value_get_pointer (value));
		break;

	    default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
gtk_source_languages_manager_get_property (GObject 	*object, 
					   guint 	 prop_id,
			    		   GValue 	*value, 
					   GParamSpec	*pspec)
{
	GtkSourceLanguagesManager *lm;

	g_return_if_fail (GTK_IS_SOURCE_LANGUAGES_MANAGER (object));

	lm = GTK_SOURCE_LANGUAGES_MANAGER (object);
    
	switch (prop_id) {	
	    case PROP_LANG_SPECS_DIRS:
		    g_value_set_pointer (value, lm->priv->language_specs_directories);
		    break;
	    default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
gtk_source_languages_manager_instance_init (GtkSourceLanguagesManager *lm)
{
	lm->priv = g_new0 (GtkSourceLanguagesManagerPrivate, 1);	
}

/**
 * gtk_source_languages_manager_new:
 *
 * Creates a new language manager.
 *
 * Return value: a #GtkSourceLanguagesManager.
 **/
GtkSourceLanguagesManager *
gtk_source_languages_manager_new (void)
{
	return GTK_SOURCE_LANGUAGES_MANAGER (
					g_object_new (GTK_TYPE_SOURCE_LANGUAGES_MANAGER,
				      	NULL));
}

static void
gtk_source_languages_manager_finalize (GObject *object)
{
	GtkSourceLanguagesManager *lm;
       
	lm = GTK_SOURCE_LANGUAGES_MANAGER (object);
	
	if (lm->priv->available_languages != NULL)
	{
		GSList *list = lm->priv->available_languages;
		
		g_slist_foreach (list, (GFunc) g_object_unref, NULL);
		g_slist_free (list);
	}

	slist_deep_free (lm->priv->language_specs_directories);

	g_free (lm->priv);

	(* G_OBJECT_CLASS (parent_class)->finalize) (object);
}

#ifdef G_OS_WIN32

#include <windows.h> /* read registry */

# define MAXLEN 512

static int read_gretldir_reg_val (HKEY tree, char *keyval)
{
    unsigned long datalen = MAXLEN;
    int error = 0;
    HKEY regkey;

    if (RegOpenKeyEx(
                     tree,                        /* handle to open key */
                     "Software\\gretl",           /* subkey name */
                     0,                           /* reserved */
                     KEY_READ,                    /* access mask */
                     &regkey                      /* key handle */
                     ) != ERROR_SUCCESS) {
        return 1;
    }

    if (RegQueryValueEx(
                        regkey,
                        "gretldir",
                        NULL,
                        NULL,
                        keyval,
                        &datalen
                        ) != ERROR_SUCCESS) {
        error = 1;
    }

    RegCloseKey(regkey);

    return error;
}

static void get_gretl_lang_dir (char *targ)
{
    size_t n;

    read_gretldir_reg_val(HKEY_LOCAL_MACHINE, targ);

    if (*targ == '\0') {
        strcpy(targ, "c:\\Program Files\\gretl");
    }

    n = strlen(targ);
    if (targ[n-1] != '\\') strcat(targ, "\\");

    strcat(targ, "share\\gtksourceview-1.0\\language-specs");
}

static void 
gtk_source_languages_manager_set_specs_dirs (GtkSourceLanguagesManager	*lm,
					     const GSList		*dirs)
{
    char gretl_lang_dir[MAXLEN];

    g_return_if_fail (GTK_IS_SOURCE_LANGUAGES_MANAGER (lm));
    g_return_if_fail (lm->priv->language_specs_directories == NULL);
			
    if (dirs == NULL) {
	get_gretl_lang_dir(gretl_lang_dir);
	lm->priv->language_specs_directories =
	    g_slist_prepend (lm->priv->language_specs_directories,
			     g_strdup (gretl_lang_dir));
	return;
    }

    while (dirs != NULL) {
	lm->priv->language_specs_directories = 
	    g_slist_prepend (lm->priv->language_specs_directories,
			     g_strdup ((const gchar*)dirs->data));
	dirs = g_slist_next (dirs);
    }
}

#else /* !G_OS_WIN32 */

#define DEFAULT_LANGUAGE_DIR	DATADIR "/gtksourceview"
#define USER_LANGUAGE_DIR	"gtksourceview-1.0/language-specs"
#define USER_CONFIG_BASE_DIR	".gnome2"

static void
gtk_source_languages_manager_set_specs_dirs (GtkSourceLanguagesManager *lm,
					     const GSList              *dirs)
{
	g_return_if_fail (GTK_IS_SOURCE_LANGUAGES_MANAGER (lm));
	g_return_if_fail (lm->priv->language_specs_directories == NULL);

	if (dirs == NULL)
	{
	        gchar *env_language_dir = getenv("GTKSOURCEVIEW_LANGUAGE_DIR");

		lm->priv->language_specs_directories =
			g_slist_prepend (lm->priv->language_specs_directories,
					g_strdup (DEFAULT_LANGUAGE_DIR));
		lm->priv->language_specs_directories = 
			g_slist_prepend (lm->priv->language_specs_directories,
					g_build_filename (g_get_home_dir(), 
						USER_CONFIG_BASE_DIR, USER_LANGUAGE_DIR, 
						NULL));

                if (env_language_dir != NULL && strcmp(DEFAULT_LANGUAGE_DIR, env_language_dir)) 
                {
                        lm->priv->language_specs_directories = 
                                g_slist_prepend (lm->priv->language_specs_directories,
                                        g_strdup (env_language_dir));
                }


		return;
	}
	
	while (dirs != NULL)
	{
		lm->priv->language_specs_directories = 
			g_slist_prepend (lm->priv->language_specs_directories,
					 g_strdup ((const gchar*)dirs->data));

		dirs = g_slist_next (dirs);
	}
}

#endif /* G_OS_WIN32 vs other */

/**
 * gtk_source_languages_manager_get_lang_files_dirs:
 * @lm: a #GtkSourceLanguagesManager.
 * 
 * Gets a list of language files directories for the given language manager.
 *
 * Return value: a list of language files directories (as strings).
 **/
const GSList *
gtk_source_languages_manager_get_lang_files_dirs (GtkSourceLanguagesManager *lm)
{
	g_return_val_if_fail (GTK_IS_SOURCE_LANGUAGES_MANAGER (lm), NULL);

	return lm->priv->language_specs_directories;
}

static void
prepend_lang (gchar *id, GtkSourceLanguage *lang, GtkSourceLanguagesManager *lm)
{
	lm->priv->available_languages = 
		g_slist_prepend (lm->priv->available_languages, lang);
}

/**
 * gtk_source_languages_manager_get_available_languages:
 * @lm: a #GtkSourceLanguagesManager.
 * 
 * Gets a list of available languages for the given language manager.
 * This function returns a pointer to a internal list, so there is no need to
 * free it after usage.
 *
 * Return value: a list of #GtkSourceLanguage.
 **/
const GSList *
gtk_source_languages_manager_get_available_languages (GtkSourceLanguagesManager *lm)
{
	GSList *filenames, *l;
	GHashTable *lang_hash;

	g_return_val_if_fail (GTK_IS_SOURCE_LANGUAGES_MANAGER (lm), NULL);

	if (lm->priv->available_languages != NULL)
	{
		return lm->priv->available_languages;
	}
	
	/* Build list of availables languages */
	filenames = get_lang_files (lm);
	
	lang_hash = g_hash_table_new (g_str_hash, g_str_equal);

	for (l = filenames; l != NULL; l = l->next)
	{
		GtkSourceLanguage *lang;

		lang = _gtk_source_language_new_from_file ((const gchar*)l->data,
							   lm);

		if (lang == NULL)
		{
			g_warning ("Error reading language specification file '%s'", 
				   (const gchar*)l->data);
			continue;
		}

		if (g_hash_table_lookup (lang_hash, lang->priv->id) == NULL)
		{	
			g_hash_table_insert (lang_hash, 
					     lang->priv->id,
					     lang);
		}
	}

	slist_deep_free (filenames);

	g_hash_table_foreach (lang_hash, (GHFunc) prepend_lang, lm);

	g_hash_table_destroy (lang_hash);

	return lm->priv->available_languages;
}

/**
 * gtk_source_languages_manager_get_language_from_mime_type:
 * @lm: a #GtkSourceLanguagesManager.
 * @mime_type: a mime type.
 * 
 * Gets the #GtkSourceLanguage which is associated with the given @mime_type
 * in the language manager.
 *
 * Return value: a #GtkSourceLanguage, or %NULL if there is no language
 * associated with the given @mime_type.
 **/
GtkSourceLanguage *
gtk_source_languages_manager_get_language_from_mime_type (GtkSourceLanguagesManager 	*lm,
							  const gchar 			*mime_type)
{
	const GSList *languages;
	g_return_val_if_fail (mime_type != NULL, NULL);

	languages = gtk_source_languages_manager_get_available_languages (lm);

	while (languages != NULL)
	{
		GSList *mime_types, *tmp;

		GtkSourceLanguage *lang = GTK_SOURCE_LANGUAGE (languages->data);
		
		tmp = mime_types = gtk_source_language_get_mime_types (lang);

		while (tmp != NULL)
		{
			if (strcmp ((const gchar*)tmp->data, mime_type) == 0)
			{		
				break;
			}

			tmp = g_slist_next (tmp);
		}

		slist_deep_free (mime_types);
		if (tmp != NULL)
			return lang;
		
		languages = g_slist_next (languages);
	}

	return NULL;
}

static GSList *
get_lang_files (GtkSourceLanguagesManager *lm)
{
	GSList *filenames = NULL;
	GSList *dirs;

	g_return_val_if_fail (lm->priv->language_specs_directories != NULL, NULL);

	dirs = lm->priv->language_specs_directories;

	while (dirs != NULL)
	{
		filenames = build_file_listing ((const gchar*)dirs->data,
						filenames);

		dirs = g_slist_next (dirs);
	}

	return filenames;
}

static GSList *
build_file_listing (const gchar *directory, GSList *filenames)
{
	GDir *dir;
	const gchar *file_name;
	
	dir = g_dir_open (directory, 0, NULL);
	
	if (dir == NULL)
		return filenames;

	file_name = g_dir_read_name (dir);
	
	while (file_name != NULL)
	{
		gchar *full_path = g_build_filename (directory, file_name, NULL);
		gchar *last_dot = strrchr (full_path, '.');
		
		if (!g_file_test (full_path, G_FILE_TEST_IS_DIR) && 
		    last_dot && 
		    (strcmp (last_dot + 1, "lang") == 0))
			filenames = g_slist_prepend (filenames, full_path);
		else
			g_free (full_path);

		file_name = g_dir_read_name (dir);
	}

	g_dir_close (dir);

	return filenames;
}

static void
slist_deep_free (GSList *list)
{
	g_slist_foreach (list, (GFunc) g_free, NULL);
	g_slist_free (list);
}

