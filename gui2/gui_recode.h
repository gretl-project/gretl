#ifndef GUI_RECODE_H
#define GUI_RECODE_H

gchar *my_filename_from_utf8 (char *fname);

gchar *my_locale_from_utf8 (const gchar *src);

gchar *my_filename_to_utf8 (const char *fname);

gchar *my_locale_to_utf8 (const gchar *src);

gchar *my_locale_to_utf8_next (const gchar *src);

int maybe_recode_gp_file_to_utf8 (const char *fname);

gchar *gp_locale_from_utf8 (const gchar *src);

#endif /* GUI_RECODE_H */
