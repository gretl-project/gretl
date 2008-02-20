
/* Generated data (by glib-mkenums) */

#include <glib-object.h>
#include "gtksourceview-typebuiltins.h"


/* enumerations from "gtksourceiter.h" */
static const GFlagsValue _gtk_source_search_flags_values[] = {
  { GTK_SOURCE_SEARCH_VISIBLE_ONLY, "GTK_SOURCE_SEARCH_VISIBLE_ONLY", "visible-only" },
  { GTK_SOURCE_SEARCH_TEXT_ONLY, "GTK_SOURCE_SEARCH_TEXT_ONLY", "text-only" },
  { GTK_SOURCE_SEARCH_CASE_INSENSITIVE, "GTK_SOURCE_SEARCH_CASE_INSENSITIVE", "case-insensitive" },
  { 0, NULL, NULL }
};

GType
gtk_source_search_flags_get_type (void)
{
  static GType type = 0;

  if (!type)
    type = g_flags_register_static ("GtkSourceSearchFlags", _gtk_source_search_flags_values);

  return type;
}


/* enumerations from "gtksourcetagstyle.h" */
static const GFlagsValue _gtk_source_tag_style_mask_values[] = {
  { GTK_SOURCE_TAG_STYLE_USE_BACKGROUND, "GTK_SOURCE_TAG_STYLE_USE_BACKGROUND", "use_background" },
  { GTK_SOURCE_TAG_STYLE_USE_FOREGROUND, "GTK_SOURCE_TAG_STYLE_USE_FOREGROUND", "use_foreground" },
  { 0, NULL, NULL }
};

GType
gtk_source_tag_style_mask_get_type (void)
{
  static GType type = 0;

  if (!type)
    type = g_flags_register_static ("GtkSourceTagStyleMask", _gtk_source_tag_style_mask_values);

  return type;
}


/* Generated data ends here */

