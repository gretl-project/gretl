#ifndef MENUSTATE_H
#define MENUSTATE_H

void refresh_data (void);

void edit_info_state (gboolean s);
void add_remove_markers_state (gboolean s);
void main_menubar_state (gboolean s);
void time_series_menu_state (gboolean s);
void panel_menu_state (gboolean s);
void ts_or_panel_menu_state (gboolean s);
void session_menu_state (gboolean s);
void restore_sample_state (gboolean s);
void compact_data_state (gboolean s);

GtkWidget *build_var_popup (void);
GtkWidget *build_selection_popup (void);

void clear_sample_label (void);
void set_sample_label (DATAINFO *pdinfo);

#endif /* MENUSTATE_H */
