#ifndef TOOLBAR_H
#define TOOLBAR_H

void show_toolbar (void);

void gretl_pdf_manual (void);

void vwin_add_viewbar (windata_t *vwin, int text_out);

void viewbar_add_edit_items (windata_t *vwin);

GtkWidget *build_text_popup (windata_t *vwin);

void gretl_stock_icons_init (void);

#endif /* TOOLBAR_H */
