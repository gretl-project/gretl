#ifndef CLIPBOARD_H
#define CLIPBOARD_H

extern gchar *clipboard_buf; 

void gretl_clipboard_free (void);

void gretl_clipboard_set (int copycode);

#endif /* CLIPBOARD_H */
