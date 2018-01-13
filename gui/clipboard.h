#ifndef CLIPBOARD_H
#define CLIPBOARD_H

int prn_to_clipboard (PRN *prn, int fmt);

int buf_to_clipboard (const char *buf);

void flag_image_available (void);

int image_file_to_clipboard (const char *fname);

#endif /* CLIPBOARD_H */
