#ifndef CLIPBOARD_H
#define CLIPBOARD_H

int prn_to_clipboard (PRN *prn, int fmt);

int buf_to_clipboard (const char *buf);

#endif /* CLIPBOARD_H */
