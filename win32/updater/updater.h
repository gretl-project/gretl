#ifndef UPDATER_H_
#define UPDATER_H_

#define GRETL_BUFSIZE 8192

extern FILE *flg;
extern int logit;

int errbox (const char *msg);
int infobox (const char *msg);
int untgz (char *fname);
void update_windows_progress_bar (int gotbytes);

#endif /* UPDATER_H_ */
