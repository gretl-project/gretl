#ifndef UPDATER_H_
#define UPDATER_H_

#define GRETL_BUFSIZE 8192

#define _(String) String

enum {
    SP_NONE, 
    SP_LOAD_INIT,
    SP_SAVE_INIT,
    SP_FONT_INIT,
    SP_FINISH 
} progress_flags;

enum {
    SP_RETURN_OK,
    SP_RETURN_DONE,
    SP_RETURN_CANCELED
} progress_return_flags;

extern FILE *flg;
extern int logit;

int errbox (const char *msg);
int infobox (const char *msg);
int untgz (char *fname);
int show_progress (long res, long expected, int flag);

#endif /* UPDATER_H_ */
