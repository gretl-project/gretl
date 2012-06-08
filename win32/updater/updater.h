#ifndef UPDATER_H_
#define UPDATER_H_

#define _(String) String
#define I_(String) String

enum progress_flags {
    SP_NONE,
    SP_TOTAL,
    SP_LOAD_INIT,
    SP_SAVE_INIT,
    SP_FONT_INIT,
    SP_UPDATER_INIT,
    SP_FINISH 
};

enum {
    SP_RETURN_OK,
    SP_RETURN_DONE,
    SP_RETURN_CANCELED
} progress_return_flags;

extern FILE *flg;
extern int logit;

char *gretl_strdup (const char *src);
FILE *gretl_fopen (const char *filename, const char *mode);
int gretl_remove (const char *filename);
int errbox (const char *msg);
int infobox (const char *msg);
int gretl_untar (const char *fname);
int show_progress (long res, long expected, int flag);

#endif /* UPDATER_H_ */
