#ifndef UPDATER_H_
#define UPDATER_H_

#define GRETL_BUFSIZE 8192

extern FILE *flg;
extern int logit;

int untgz (char *fname);

#endif /* UPDATER_H_ */
