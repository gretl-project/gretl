#ifndef UPDATER_H_
#define UPDATER_H_

int retrieve_url (int opt, char *fname, char **savebuf, char *localfile,
		  char *errbuf, time_t filedate);
void clear (char *str, const int len);
int untgz (char *fname);
void *mymalloc (size_t size);
#ifdef OS_WIN32
int read_reg_val (HKEY tree, char *keyname, char *keyval);
#endif

#endif /* UPDATER_H_ */
