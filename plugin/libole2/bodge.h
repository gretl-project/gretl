/* POSIX/Non-ANSI names for increased portability */
#define O_RDONLY        _O_RDONLY
#define O_WRONLY        _O_WRONLY
#define O_RDWR          _O_RDWR
#define O_ACCMODE       _O_ACCMODE
#define O_APPEND        _O_APPEND
#define O_CREAT         _O_CREAT
#define O_TRUNC         _O_TRUNC
#define O_EXCL          _O_EXCL
#define O_TEXT          _O_TEXT
#define O_BINARY        _O_BINARY
#define O_TEMPORARY     _O_TEMPORARY
#define O_NOINHERIT     _O_NOINHERIT
#define O_SEQENTIAL     _O_SEQUENTIAL
#define O_RANDOM        _O_RANDOM

int             access (const char*, int);
int             chsize (int, long );
int             close (int);
int             creat (const char*, int);
int             dup (int);
int             dup2 (int, int);
int             eof (int);
long            filelength (int);
int             fileno (FILE*);
int             isatty (int);
long            lseek (int, long, int);
int             open (const char*, int, ...);
int             read (int, void*, unsigned int);
int             sopen (const char*, int, int, ...);
long            tell (int);
int             umask (int);
int             unlink (const char*);
int             write (int, const void*, unsigned int);

