/* private stuff for gretl CMD structure */

#ifndef CMD_PRIVATE_H
#define CMD_PRIVATE_H

#undef USE_LAGINFO

#if USE_LAGINFO
typedef struct Laginfo_ Laginfo;
#endif

struct _CMD {
    char word[9];               /* command word */
    int ci;                     /* command index number */
    int context;                /* context for subsetted commands */
    gretlopt opt;               /* option flags */
    char savename[MAXSAVENAME]; /* name used to save an object from the command */
    char str[4];                /* used, e.g., in "multiply" command */
    int nolist;                 /* = 1 if the command does not take a list */
    int *list;                  /* list of variables by ID number */
    char *param;                /* general-purpose parameter to command */
    int errcode;                /* error code */
#if USE_LAGINFO
    Laginfo *linfo;
#endif
};

#endif /* CMD_PRIVATE_H */
