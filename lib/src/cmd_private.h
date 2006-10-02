/* private stuff for gretl CMD structure */

#ifndef CMD_PRIVATE_H
#define CMD_PRIVATE_H

#include "gretl_restrict.h"

typedef struct Laginfo_ Laginfo;

struct CMD_ {
    char word[9];               /* command word */
    int ci;                     /* command index number */
    int context;                /* context for subsetted commands */
    int order;                  /* lag order, for various commands */
    gretlopt opt;               /* option flags */
    char savename[MAXSAVENAME]; /* name used to save an object from the command */
    int nolist;                 /* = 1 if the command does not take a list */
    int *list;                  /* list of variables by ID number */
    char *param;                /* general-purpose parameter to command */
    char *extra;                /* second parameter for some special uses */
    int ignore;                 /* flag set to ignore comment lines in input */
    int errcode;                /* error code */
    Laginfo *linfo;             /* struct for recording info on automatically
                                   generated lags */
};

typedef void (*EXEC_CALLBACK) (ExecState *, double ***, DATAINFO *);

struct ExecState_ {
    ExecFlags flags;
    CMD *cmd;
    PRN *prn;
    char *line;
    char runfile[MAXLEN];
    MODEL **models;
    gretl_equation_system *sys;
    gretl_restriction_set *rset;
    GRETL_VAR *var;
    int alt_model;
    int in_comment;
    EXEC_CALLBACK callback;
};

void gretl_exec_state_init (ExecState *s,
			    ExecFlags flags,
			    char *line,
			    CMD *cmd,
			    MODEL **models, 
			    PRN *prn);

void gretl_exec_state_clear (ExecState *s);

#endif /* CMD_PRIVATE_H */
