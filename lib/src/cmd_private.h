/* private stuff for gretl CMD structure */

#ifndef CMD_PRIVATE_H
#define CMD_PRIVATE_H

#include "gretl_restrict.h"
#include "gretl_func.h"

typedef struct Laginfo_ Laginfo;

enum {
    CMD_NOLIST  = 1 << 0, /* command doesn't have a list of variables */
    CMD_IGNORE  = 1 << 1, /* line should be ignored */
    CMD_NULLIST = 1 << 2, /* command has been given a null list on input */
    CMD_SUBST   = 1 << 3, /* string substitution has been done on command */
    CMD_PROG    = 1 << 4  /* command is in context of progressive loop */
};

#define cmd_nolist(c)    (c->flags & CMD_NOLIST)
#define cmd_ignore(c)    (c->flags & CMD_IGNORE)
#define cmd_subst(c)     (c->flags & CMD_SUBST)

struct CMD_ {
    char word[FN_NAMELEN];      /* command word */
    int ci;                     /* command index number */
    int context;                /* context for subsetted commands */
    int order;                  /* lag order, for various commands */
    int aux;                    /* auxiliary int (e.g. for VECM rank) */
    gretlopt opt;               /* option flags */
    char flags;                 /* internal flags */
    char savename[MAXSAVENAME]; /* name used to save an object from the command */
    int *list;                  /* list of variables by ID number */
    char *param;                /* general-purpose parameter to command */
    char *extra;                /* second parameter for some special uses */
    int err;                    /* error code */
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
    MODEL **models;        /* "workspace" models */
    MODEL *pmod;           /* set if new model is estimated */
    equation_system *sys;
    gretl_restriction *rset;
    GRETL_VAR *var;
    void *prev_model;
    GretlObjType prev_type;
    char *submask;        /* record of incoming sub-sample for functions */
    int in_comment;
    int funcerr;
    EXEC_CALLBACK callback;
};

void gretl_exec_state_init (ExecState *s,
			    ExecFlags flags,
			    char *line,
			    CMD *cmd,
			    MODEL **models, 
			    PRN *prn);

void gretl_exec_state_set_callback (ExecState *s, EXEC_CALLBACK callback);

void gretl_exec_state_clear (ExecState *s);

void gretl_exec_state_uncomment (ExecState *s);

int maybe_exec_line (ExecState *s, double ***pZ, DATAINFO *pdinfo);

int plausible_genr_start (const char *s, const DATAINFO *pdinfo);

#endif /* CMD_PRIVATE_H */
