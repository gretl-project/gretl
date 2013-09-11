/* private stuff for gretl CMD structure */

#ifndef CMD_PRIVATE_H
#define CMD_PRIVATE_H

#include "gretl_restrict.h"
#include "gretl_func.h"

typedef struct Laginfo_ Laginfo;

typedef enum {
    CMD_NOLIST  = 1 << 0, /* command doesn't have a list of variables */
    CMD_IGNORE  = 1 << 1, /* line should be ignored */
    CMD_NULLIST = 1 << 2, /* command has been given a null list on input */
    CMD_SUBST   = 1 << 3, /* string substitution has been done on command */
    CMD_PROG    = 1 << 4, /* command is in context of progressive loop */
    CMD_CATCH   = 1 << 5, /* error from command should be "caught" */
    CMD_NOSUB   = 1 << 6, /* no @-substitution wanted (pre-checked) */
    CMD_NOOPT   = 1 << 7  /* no options present (pre-checked) */
} CmdFlags;

#define cmd_nolist(c)  (c->flags & CMD_NOLIST)
#define cmd_ignore(c)  (c->flags & CMD_IGNORE)
#define cmd_subst(c)   (c->flags & CMD_SUBST)
#define cmd_nosub(c)   (c->flags & CMD_NOSUB)
#define cmd_noopt(c)   (c->flags & CMD_NOOPT)

struct CMD_ {
    char word[FN_NAMELEN];      /* command word */
    int ci;                     /* command index number */
    int err;                    /* error code */
    int context;                /* context for subsetted commands */
    gretlopt opt;               /* option flags */
    int order;                  /* lag order, for various commands */
    int aux;                    /* auxiliary int (e.g. for VECM rank) */
    CmdFlags flags;             /* internal flags */
    char *param;                /* general-purpose parameter to command */
    char *parm2;                /* second parameter for some special uses */
    int *list;                  /* list of variables by ID number */
    int *auxlist;               /* auxiliary list for some uses */
    char savename[MAXSAVENAME]; /* for object-saving mechanism */
    Laginfo *linfo;             /* struct for recording info on automatically
                                   generated lags */
};

typedef void (*EXEC_CALLBACK) (ExecState *, void *, GretlObjType type);

struct ExecState_ {
    ExecFlags flags;
    CMD *cmd;
    PRN *prn;
    char *line;
    char runfile[MAXLEN];
    MODEL *model;          /* "workspace" model */
    MODEL *pmod;           /* set if new model is estimated */
    equation_system *sys;
    gretl_restriction *rset;
    GRETL_VAR *var;
    void *prev_model;
    GretlObjType prev_type;
    int prev_model_count;
    char *submask;        /* record of incoming sub-sample for functions */
    int padded;           /* record of incoming panel padding, if any */
    int in_comment;
    EXEC_CALLBACK callback;
};

void gretl_exec_state_init (ExecState *s,
			    ExecFlags flags,
			    char *line,
			    CMD *cmd,
			    MODEL *model, 
			    PRN *prn);

void function_state_init (CMD *cmd, ExecState *state, int *indent0);

void gretl_exec_state_set_callback (ExecState *s, EXEC_CALLBACK callback,
				    gretlopt opt);

EXEC_CALLBACK get_gui_callback (void);

void gretl_exec_state_clear (ExecState *s);

void gretl_exec_state_uncomment (ExecState *s);

void gretl_exec_state_transcribe_flags (ExecState *s, CMD *cmd);

void gretl_exec_state_set_model (ExecState *s, MODEL *pmod);

int process_command_error (CMD *cmd, int err);

int maybe_exec_line (ExecState *s, DATASET *dset);

int plausible_genr_start (const char *s, const DATASET *dset);

#endif /* CMD_PRIVATE_H */
