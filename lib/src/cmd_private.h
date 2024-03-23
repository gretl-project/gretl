/* private stuff for gretl CMD structure */

#ifndef CMD_PRIVATE_H
#define CMD_PRIVATE_H

#include "gretl_restrict.h"
#include "gretl_func.h"
#include "monte_carlo.h"

/* As of 2024-03-23, compilation seems to be working correctly for the
   special functions feval() and fevalb(), which had previously been
   leaking memory when compiled (though not otherwise). For now we're
   leaving this define here in case of trouble. The symbol
   COMPILE_FEVAL is referenced in cmd_private.c and monte_carlo.c:
   setting it to 0 will prevent compilation in both contexts.
*/
#define COMPILE_FEVAL 1

typedef enum {
    CMD_CCMT    = 1 << 0, /* line is in a C-style multi-line comment */
    CMD_SUBST   = 1 << 1, /* string substitution has been done on command */
    CMD_PROG    = 1 << 2, /* command is in context of progressive loop */
    CMD_CATCH   = 1 << 3, /* error from command should be "caught" */
    CMD_NOSUB   = 1 << 4, /* no @-substitution called for (pre-checked) */
    CMD_ENDFUN  = 1 << 5  /* line terminates a function definition */
} CmdFlags;

#define cmd_ccmt(c)  (c->flags & CMD_CCMT)
#define cmd_subst(c) (c->flags & CMD_SUBST)
#define cmd_nosub(c) (c->flags & CMD_NOSUB)

typedef struct cmd_token_ cmd_token;

struct CMD_ {
    int ci;          /* current command index */
    int err;         /* error code */
    int context;     /* for block commands, index of current context */
    int ciflags;     /* see CIFlags in tokenizer.c */
    gretlopt opt;    /* option(s) for command */
    CmdFlags flags;  /* status flags for command invocation */
    GretlType gtype; /* specified type for "genr" */
    int order;       /* lag order, where appropriate */
    int auxint;      /* auxiliary int (e.g. VECM rank) */
    int cstart;      /* token index of start of 'real' command */
    int ntoks;       /* number of tokens actually used */
    int nt_alloced;  /* number of tokens allocated */
    cmd_token *toks;    /* tokens */
    const char *vstart; /* pointer to where in line varargs or expr start */
    char *param;     /* basic parameter string */
    char *parm2;     /* second parameter string */
    int *list;       /* list of series and/or control integers */
    int *auxlist;    /* needed for "gappy" lag lists */
    char savename[MAXSAVENAME]; /* for object-saving mechanism */
};

typedef int (*EXEC_CALLBACK) (ExecState *, void *, GretlObjType type);

struct ExecState_ {
    ExecFlags flags;
    CMD *cmd;
    PRN *prn;
    char *line;
    char runfile[MAXLEN];
    MODEL *model;          /* "workspace" model */
    MODEL *pmod;           /* set if new model is estimated */
    LOOPSET *loop;
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

struct OpenOp_ {
    char fname[MAXLEN];
    int quiet;
    int http;
    int dbdata;
    int ftype;
};

typedef struct OpenOp_ OpenOp;

/* structure representing a line of a function or loop */
struct stmt_ {
    char *s;       /* text of command line */
    gint16 ci;     /* command index */
    gint16 next;   /* line index to skip to after loop or conditional */
    gint16 idx;    /* 1-based line index allowing for blanks (functions only) */
    gint16 flags;  /* state flags (different for functions vs loops) */
    void *ptr;     /* attachment for LOOPSET or GENERATOR */
};

typedef struct stmt_ stmt;

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

void set_gui_callback (EXEC_CALLBACK callback);

void gretl_exec_state_clear (ExecState *s);

void gretl_exec_state_destroy (ExecState *s);

void gretl_exec_state_uncomment (ExecState *s);

void gretl_exec_state_transcribe_flags (ExecState *s, CMD *cmd);

void gretl_exec_state_set_model (ExecState *s, MODEL *pmod);

int process_command_error (ExecState *s, int err);

int maybe_exec_line (ExecState *s, DATASET *dset, void *ptr);

int plausible_genr_start (const char *s, const DATASET *dset);

int statements_get_structure (stmt *lines,
			      int n_lines,
			      int context,
			      const char *name);

int cmd_arg1_quoted (CMD *cmd);

#endif /* CMD_PRIVATE_H */
