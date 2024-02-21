/* Code for the "new tokenizer", which replaced the old command
   parser in interact.c in September 2014, following the 1.9.92
   release of gretl.
*/

#define CDEBUG 0
#define TDEBUG 0

typedef enum {
    CI_LIST  = 1 << 0,  /* list may be present */
    CI_LLEN1 = 1 << 1,  /* list must contain exactly 1 member */
    CI_LLEN2 = 1 << 2,  /* list must contain exactly 2 members */
    CI_ORD1  = 1 << 3,  /* args start with order specifier */
    CI_ORD2  = 1 << 4,  /* args (may) end with order */
    CI_L1INT = 1 << 5,  /* first portion of list contains ints */
    CI_PARM1 = 1 << 6,  /* args start with (non-list) parameter */
    CI_PARM2 = 1 << 7,  /* may use second param */
    CI_FNAME = 1 << 8,  /* (first) parameter is filename */
    CI_EXPR  = 1 << 9,  /* uses "genr"-type expression */
    CI_VARGS = 1 << 10, /* ends with varargs field */
    CI_EXTRA = 1 << 11, /* uses some special "extra" feature */
    CI_ADHOC = 1 << 12, /* needs special-purpose parser */
    CI_DOALL = 1 << 13, /* operates on all series if no list given */
    CI_NOOPT = 1 << 14, /* command never takes any options */
    CI_BLOCK = 1 << 15, /* command starts a block */
    CI_FFORM = 1 << 16, /* command also has function-form */
    CI_LCHK  = 1 << 17, /* needs checking for "list" specials */
    CI_INFL  = 1 << 18, /* command arglist "inflected" by options */
    CI_FCMIN = 1 << 19, /* minimal (single word) flow control */
    CI_LGEN  = 1 << 20, /* command generates a named list */
    CI_OBSOL = 1 << 21  /* command is deprecated */
} CIFlags;

struct gretl_cmd {
    int cnum;
    const char *cword;
    CIFlags flags;
};

/* Note: the flags CI_EXPR, CI_VARGS and CI_ADHOC all indicate that
   the command in question carries (or may carry) material that gets
   passed on for parsing elsewhere (i.e. not in this translation
   unit).

   In the case of CI_EXPR and CI_VARGS, this material (i) is a
   strictly trailing portion of the command line and (ii) is destined
   for handling via "genr". Such commands never take options (so the
   CI_NOOPT flag would be redundant). The unparsed material can
   therefore simply take the form of a const char pointer into the
   command line (namely, the @vstart member of the CMD struct). The
   CI_VARGS case differs from CI_EXPR in that the former requires the
   parsing here of either one or two leading parameters; there are no
   such parameters with CI_EXPR.

   CI_ADHOC is used for cases where a certain chunk of the command
   line needs to be packaged up for specialized parsing
   elsewhere. Options may follow the chunk in question; it's not a
   strictly trailing portion of the line. The ad hoc material is
   packaged into either the @param or @parm2 member of the CMD struct.

   CI_EXTRA is used where some (possibly optional) element of a
   command does not fit neatly into the slots of @list, @param and
   @parm2. The distinction between CI_ADHOC and CI_EXTRA is not
   totally rigorous, but it seems to be useful to have both flags;
   the idea is that CI_EXTRA is used when the "special" portion of
   a command line is a limited subset of the line.
*/

static struct gretl_cmd gretl_cmds[] = {
    { 0,        "",         0},
    { ADD,      "add",      CI_LIST },
    { ADF,      "adf",      CI_ORD1 | CI_LIST },
    { ANOVA,    "anova",    CI_LIST },
    { APPEND,   "append",   CI_PARM1 | CI_FNAME },
    { AR,       "ar",       CI_LIST | CI_L1INT },
    { AR1,      "ar1",      CI_LIST },
    { ARCH,     "arch",     CI_ORD1 | CI_LIST },
    { ARMA,     "arima",    CI_LIST | CI_L1INT },
    { BDS,      "bds",      CI_ORD1 | CI_LIST | CI_LLEN1 },
    { BIPROBIT, "biprobit", CI_LIST },
    { BKW,      "bkw",      0 },
    { BREAK,    "break",    CI_NOOPT | CI_FCMIN },
    { BXPLOT,   "boxplot",  CI_LIST | CI_EXTRA | CI_INFL },
    { CHOW,     "chow",     CI_PARM1 },
    { CLEAR,    "clear",    0 },
    { COEFFSUM, "coeffsum", CI_LIST  },
    { COINT,    "coint",    CI_ORD1 | CI_LIST },
    { COINT2,   "johansen", CI_ORD1 | CI_LIST },
    { CONTINUE, "continue", CI_NOOPT | CI_FCMIN },
    { CORR,     "corr",     CI_LIST | CI_DOALL },
    { CORRGM,   "corrgm",   CI_LIST | CI_LLEN1 | CI_ORD2 },
    { CUSUM,    "cusum",    0 },
    { DATA,     "data",     CI_ADHOC }, /* special: needs whole line */
    { DATAMOD,  "dataset",  CI_PARM1 | CI_LIST | CI_PARM2 },
    { DELEET,   "delete",   CI_PARM1 | CI_INFL }, /* or CI_LIST */
    { DIFF,     "diff",     CI_LIST },
    { DIFFTEST, "difftest", CI_LIST | CI_LLEN2 },
    { DISCRETE, "discrete", CI_LIST },
    { DPANEL  , "dpanel",   CI_LIST | CI_L1INT },
    { DUMMIFY,  "dummify",  CI_LIST },
    { DURATION, "duration", CI_LIST },
    { ELIF,     "elif",     CI_EXPR },
    { ELSE,     "else",     CI_NOOPT | CI_FCMIN },
    { END,      "end",      CI_PARM1 },
    { ENDIF,    "endif",    CI_NOOPT | CI_FCMIN },
    { ENDLOOP,  "endloop",  CI_NOOPT | CI_FCMIN },
    { EQNPRINT, "eqnprint", 0 }, /* special, handled later */
    { EQUATION, "equation", CI_LIST },
    { ESTIMATE, "estimate", CI_PARM1 | CI_PARM2 }, /* params optional */
    { EVAL,     "eval",     CI_EXPR },
    { FCAST,    "fcast",    CI_ADHOC },
    { FLUSH,    "flush",    CI_NOOPT },
    { FOREIGN,  "foreign",  CI_PARM1 | CI_BLOCK },
    { FRACTINT, "fractint", CI_LIST | CI_LLEN1 | CI_ORD2 },
    { FREQ,     "freq",     CI_LIST | CI_LLEN1 },
    { FUNC,     "function", CI_ADHOC | CI_NOOPT | CI_BLOCK },
    { FUNCERR,  "funcerr",  CI_PARM1 | CI_NOOPT },
    { GARCH,    "garch",    CI_LIST | CI_L1INT },
    { GENR,     "genr",     CI_EXPR },
    { GMM,      "gmm",      CI_EXPR | CI_BLOCK },
    { GNUPLOT,  "gnuplot",  CI_LIST | CI_EXTRA | CI_INFL },
    { GPBUILD,  "gpbuild",  CI_PARM1 | CI_BLOCK },
    { GRAPHPG,  "graphpg",  CI_PARM1 | CI_PARM2 }, /* params optional */
    { GRIDPLOT, "gridplot", CI_PARM1 },
    { HECKIT,   "heckit",   CI_LIST },
    { HELP,     "help",     CI_PARM1 },
    { HFPLOT,   "hfplot",   CI_LIST | CI_EXTRA },
    { HSK,      "hsk",      CI_LIST },
    { HURST,    "hurst",    CI_LIST | CI_LLEN1 },
    { IF,       "if",       CI_EXPR },
    { INCLUDE,  "include",  CI_PARM1 | CI_FNAME },
    { INFO,     "info",     0 },
    { INTREG,   "intreg",   CI_LIST },
    { JOIN,     "join",     CI_PARM1 | CI_FNAME | CI_EXTRA },
    { KDPLOT,   "kdplot",   CI_LIST | CI_LLEN1 },
    { KPSS,     "kpss",     CI_ORD1 | CI_LIST },
    { LABELS,   "labels",   CI_LIST | CI_DOALL },
    { LAD,      "lad",      CI_LIST },
    { LAGS,     "lags",     CI_ORD1 | CI_LIST },
    { LDIFF,    "ldiff",    CI_LIST | CI_NOOPT },
    { LEVERAGE, "leverage", 0 },
    { LEVINLIN, "levinlin", CI_PARM1 | CI_LIST | CI_LLEN1 },
    { LOGISTIC, "logistic", CI_LIST },
    { LOGIT,    "logit",    CI_LIST },
    { LOGS,     "logs",     CI_LIST | CI_NOOPT },
    { LOOP,     "loop",     CI_ADHOC }, /* ? */
    { MAHAL,    "mahal",    CI_LIST },
    { MAKEPKG,  "makepkg",  CI_PARM1 },
    { MARKERS,  "markers",  0 },
    { MEANTEST, "meantest", CI_LIST | CI_LLEN2 },
    { MIDASREG, "midasreg", CI_LIST },
    { MLE,      "mle",      CI_EXPR | CI_BLOCK },
    { MODELTAB, "modeltab", CI_PARM1 | CI_INFL },
    { MODPRINT, "modprint", CI_PARM1 | CI_PARM2 | CI_EXTRA },
    { MODTEST,  "modtest",  CI_ORD1 },
    { MPI,      "mpi",      CI_BLOCK },
    { MPOLS,    "mpols",    CI_LIST },
    { NEGBIN,   "negbin",   CI_LIST },
    { NLS,      "nls",      CI_EXPR | CI_BLOCK },
    { NORMTEST, "normtest", CI_LIST | CI_LLEN1 },
    { NULLDATA, "nulldata", CI_ORD1 },
    { OLS,      "ols",      CI_LIST },
    { OMIT,     "omit",     CI_LIST },
    { OPEN,     "open",     CI_PARM1 | CI_FNAME | CI_INFL }, /* + ODBC specials */
    { ORTHDEV,  "orthdev",  CI_LIST | CI_NOOPT },
    { OUTFILE,  "outfile",  CI_PARM1 | CI_FNAME | CI_INFL | CI_BLOCK },
    { PANEL,    "panel",    CI_LIST },
    { PANPLOT,  "panplot",  CI_LIST | CI_LLEN1 | CI_EXTRA },
    { PANSPEC,  "panspec",  0 },
    { PCA,      "pca",      CI_LIST | CI_DOALL },
    { PERGM,    "pergm",    CI_LIST | CI_LLEN1 | CI_ORD2 },
    { PLOT,     "plot",     CI_BLOCK | CI_PARM1 },
    { POISSON,  "poisson",  CI_LIST },
    { PRINT,    "print",    CI_INFL }, /* special: handled later */
    { PRINTF,   "printf",   CI_PARM1 | CI_VARGS },
    { PROBIT,   "probit",   CI_LIST },
    { PVAL,     "pvalue",   CI_ADHOC | CI_NOOPT },
    { QUANTREG, "quantreg", CI_PARM1 | CI_LIST },
    { QLRTEST,  "qlrtest",  0 },
    { QQPLOT,   "qqplot",   CI_LIST },
    { QUIT,     "quit",     CI_NOOPT },
    { RENAME,   "rename",   CI_PARM1 | CI_PARM2 },
    { RESET,    "reset",    0 },
    { RESTRICT, "restrict", CI_PARM1 | CI_BLOCK },
    { RMPLOT,   "rmplot",   CI_LIST | CI_LLEN1 },
    { RUN,      "run",      CI_PARM1 | CI_FNAME },
    { RUNS,     "runs",     CI_LIST | CI_LLEN1 },
    { SCATTERS, "scatters", CI_LIST | CI_DOALL },
    { SDIFF,    "sdiff",    CI_LIST | CI_NOOPT },
    { SET,      "set",      CI_PARM1 | CI_PARM2 | CI_INFL },
    { SETINFO,  "setinfo",  CI_LIST | CI_LLEN1 | CI_INFL }, /* + special: handled later */
    { SETOBS,   "setobs",   CI_PARM1 | CI_PARM2 },
    { SETOPT,   "setopt",   CI_PARM1 | CI_PARM2 },
    { SETMISS,  "setmiss",  CI_PARM1 | CI_LIST | CI_DOALL },
    { SHELL,    "shell",    CI_EXPR },
    { SMPL,     "smpl",     CI_PARM1 | CI_PARM2 | CI_INFL }, /* alternate forms */
    { SPEARMAN, "spearman", CI_LIST | CI_LLEN2 },
    { SQUARE,   "square",   CI_LIST },
    { STDIZE,   "stdize",   CI_LIST },
    { STORE,    "store",    CI_PARM1 | CI_FNAME | CI_LIST | CI_DOALL },
    { SUMMARY,  "summary",  CI_LIST | CI_DOALL },
    { SYSTEM,   "system",   CI_PARM1 | CI_BLOCK },
    { TABPRINT, "tabprint", 0 }, /* special, handled later */
    { TEXTPLOT, "textplot", CI_LIST },
    { TOBIT,    "tobit",    CI_LIST },
    { IVREG,    "tsls",     CI_LIST },
    { VAR,      "var",      CI_ORD1 | CI_LIST },
    { VARLIST,  "varlist",  0 },
    { VARTEST,  "vartest",  CI_LIST | CI_LLEN2 | CI_NOOPT },
    { VECM,     "vecm",     CI_ORD1 | CI_LIST },
    { VIF,      "vif",      0 },
    { WLS,      "wls",      CI_LIST },
    { XCORRGM,  "xcorrgm",  CI_LIST | CI_LLEN2 | CI_ORD2 },
    { XTAB,     "xtab",     CI_LIST | CI_INFL },
    { FUNCRET,  "return",   CI_EXPR },
    { CATCH,    "catch",    0 },
    { PKG,      "pkg",      CI_PARM1 | CI_PARM2 },
    { NC,       NULL,       0 }
};

#define not_catchable(c) (c == IF || c == ENDIF || c == ELIF || \
			  c == FUNC)

#define param_optional(c) (c == SET || c == HELP || c == RESTRICT || \
			   c == SMPL || c == SYSTEM || c == FUNCERR || \
			   c == GRAPHPG || c == PLOT || c == OUTFILE)

#define parm2_optional(c) (c == SET || c == SETOPT || c == SETOBS || \
			   c == ESTIMATE || c == HELP || c == GRAPHPG || \
			   c == EQUATION || c == MODPRINT)

#define vargs_optional(c) (c == PRINTF)

#define expr_keep_cmdword(c) (c == GMM || c == MLE || c == NLS)

#define has_function_form(c) (gretl_cmds[c].flags & CI_FFORM)

#define option_inflected(c) (c->ciflags & CI_INFL)

#define simple_flow_control(c) (c->ciflags & CI_FCMIN)

static int command_get_flags (int ci)
{
    if (ci >= 0 && ci < NC) {
	return gretl_cmds[ci].flags;
    } else {
	return 0;
    }
}

static int never_takes_options (CMD *c)
{
    return c->ciflags & (CI_EXPR | CI_VARGS | CI_NOOPT);
}

static void check_for_shadowed_commands (void)
{
    int i;

    for (i=1; i<NC; i++) {
	if (function_lookup(gretl_cmds[i].cword)) {
	    gretl_cmds[i].flags |= CI_FFORM;
	}
    }
}

/* Get the maximum number of (semicolon) separators
   supported by the command with index @ci. For most
   commands this is zero. If @pmin is non-NULL it
   gets the minimum number of such separators needed
   by the command.
*/

static int get_sep_max (int ci, int *pmin)
{
    /* default: semicolon separator neither required nor allowed */
    int minsep = 0, maxsep = 0;

    switch (ci) {
    case AR:
    case GARCH:
    case HECKIT:
    case IVREG:
    case MIDASREG:
	minsep = maxsep = 1;
	break;
    case DPANEL:
    case ARMA:
	minsep = 1;
	maxsep = 2;
	break;
    case COINT2:
    case VECM:
	maxsep = 2;
	break;
    case BIPROBIT:
    case DURATION:
    case EQUATION:
    case MPOLS:
    case NEGBIN:
    case POISSON:
    case VAR:
    case XTAB:
    case LAGS:
    case SCATTERS:
    case HFPLOT:
    case SMPL:
	maxsep = 1;
	break;
    default:
	break;
    }

    if (pmin != NULL) {
	*pmin = minsep;
    }

    return maxsep;
}

/* The difference between TOK_JOINED and TOK_NOGAP below is that we
   don't count a comma with no gap to its left as JOINED, since it's
   generally punctuation, but we do record it as NOGAP, since this
   info can be relevant when reconstructing an option parameter.
*/

enum {
    TOK_JOINED = 1 << 0, /* token is joined on the left (no space) */
    TOK_NOGAP  = 1 << 1, /* as TOK_JOINED but including joined comma */
    TOK_DONE   = 1 << 2, /* token has been handled */
    TOK_QUOTED = 1 << 3, /* token was found in double quotes */
    TOK_IGNORE = 1 << 4, /* token not actually wanted, ignored */
    TOK_LSTR   = 1 << 5  /* token provisionally added to list string */
} TokenFlags;

enum {
    TOK_NAME,    /* potentially valid identifier (not quoted) */
    TOK_DOLSTR,  /* '$' plus potentially valid identifier */
    TOK_OPT,     /* long-form option flag */
    TOK_SOPT,    /* short-form option flag */
    TOK_CATCH,   /* "catch" keyword */
    TOK_STRING,  /* string in double quotes */
    TOK_NUMBER,  /* numeric, plus 'colonized' dates such as 1990:1 */
    TOK_INT,     /* integer (may start with '+' or '-') */
    TOK_DASH,    /* single dash */
    TOK_DDASH,   /* double dash */
    TOK_DPLUS,   /* double plus */
    TOK_EQUALS,  /* equals sign by itself */
    TOK_EQMOD,   /* "+=", "-=", etc. */
    TOK_ASSIGN,  /* assignment to object, "<-" */
    TOK_SEMIC,   /* semicolon */
    TOK_COLON,   /* colon */
    TOK_DOT,     /* single dot */
    TOK_COMMA,   /* single comma */
    TOK_DDOT,    /* double dot */
    TOK_PRSTR,   /* string in parentheses */
    TOK_BRSTR,   /* string in square brackets */
    TOK_CBSTR,   /* string in curly braces */
    TOK_OPTDASH, /* dash preceding an option flag */
    TOK_OPTEQ,   /* '=' that joins option flag to value */
    TOK_OPTVAL,  /* value attached to option flag */
    TOK_AST,     /* single asterisk */
    TOK_SYMB,    /* symbols, not otherwise handled */
    TOK_EVAL     /* string that needs to be eval'd */
} TokenTypes;

struct cmd_token_ {
    char *s;         /* allocated token string */
    const char *lp;  /* pointer to line position */
    guint8 type;     /* one of TokenTypes */
    guint8 flag;     /* zero or more of TokenFlags */
};

#define token_joined(t)  (t->flag & TOK_JOINED)
#define token_done(t)    (t->flag & TOK_DONE)
#define token_ignored(t) (t->flag & TOK_IGNORE)

#define mark_token_done(t) (t.flag |= TOK_DONE)
#define mark_list_token_done(t) (t.flag |= (TOK_DONE|TOK_LSTR))
#define mark_token_ignored(t) (t.flag |= (TOK_DONE|TOK_IGNORE))

#define option_type(t) (t == TOK_OPT || t == TOK_SOPT || \
			t == TOK_OPTDASH || t == TOK_OPTEQ || \
			t == TOK_OPTVAL)

#define delimited_type(t) (t == TOK_STRING || \
			   t == TOK_PRSTR ||  \
			   t == TOK_CBSTR ||  \
			   t == TOK_BRSTR)

#define bracketed_type(t) (t == TOK_PRSTR || t == TOK_BRSTR)

#define wildsym(t) (t->type == TOK_AST || \
		    (t->type == TOK_SYMB && t->s[0] == '?'))

static void cmd_token_init (cmd_token *t)
{
    t->s = NULL;
    t->lp = NULL;
    t->type = 0;
    t->flag = 0;
}

static void cmd_token_clear (cmd_token *t)
{
    if (t != NULL) {
	free(t->s);
	cmd_token_init(t);
    }
}

void gretl_cmd_free (CMD *cmd)
{
    char *s;
    int i;

    for (i=0; i<cmd->ntoks; i++) {
	s = cmd->toks[i].s;
	if (s != cmd->param && s != cmd->parm2) {
	    free(s);
	}
    }

    free(cmd->list);
    free(cmd->param);
    free(cmd->parm2);
    free(cmd->auxlist);

    free(cmd->toks);
}

void gretl_cmd_destroy (CMD *cmd)
{
    gretl_cmd_free(cmd);
    free(cmd);
}

CMD *gretl_cmd_new (void)
{
    CMD *cmd = malloc(sizeof *cmd);

    if (cmd != NULL) {
	gretl_cmd_init(cmd);
    }

    return cmd;
}

/* Tells whether the second token (typically the first
   following a command-word) of a command was quoted
   on input. This is something we may want to know when
   executing the printf command.
*/

int cmd_arg1_quoted (CMD *cmd)
{
    if (cmd->toks != NULL && cmd->ntoks > 1) {
	return (cmd->toks[1].flag & TOK_QUOTED)? 1 : 0;
    } else {
	return 0;
    }
}

static void cmd_set_vstart (CMD *c, const char *s)
{
    c->vstart = s;
}

int gretl_cmd_init (CMD *c)
{
    int i, n = 16;
    int err = 0;

    c->ci = 0;
    c->err = 0;
    c->context = 0;
    c->ciflags = 0;
    c->opt = 0;
    c->flags = 0;
    c->order = 0;
    c->auxint = 0;
    c->cstart = 0;
    c->ntoks = 0;
    c->nt_alloced = 0;
    c->toks = NULL;
    c->vstart = NULL;
    c->param = NULL;
    c->parm2 = NULL;
    c->list = NULL;
    c->auxlist = NULL;

    *c->savename = '\0';
    c->gtype = GRETL_TYPE_ANY;

    c->toks = malloc(n * sizeof *c->toks);
    if (c->toks == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<n; i++) {
	cmd_token_init(&c->toks[i]);
    }

    if (err) {
	gretl_cmd_destroy(c);
    } else {
	c->nt_alloced = n;
    }

    return err;
}

static void gretl_cmd_clear (CMD *c)
{
    cmd_token *tok;
    int i, ci = c->ci;

#if CDEBUG
    fprintf(stderr, "gretl_cmd_clear: ci = %d (%s), context = %d\n\n", ci,
	    gretl_command_word(ci), c->context);
#endif

    if (ci == END && c->param != NULL) {
	ci = gretl_command_number(c->param);
    }

    for (i=0; i<c->ntoks; i++) {
	tok = &c->toks[i];
	if (tok->s == c->param) {
	    /* avoid double-freeing */
	    c->param = NULL;
	} else if (tok->s == c->parm2) {
	    c->parm2 = NULL;
	}
	cmd_token_clear(tok);
    }

    /* FIXME: do the next step only if the CMD has actually
       been parsed/assembled, not just tokenized (as in
       get_command_index)?
    */

    if (c->context > 0) {
	; /* don't clear any pending options */
    } else if (ci > 0 && ci != SETOPT &&
	!(c->ciflags & (CI_NOOPT | CI_EXPR | CI_VARGS))) {
	clear_stored_options_for_command(ci);
    }

    c->ci = 0;
    c->err = 0;
    c->ciflags = 0;
    c->opt = 0;
    c->order = 0;
    c->auxint = 0;
    c->cstart = 0;
    c->ntoks = 0;
    c->vstart = NULL;
    c->gtype = GRETL_TYPE_ANY;

    /* Note: c->context, c->savename and the flag CMD_CATCH should
       persist across in-block commands until end-of-block is
       reached. But once we've exited a block (which is signaled by
       context == 0), these elements should be cleared.
       FIXME: same issue as above?
    */

    if (c->context == 0) {
	*c->savename = '\0';
	c->flags &= ~CMD_CATCH;
    }

    free(c->param);
    free(c->parm2);
    c->param = c->parm2 = NULL;

    free(c->list);
    free(c->auxlist);
    c->list = c->auxlist = NULL;
}

static int real_add_token (CMD *c, const char *tok,
			   const char *lp, char type,
			   char flag)
{
    int n = c->ntoks;
    int err = 0;

#if TDEBUG
    fprintf(stderr, "real_add_token: '%s' (%s)\n", tok,
	    flag & TOK_QUOTED ? "quoted" : "not quoted");
#endif

    if (n == c->nt_alloced - 1) {
	/* we've used all the existing slots */
	int i, nt_new = c->nt_alloced * 2;
	cmd_token *toks;

	toks = realloc(c->toks, nt_new * sizeof *toks);

	if (toks == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=c->nt_alloced; i<nt_new; i++) {
		cmd_token_init(&toks[i]);
	    }
	    c->toks = toks;
	    c->nt_alloced = nt_new;
	}
    }

    if (!err) {
	c->toks[n].s = gretl_strdup(tok);
	c->toks[n].lp = lp;
	c->toks[n].type = type;
	c->toks[n].flag = flag;
	c->ntoks += 1;
    }

    return err;
}

static int push_token (CMD *c, const char *tok, const char *s,
		       int pos, char type, char flag)
{
    if (pos > 0 && !isspace(*(s-1))) {
	flag |= TOK_NOGAP;
	if (type != TOK_COMMA) {
	    flag |= TOK_JOINED;
	}
    }

    return real_add_token(c, tok, s, type, flag);
}

static int push_string_token (CMD *c, const char *tok,
			      const char *s, int pos)
{
    char type = TOK_NAME;

    if (c->ntoks == 0 && !strcmp(tok, "catch")) {
	type = TOK_CATCH;
    } else if (*tok == '$') {
	type = TOK_DOLSTR;
    } else if (!strncmp(tok, "eval(", 5)) {
	char strvar[VNAMELEN];

	if (sscanf(tok, "eval(%31[^)]", strvar) == 1 &&
	    tok[strlen(tok)-1] == ')') {
	    return push_token(c, strvar, s, pos, TOK_EVAL, 0);
	}
    }

    return push_token(c, tok, s, pos, type, 0);
}

static int push_symbol_token (CMD *c, const char *tok,
			      char *s, int pos)
{
    char type = TOK_SYMB;

    if (!strcmp(tok, "-")) {
	type = TOK_DASH;
    } else if (!strcmp(tok, ".")) {
	type = TOK_DOT;
    } else if (!strcmp(tok, ",")) {
	type = TOK_COMMA;
    } else if (!strcmp(tok, ";")) {
	type = TOK_SEMIC;
    } else if (!strcmp(tok, ":")) {
	type = TOK_COLON;
    } else if (!strcmp(tok, "--")) {
	type = TOK_DDASH;
    } else if (!strcmp(tok, "++")) {
	type = TOK_DPLUS;
    } else if (!strcmp(tok, "..")) {
	type = TOK_DDOT;
    } else if (!strcmp(tok, "=")) {
	type = TOK_EQUALS;
    } else if (!strcmp(tok, "*")) {
	type = TOK_AST;
    } else if (strlen(tok) == 2 && tok[1] == '=') {
	type = TOK_EQMOD;
    } else if (c->ntoks == 1 && !strcmp(tok, "<-")) {
	/* FIXME allow for "catch" */
	type = TOK_ASSIGN;
    }

    if (type == TOK_SEMIC) {
	if (c->ci > 0 && get_sep_max(c->ci, NULL) == 0) {
	    gretl_errmsg_sprintf(_("The symbol '%c' is not valid in this context\n"),
				 ';');
	    return E_PARSE;
	}
    }

    return push_token(c, tok, s, pos, type, 0);
}

static unsigned int digit_spn (const char *s)
{
    const char *digits = "0123456789";

    return strspn(s, digits);
}

static int push_numeric_token (CMD *c, const char *tok,
			       const char *s, int pos)
{
    char type = TOK_NUMBER;
    const char *test = tok;

    if (*tok == '+' || *tok == '-') {
	test++;
    }

    if (strlen(test) == digit_spn(test)) {
	type = TOK_INT;
    }

    return push_token(c, tok, s, pos, type, 0);
}

#define ldelim(c) (c == '(' || c == '{' || c == '[')

static int push_delimited_token (CMD *c, const char *tok,
				 const char *s, int pos)
{
    char type = TOK_PRSTR;

    if (*s == '{') {
	type = TOK_CBSTR;
    } else if (*s == '[') {
	type = TOK_BRSTR;
    } else if (*s == '"') {
	type = TOK_STRING;
    }

    return push_token(c, tok, s, pos, type, 0);
}

static int push_quoted_token (CMD *c, const char *s,
			      int len, int pos,
			      int esc)
{
    char *tok = malloc(len + 1);
    int err = 0;

    if (tok == NULL) {
	err = E_ALLOC;
    } else {
	const char *p = s + 1;
	int i = 0;

	if (c->ci == PRINTF || (c->ci == PRINT && !esc)) {
	    /* either a format string, or "print" with no
	       backslash-quote escapes
	    */
	    *tok = '\0';
	    strncat(tok, p, len);
	} else {
	    /* unescape escaped quotes */
	    while (*p) {
		if (*p == '\\' && *(p+1) == '"') {
		    if (*(p+2) == '\0') {
			tok[i++] = *p++;
		    } else {
			tok[i++] = '"';
			p += 2;
		    }
		} else if (*p == '"') {
		    tok[i] = '\0';
		    break;
		} else {
		    tok[i++] = *p++;
		}
	    }
	}

	err = push_token(c, tok, s, pos, TOK_STRING,
			 TOK_QUOTED);
	free(tok);
    }

    return err;
}

static int symbol_spn (const char *s)
{
    const char *ok = "=+-/*<>?|~^!%&.,:;\\'";

    if (*s == '=' && *(s+1) != '=') {
	return 1;
    }

    return strspn(s, ok);
}

/* We'll treat observation identifiers such as "1995:04"
   as numeric in this context, provided the string
   starts with a digit.
*/

static int numeric_spn (const char *s, int digstart)
{
    char *endptr = NULL;
    int n;

    strtod(s, &endptr);
    n = endptr - s;

    if (n > 1 && s[n-1] == '.' && *endptr == '.') {
	/* trailing double-dot */
	return n - 1;
    }

    if (digstart && *endptr == ':') {
	int m = digit_spn(endptr + 1);

	if (m > 0) {
	    n += m + 1;
	}
    }

    return n;
}

static int namechar_spn (const char *s)
{
    const char *ok = "abcdefghijklmnopqrstuvwxyz"
	"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	"0123456789_";

    return strspn(s, ok);
}

static int wild_spn (const char *s)
{
    const char *ok = "abcdefghijklmnopqrstuvwxyz"
	"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	"0123456789_*";

    return strspn(s, ok);
}

/* We allow here for special behavior with PRINT, whereby
   '\"' is interpreted as an escaped double-quote, unless
   it occurs at the end of the string literal. We set *esc
   to non-zero if @s contains any occurrences of the
   special combination of bytes.
*/

static int closing_quote_pos (const char *s, int ci, int *esc)
{
    int i, j, bsl;

    for (i=1; s[i]; i++) {
	if (s[i] == '"' && ci == PRINT) {
	    if (s[i-1] == '\\') {
		if (s[i+1] == '\0') {
		    /* got the closer */
		    return i-1;
		} else {
		    /* the quote is escaped */
		    *esc = 1;
		}
	    } else {
		/* not preceded by backslash */
		return i-1;
	    }
	} else if (s[i] == '"') {
	    /* not plain PRINT */
	    bsl = 0;
	    for (j=i-1; j>=1; j--) {
		if (s[j] == '\\') {
		    bsl++;
		} else {
		    break;
		}
	    }
	    if (bsl % 2 == 0) {
		return i-1;
	    }
	}
    }

    return -1;
}

static int closing_delimiter_pos (const char *s)
{
    const char *ldelims = "({[";
    const char *rdelims = ")}]";
    int ltype = *s;
    int targ = 0;
    int quoted = 0;
    int i, net = 1, n = 0;

    for (i=0; i<3; i++) {
	if (ldelims[i] == ltype) {
	    targ = rdelims[i];
	    break;
	}
    }

    s++;
    while (*s) {
	if (*s == '"' && *(s-1) != '\\') {
	    quoted = !quoted;
	} else if (!quoted) {
	    if (*s == ltype) {
		net++;
	    } else if (*s == targ) {
		net--;
		if (net == 0) {
		    return n;
		}
	    }
	}
	s++;
	n++;
    }

    return -1;
}

/* Determine the index of the first 'real command' token, beyond
   "catch" or assignment to an object.
*/

static int min_token_index (CMD *c, int compmode)
{
    int pos = 0, apos = 1;

    if (c->toks[0].type == TOK_CATCH) {
	/* advance everything by one place */
	pos = 1;
	apos = 2;
    }

    if (c->ntoks > apos && c->toks[apos].type == TOK_ASSIGN) {
	/* advance by assignment target and operator */
	pos += 2;
    } else if (compmode && c->ntoks == 2 && c->toks[0].s[0] == '@') {
	/* got a @-term in pos 0, not a command word */
	pos = 1;
    }

    c->cstart = pos;

    return c->cstart;
}

/* old-style "-%c param" settings, etc. */

static int pseudo_option (const char *s, int ci)
{
    if (ci == EQNPRINT || ci == TABPRINT) {
	if (!strcmp(s, "f")) {
	    return 1;
	}
    } else if (ci == SETINFO) {
	if (!strcmp(s, "d") || !strcmp(s, "n")) {
	    return 1;
	}
    } else if (ci == SMPL || ci == SET) {
	/* short opts confusable with params */
	return 1;
    }

    return 0;
}

/* Find the number of bytes a token takes up in
   the input line, including delimiters if
   applicable.
*/

static int real_toklen (cmd_token *tok)
{
    int n = strlen(tok->s);

    if (delimited_type(tok->type)) {
	n += 2;
    }

    return n;
}

/* Build a string composed of tokens k1 to k2. */

static char *fuse_tokens (CMD *c, int k1, int k2, int n)
{
    char *ret = malloc(n);
    int i;

    if (ret == NULL) {
	c->err = E_ALLOC;
    } else {
	*ret = '\0';
	for (i=k1; i<=k2; i++) {
	    if (c->toks[i].type == TOK_BRSTR) {
		strcat(ret, "[");
		strcat(ret, c->toks[i].s);
		strcat(ret, "]");
	    } else if (c->toks[i].type == TOK_PRSTR) {
		strcat(ret, "(");
		strcat(ret, c->toks[i].s);
		strcat(ret, ")");
	    } else {
		strcat(ret, c->toks[i].s);
	    }
	    mark_token_done(c->toks[i]);
	}
    }

    return ret;
}

static char *merge_option_toks_l_to_r (CMD *c, int k1)
{
    cmd_token *tok = &c->toks[k1];
    int n = real_toklen(tok) + 1;
    int i, k2 = k1;

    for (i=k1+1; i<c->ntoks; i++) {
	tok = &c->toks[i];
	if ((tok->flag & TOK_NOGAP) && !token_done(tok)) {
	    n += real_toklen(tok);
	    k2 = i;
	} else {
	    break;
	}
    }

    return fuse_tokens(c, k1, k2, n);
}

/* Merge tokens from position k1 rightward, stopping at the
   first token that is not left-joined.
*/

static char *merge_toks_l_to_r (CMD *c, int k1)
{
    cmd_token *tok = &c->toks[k1];
    int n = real_toklen(tok) + 1;
    int i, k2 = k1;

    for (i=k1+1; i<c->ntoks; i++) {
	tok = &c->toks[i];
	if (token_joined(tok) && !token_done(tok)) {
	    n += real_toklen(tok);
	    k2 = i;
	} else {
	    break;
	}
    }

    return fuse_tokens(c, k1, k2, n);
}

/* Merge tokens from position k2 leftward, stopping at the
   first token that is not left-joined.
*/

static char *merge_toks_r_to_l (CMD *c, int k2)
{
    cmd_token *prevtok, *tok = &c->toks[k2];
    int n = real_toklen(tok) + 1;
    int i, k1 = k2;

    for (i=k2; i>c->cstart+1; i--) {
	tok = &c->toks[i];
	prevtok = &c->toks[i-1];
	if (token_joined(tok) && !token_done(prevtok)) {
	    n += real_toklen(prevtok);
	    k1 = i - 1;
	} else {
	    break;
	}
    }

    return fuse_tokens(c, k1, k2, n);
}

/* Mark short options (e.g. "-o") and long options (e.g. "--robust"),
   which may be followed by "=".  At this step we're just marking the
   tokens on a purely syntactical basis.
*/

static void mark_option_tokens (CMD *c)
{
    cmd_token *tok, *prevtok;
    int i;

    for (i=1; i<c->ntoks; i++) {
	tok = &c->toks[i];
	if (tok->flag & TOK_NOGAP) {
	    prevtok = &c->toks[i-1];
	    if (!token_joined(prevtok)) {
		if (prevtok->type == TOK_DDASH) {
		    tok->type = TOK_OPT;
		} else if (prevtok->type == TOK_DASH) {
		    if (!pseudo_option(tok->s, c->ci)) {
			tok->type = TOK_SOPT;
		    }
		}
		if (tok->type == TOK_OPT || tok->type == TOK_SOPT) {
		    prevtok->type = TOK_OPTDASH;
		    prevtok->flag |= TOK_DONE;
		}
	    } else if (prevtok->type == TOK_OPT) {
		if (tok->type == TOK_EQUALS) {
		    tok->type = TOK_OPTEQ;
		} else if (tok->type == TOK_DASH ||
			   tok->type == TOK_NAME ||
			   tok->type == TOK_INT) {
		    tok->type = TOK_OPT; /* continuation of option flag */
		}
	    } else if (prevtok->type == TOK_OPTEQ) {
		if (tok->type == TOK_NAME || tok->type == TOK_STRING) {
		    tok->type = TOK_OPTVAL;
		}
	    } else if (prevtok->type == TOK_OPTVAL) {
		if (tok->type == TOK_NAME || tok->type == TOK_DOT) {
		    tok->type = TOK_OPTVAL;
		}
	    }
	}
    }
}

static cmd_token *next_joined_token (CMD *c, int i)
{
    if (i >= 0 && i < c->ntoks - 1) {
	if (c->toks[i+1].flag & TOK_JOINED) {
	    return &c->toks[i+1];
	}
    }

    return NULL;
}

static cmd_token *next_nogap_token (CMD *c, int i)
{
    if (i >= 0 && i < c->ntoks - 1) {
	if (c->toks[i+1].flag & TOK_NOGAP) {
	    return &c->toks[i+1];
	}
    }

    return NULL;
}

static int handle_option_value (CMD *c, int i, gretlopt opt,
				OptStatus status)
{
    cmd_token *tok = &c->toks[i];
    cmd_token *nexttok = next_joined_token(c, i);
    int getval = 0;

    if (nexttok != NULL) {
	if (status == OPT_NO_PARM) {
	    /* nothing should be glued onto the option flag */
	    c->err = E_PARSE;
	} else if (nexttok->type == TOK_OPTEQ) {
	    nexttok->flag |= TOK_DONE;
	    getval = 1;
	} else {
	    /* the only acceptable following field is '=' */
	    c->err = E_PARSE;
	}
    }

    if (status == OPT_NEEDS_PARM && !getval && !c->err) {
	/* we need but didn't get a value */
	gretl_errmsg_sprintf(_("The option '--%s' requires a parameter"),
			     tok->s);
	c->err = E_BADOPT;
    }

    if (getval) {
	char *val = NULL;

	nexttok = next_nogap_token(c, i + 1);
	if (nexttok == NULL) {
	    c->err = E_PARSE;
	} else if (delimited_type(nexttok->type)) {
	    val = gretl_strdup(nexttok->s);
	    nexttok->flag |= TOK_DONE;
	} else {
	    val = merge_option_toks_l_to_r(c, i + 2);
	}

#if CDEBUG > 1
	fprintf(stderr, "option '--%s': param='%s'\n", tok->s,
		(val == NULL)? "NULL" : val);
#endif

	if (val != NULL) {
	    c->err = push_option_param(c->ci, opt, val);
	    if (c->err) {
		free(val);
	    }
	} else if (!c->err) {
	    c->err = E_ALLOC;
	}
    }

    return c->err;
}

#define OPTLEN 32

static int assemble_option_flag (CMD *c, cmd_token *tok,
				 char *flag, int *pk,
				 int dryrun)
{
    int i, n = 0, added = 0;

    *flag = '\0';

    for (i=*pk; i<c->ntoks; i++) {
	tok = &c->toks[i];
	if (tok->type == TOK_OPT) {
	    n += strlen(tok->s);
	    if (n >= OPTLEN) {
		fprintf(stderr, "option string too long\n");
		c->err = E_PARSE;
		break;
	    } else {
		if (!dryrun) {
		    tok->flag |= TOK_DONE;
		}
		strcat(flag, tok->s);
		if (added) {
		    *pk += 1;
		} else {
		    added = 1;
		}
	    }
	} else {
	    break;
	}
    }

#if CDEBUG > 2
    fprintf(stderr, "option flag: '%s'\n", flag);
#endif

    return c->err;
}

/* For now, handle the old fit-type options whose strings
   are listed below. If would be nice to get rid of this
   altogether (in favour of the newer --fit=whatever
   option).
*/

static void handle_legacy_gnuplot_options (CMD *c)
{
    /* these used to be options in their own right */
    const char *old_opts[] = {
	"inverse-fit",
	"loess-fit",
	"quadratic-fit",
	"linear-fit",
	"semilog-fit",
	"suppress-fitted"
    };
    /* but are replaced by these flags for --fit */
    const char *repl[] = {
	"inverse",
	"loess",
	"quadratic",
	"linear",
	"semilog",
	"none"
    };
    char optflag[OPTLEN];
    cmd_token *tok;
    int i, j, n, pos, len;
    int err, done = 0;

    n = G_N_ELEMENTS(old_opts);

    for (i=1; i<c->ntoks && !done; i++) {
	tok = &c->toks[i];
	if (token_done(tok)) {
	    continue;
	}
	if (tok->type == TOK_OPT) {
	    pos = i;
	    err = assemble_option_flag(c, tok, optflag, &i, 1);
	    if (err) {
		break;
	    }
	    for (j=0; j<n && !done; j++) {
		len = strlen(optflag);
		if (len > 2 && !strncmp(optflag, old_opts[j], len)) {
		    /* found an obsolete fit option */
		    c->opt |= OPT_F;
		    set_optval_string(GNUPLOT, OPT_F, repl[j]);
		    for (j=pos; j<=i; j++) {
			c->toks[j].flag |= TOK_DONE;
		    }
		    done = 1;
		}
	    }
	}
    }
}

/* handle --remove and --purge, which used to be options
   to the defunct "install" command
*/

static void handle_legacy_install_options (CMD *c)
{
    char optflag[OPTLEN];
    cmd_token *tok;
    int i, j, pos;
    int err, done = 0;

    for (i=1; i<c->ntoks && !done; i++) {
	tok = &c->toks[i];
	if (token_done(tok)) {
	    continue;
	}
	if (tok->type == TOK_OPT) {
	    pos = i;
	    err = assemble_option_flag(c, tok, optflag, &i, 1);
	    if (err) {
		break;
	    }
	    if (!strcmp(optflag, "remove")) {
		c->opt |= OPT_R;
		done = 1;
	    } else if (!strcmp(optflag, "purge")) {
		c->opt |= OPT_P;
		done = 1;
	    }
	    if (done) {
		for (j=pos; j<=i; j++) {
		    c->toks[j].flag |= TOK_DONE;
		}
	    }
	}
    }
}

static int check_command_options (CMD *c)
{
    cmd_token *tok;
    char optflag[OPTLEN];
    OptStatus status;
    gretlopt opt;
    int save_ci = c->ci;
    int i, j, n;
    int err = 0;

    mark_option_tokens(c);

    if (c->ci == END && c->context) {
	c->ci = c->context;
    } else if (c->ci == SETOPT) {
	c->ci = gretl_command_number(c->param);
    } else if (c->ci == GNUPLOT) {
	handle_legacy_gnuplot_options(c);
    } else if (c->ci == PKG && (c->opt & OPT_B)) {
	handle_legacy_install_options(c);
    } else if (c->ci == END && c->ntoks > 1) {
        if (gretl_command_number(c->toks[1].s) == GPBUILD) {
            /* support the --output option */
            c->ci = GPBUILD;
        }
    }

    for (i=1; i<c->ntoks && !err; i++) {
	tok = &c->toks[i];
	if (token_done(tok)) {
	    continue;
	}
	if (tok->type == TOK_OPT) {
	    /* long-form option, possibly with attached value */
	    err = assemble_option_flag(c, tok, optflag, &i, 0);
	    if (!err) {
		opt = valid_long_opt(c->ci, optflag, &status);
		if (opt == OPT_NONE) {
		    gretl_errmsg_sprintf(_("Invalid option '--%s'"), optflag);
		    err = E_BADOPT;
		} else {
		    err = handle_option_value(c, i, opt, status);
		    if (!err) {
			c->opt |= opt;
		    }
		}
	    }
	} else if (tok->type == TOK_SOPT) {
	    /* short-form option(s) */
	    tok->flag |= TOK_DONE;
	    n = strlen(tok->s);
	    for (j=0; j<n && !err; j++) {
		opt = valid_short_opt(c->ci, tok->s[j]);
		if (opt == OPT_NONE) {
		    /* error message handled in options.c */
		    err = E_BADOPT;
		} else {
		    c->opt |= opt;
		}
	    }
	}
    }

    if (!c->err && save_ci != SETOPT) {
	/* Retrieve any options put in place via "setopt" */
	maybe_get_stored_options(c->ci, &c->opt);
    }

    c->ci = save_ci;

    return err;
}

/* Get the 0-based token position of the 'real' command argument
   in 1-based position @k (i.e. k = 1 gets first real arg),
   or -1 if none.
*/

static int real_arg_index (CMD *c, int k)
{
    int i, i0 = c->cstart + k;

    for (i=i0; i<c->ntoks; i++) {
	if (!(c->toks[i].flag & TOK_DONE)) {
	    return i;
	}
    }

    return -1;
}

/* Get the 0-based token position of the last command token
   that is not already handled (or -1 if none).
*/

static int last_arg_index (CMD *c)
{
    int i;

    for (i=c->ntoks-1; i>c->cstart; i--) {
	if (!(c->toks[i].flag & TOK_DONE)) {
	    return i;
	}
    }

    return -1;
}

/* Get the 0-based token position of the first command token
   that is not already handled (or -1 if none).
*/

static int first_unused_arg_index (CMD *c)
{
    int i;

    for (i=c->cstart; i<c->ntoks; i++) {
	if (!(c->toks[i].flag & TOK_DONE)) {
	    return i;
	}
    }

    return -1;
}

/* Record the content of the token at @pos as either
   param or parm2, depending on @i.
*/

static int token_to_param (CMD *c, int pos, int i)
{
    char *s = c->toks[pos].s;

    if (i == 1) {
	c->param = s;
    } else if (i == 2) {
	c->parm2 = s;
    } else {
	c->err = E_PARSE;
    }

    return c->err;
}

static int may_be_dash_parm (cmd_token *tok, CMD *c)
{
    if (tok->type == TOK_NAME || tok->type == TOK_STRING) {
	return 1;
    } else if ((c->ci == TABPRINT || c->ci == EQNPRINT) &&
	       tok->type == TOK_SYMB && !strcmp(tok->s, "/")) {
	return 1;
    } else {
	return 0;
    }
}

/* legacy: look for, e.g., "-f filename" */

static int dash_char_index (CMD *c, const char *s)
{
    cmd_token *tok;
    int step = 0;
    int i;

    for (i=c->cstart+1; i<c->ntoks; i++) {
	tok = &c->toks[i];
	if (step == 0) {
	    /* dash, freestanding to left */
	    if (!token_done(tok) && tok->type == TOK_DASH &&
		!token_joined(tok)) {
		step = 1;
	    }
	} else if (step == 1) {
	    /* suitable char token stuck to dash */
	    if (!token_done(tok) && tok->type == TOK_NAME &&
		token_joined(tok) && !strcmp(tok->s, s)) {
		step = 2;
	    } else {
		step = 0;
	    }
	} else if (step == 2) {
	    /* suitable following string token */
	    if (!token_done(tok) && !token_joined(tok) &&
		may_be_dash_parm(tok, c)) {
		mark_token_done(c->toks[i-1]);
		mark_token_done(c->toks[i-2]);
		return i;
	    } else {
		step = 0;
	    }
	}
    }

    return -1;
}

/* handle TABPRINT and EQNPRINT */

static int legacy_get_filename (CMD *c)
{
    int pos = dash_char_index(c, "f");

    if (pos > 0) {
	if (c->toks[pos].type == TOK_STRING) {
	    c->param = c->toks[pos].s;
	    mark_token_done(c->toks[pos]);
	} else if (next_joined_token(c, pos) == NULL) {
	    c->param = c->toks[pos].s;
	    mark_token_done(c->toks[pos]);
	} else {
	    c->param = merge_toks_l_to_r(c, pos);
	}
    }

    return c->err;
}

static int handle_datamod_param (CMD *c)
{
    int op = dataset_op_from_string(c->param);

    if (op == DS_NONE) {
	c->err = E_PARSE;
    } else {
	if (op == DS_SORTBY || op == DS_DSORTBY) {
	    /* we now need a list, only */
	    c->ciflags = CI_PARM1 | CI_LIST;
	} else if (op == DS_RENUMBER) {
	    /* we need a one-member list plus param */
	    c->ciflags = CI_PARM1 | CI_LIST | CI_LLEN1 | CI_PARM2;
	} else if (op == DS_COMPACT) {
	    /* "dataset compact": the second param may contain
	       two fields, as in "4 last", so parse it as a
	       case of 'extra'
	    */
	    c->ciflags = CI_PARM1 | CI_EXTRA;
	} else if (op == DS_TRANSPOSE || op == DS_CLEAR ||
		   op == DS_UNPAD_DAILY) {
	    /* no more fields wanted */
	    c->ciflags = CI_PARM1;
	} else {
	    /* all other cases: no list wanted */
	    c->ciflags &= ~CI_LIST;
	}
    }

    if (!c->err) {
	c->auxint = op;
    }

    return c->err;
}

static char *rebrace_string (const char *s, int *err)
{
    char *ret = malloc(strlen(s) + 3);

    if (ret == NULL) {
	*err = E_ALLOC;
    } else {
	sprintf(ret, "{%s}", s);
    }

    return ret;
}

static int
looks_like_list_token (CMD *c, cmd_token *tok, const DATASET *dset)
{
    if (tok->type != TOK_QUOTED) {
	/* heuristic for a token that forms part of a series
	   list: starts with a digit (series ID?), contains
	   '*' (wildcard spec), or is the name of a series
	*/
	if (isdigit(*tok->s) ||
	    strchr(tok->s, '*') != NULL ||
	    current_series_index(dset, tok->s) >= 0) {
	    return 1;
	} else if (get_list_by_name(tok->s)) {
	    /* should be treated as list, not param, iff
	       this is a genuine "delete" command, not
	       "list foo delete"
	    */
	    if (!strcmp(c->toks[c->cstart].s, "delete")) {
		return 1;
	    }
	}
    }

    return 0;
}

/* Get command parameter in first position; may involve
   compositing tokens.
*/

static int get_param (CMD *c, const DATASET *dset)
{
    int pos = real_arg_index(c, 1);
    cmd_token *tok;

    if (pos < 0) {
	if (!param_optional(c->ci)) {
	    c->err = E_ARGS;
	    fprintf(stderr, "%s: required param is missing\n",
		    c->toks[c->cstart].s);
	} else if (c->ci == SMPL) {
	    /* allow the null form of "smpl" */
	    c->ciflags &= ~CI_PARM2;
	}
	return c->err;
    }

    tok = &c->toks[pos];

    if (c->ci == DELEET && !(c->opt & OPT_L)) {
	/* experimental */
	if (looks_like_list_token(c, tok, dset)) {
	    c->ciflags &= ~CI_PARM1;
	    c->ciflags |= CI_LIST;
	    return 0;
	}
    }

    if (tok->type == TOK_CBSTR) {
	/* if param was found in braces, it should
	   probably be passed in braces, but FIXME
	   check for exceptions?
	*/
	c->param = rebrace_string(tok->s, &c->err);
	mark_token_done(c->toks[pos]);
    } else if (tok->type == TOK_EVAL) {
	c->param = generate_string(tok->s, NULL, &c->err);
	mark_token_done(c->toks[pos]);
    } else if (delimited_type(tok->type)) {
	c->param = tok->s;
	mark_token_done(c->toks[pos]);
    } else if (next_joined_token(c, pos) == NULL) {
	c->param = tok->s;
	mark_token_done(c->toks[pos]);
    } else {
	c->param = merge_toks_l_to_r(c, pos);
    }

    if (c->ci == DATAMOD) {
	/* the "dataset" command */
	handle_datamod_param(c);
    } else if (c->ci == SMPL) {
	/* "smpl" command: drop the requirement for a second param
	   if the first is "full" */
	if (!strcmp(c->param, "full")) {
	    c->opt |= OPT_F;
	    c->ciflags ^= CI_PARM2;
	}
    } else if (c->ci == HELP) {
	/* allow a second param for "help set ..." */
	if (!strcmp(c->param, "set")) {
	    c->ciflags |= CI_PARM2;
	}
    }

    return c->err;
}

static int pkg_params_compat (CMD *c)
{
    /* swap pkgname into second position */
    c->parm2 = c->param;

    if (c->opt & OPT_R) {
	/* compat for old --remove */
	c->param = gretl_strdup("unload");
	c->opt ^= OPT_R;
    } else if (c->opt & OPT_P) {
	/* compat for old --purge */
	c->param = gretl_strdup("remove");
	c->opt ^= OPT_P;
    } else {
	/* implicit alternative */
	c->param = gretl_strdup("install");
    }

    return 0;
}

/* Get command parameter in last position; may involve
   compositing tokens.
*/

static int get_parm2 (CMD *c, int options_later)
{
    cmd_token *tok;
    int pos;

    pos = options_later ? first_unused_arg_index(c) :
	last_arg_index(c);

    if (pos < 0) {
	if (c->ci == SMPL && c->opt == 0) {
	    /* backward-compatible slop factor: allow missing ';'
	       in second place?
	    */
	    ;
	} else if (c->ci == PKG && (c->opt & OPT_B)) {
	    /* "install" emulating "pkg" */
	    c->err = pkg_params_compat(c);
	    c->opt ^= OPT_B; /* scrub temporary option */
	} else if (c->ci == DATAMOD && c->auxint == DS_EXPAND) {
	    ; /* OK, parm2 is optional in this case */
	} else if (!parm2_optional(c->ci)) {
	    c->err = E_ARGS;
	    fprintf(stderr, "%s: required parm2 is missing\n",
		    c->toks[c->cstart].s);
	}
	return c->err;
    }

    tok = &c->toks[pos];

    if (c->ciflags & CI_VARGS) {
	/* check for trailing comma: if it's present, the
	   token we want will precede it */
	if (tok->type == TOK_COMMA) {
	    tok->flag |= TOK_DONE;
	    pos--;
	} else if (vargs_optional(c->ci)) {
	    c->ciflags ^= CI_VARGS;
	} else {
	    c->err = E_PARSE;
	    return c->err;
	}
    }

    tok = &c->toks[pos];

    /* revised 2014-10-11 */
    if (token_joined(tok)) {
	c->parm2 = merge_toks_r_to_l(c, pos);
    } else {
	c->parm2 = tok->s;
	if ((tok->flag & TOK_QUOTED) && c->ci == MODPRINT) {
	    /* signal quoted status of string literal */
	    c->opt |= OPT_L;
	}
	tok->flag |= TOK_DONE;
    }

    return c->err;
}

/* legacy: handle SETINFO fields */

static int get_quoted_dash_fields (CMD *c, const char *s)
{
    char test[2] = {0};
    int pos = real_arg_index(c, 2);
    int i;

    if (pos < 0) {
	return 0;
    }

    for (i=0; s[i]; i++) {
	/* loop across "dash fields" */
	test[0] = s[i];
	pos = dash_char_index(c, test);
	if (pos > 0) {
	    if (c->toks[pos].type == TOK_STRING) {
		token_to_param(c, pos, i+1);
		mark_token_done(c->toks[pos]);
	    } else {
		c->err = E_PARSE;
	    }
	}
    }

    return c->err;
}

static int first_arg_quoted (CMD *c)
{
    int pos = real_arg_index(c, 1);

    if (pos < 0) {
	return 0;
    } else {
	return c->toks[pos].type == TOK_STRING;
    }
}

/* Count instances of list separator, ';', in the
   command line.
*/

static int cmd_get_sepcount (CMD *c)
{
    int i, n = 0;

    for (i=c->cstart+1; i<c->ntoks; i++) {
	if (c->toks[i].type == TOK_SEMIC) {
	    n++;
	}
    }

    return n;
}

/* Convert token @k to an integer. */

static int token_to_int (CMD *c, int k)
{
    cmd_token *tok = &c->toks[k];
    int ret = -1;

    if (tok->type == TOK_INT) {
	ret = atoi(tok->s);
	tok->flag |= TOK_DONE;
    } else if (tok->type == TOK_NAME) {
	double x = get_scalar_value_by_name(tok->s, &c->err);

	if (!c->err) {
	    if (x > 0 && x < INT_MAX) {
		tok->flag |= TOK_DONE;
		ret = x;
	    } else {
		c->err = E_INVARG;
	    }
	}
    } else {
	c->err = E_INVARG;
    }

    return ret;
}

static int list_max (const int *list)
{
    int i, lmax = list[1];

    for (i=2; i<=list[0]; i++) {
	if (list[i] > lmax) {
	    lmax = list[i];
	}
    }

    return lmax;
}

static int *get_auxlist (cmd_token *tok, int *err)
{
    int *alist = NULL;

    if (tok->type == TOK_CBSTR) {
	alist = gretl_list_from_string(tok->s, err);
    } else {
	gretl_matrix *m = get_matrix_by_name(tok->s);

	alist = gretl_auxlist_from_vector(m, err);
    }

    return alist;
}

static int get_VAR_order (CMD *c, int k)
{
    cmd_token *tok = &c->toks[k];
    int ret = -1;

    if (tok->type == TOK_INT) {
	ret = atoi(tok->s);
    } else if (tok->type == TOK_NAME) {
	double x = get_scalar_value_by_name(tok->s, &c->err);

	if (!c->err) {
	    if (x > 0 && x < INT_MAX) {
		ret = x;
	    } else {
		gretl_errmsg_sprintf(_("Invalid lag order %g"), x);
		c->err = E_INVARG;
	    }
	} else {
	    /* could be a non-scalar matrix */
	    c->err = gretl_error_clear();
	    c->auxlist = get_auxlist(tok, &c->err);
	}
    } else if (tok->type == TOK_CBSTR) {
	c->auxlist = get_auxlist(tok, &c->err);
    } else {
	c->err = E_INVARG;
    }

    if (!c->err) {
	tok->flag |= TOK_DONE;
	if (c->auxlist != NULL) {
	    ret = list_max(c->auxlist);
	}
    }

    return ret;
}

static int get_bundled_int (CMD *c, int k)
{
    char *s = merge_toks_l_to_r(c, k);
    int ret = -1;

    if (s != NULL) {
	double x = get_scalar_value_by_name(s, &c->err);

	if (!c->err) {
	    if (x > 0 && x < INT_MAX) {
		ret = x;
	    } else {
		c->err = E_INVARG;
	    }
	}
	free(s);
    }

    return ret;
}

/* Some commands require an integer order as the first
   argument, but we also have the cases where an order
   in first argument position is optional, viz:

   var order ...
   modtest [ order ]
   lags [ order ; ] list
*/

static int get_command_order (CMD *c)
{
    int pos = real_arg_index(c, 1);
    int try_bundle = 0;

    if (pos >= 0) {
	cmd_token *ntok = next_joined_token(c, pos);

	/* In general, nothing should be "stuck onto" the
	   order specifier; however, we have accepted ";"
	   without an intervening space, and we might be
	   given a bundle member in dot notation.
	*/
	if (ntok != NULL) {
	    if (ntok->type == TOK_SEMIC) {
		; /* OK */
	    } else if (ntok->type == TOK_DOT) {
		try_bundle = 1; /* try handling this */
	    } else {
		c->err = E_PARSE;
		return c->err;
	    }
	}
    }

    if (c->ci == MODTEST && pos < 0) {
	/* order is optional, not present, OK */
	return 0;
    }

    if (c->ci == LAGS && cmd_get_sepcount(c) == 0) {
	/* order is optional, not present, OK */
	return 0;
    }

    if (pos < 0) {
	c->err = E_ARGS;
    } else if (try_bundle) {
	c->order = get_bundled_int(c, pos);
    } else if (c->ci == VAR) {
	/* order can be special, "gappy" */
	c->order = get_VAR_order(c, pos);
    } else {
	c->order = token_to_int(c, pos);
    }

    return c->err;
}

/* special for VECM only, so far */

static int get_vecm_rank (CMD *c)
{
    int pos = real_arg_index(c, 2);

    if (pos < 0) {
	c->err = E_ARGS;
    } else {
	cmd_token *nt = next_joined_token(c, pos);

	if (nt != NULL && nt->type == TOK_DOT) {
	    c->auxint = get_bundled_int(c, pos);
	} else {
	    c->auxint = token_to_int(c, pos);
	}
    }

    return c->err;
}

/* get an optional, trailing "order" field */

static int get_optional_order (CMD *c)
{
    int pos = first_unused_arg_index(c);

    if (pos > 0) {
	cmd_token *nt = next_joined_token(c, pos);

	if (nt != NULL && nt->type == TOK_DOT) {
	    c->order = get_bundled_int(c, pos);
	} else {
	    c->order = token_to_int(c, pos);
	}
    }

    return c->err;
}

/* Stuff that can come before a command proper: the "catch"
   keyword or assignment to a named object.
*/

static int handle_command_preamble (CMD *c)
{
    int pos = 0;

    if (c->toks[0].type == TOK_CATCH) {
	if (not_catchable(c->ci)) {
	    gretl_errmsg_set(_("catch: cannot be applied to this command"));
	    c->err = E_DATA;
	    return c->err;
	} else {
	    set_gretl_errno(0);
	    gretl_error_clear();
	    c->flags |= CMD_CATCH;
	    mark_token_done(c->toks[0]);
	    pos = 1;
	}
    }

    if (c->ntoks > pos + 1 && c->toks[pos+1].type == TOK_ASSIGN) {
	char *s = c->toks[pos].s;
	int n = strlen(s);

	if (n >= MAXSAVENAME) {
	    gretl_errmsg_set(_("savename is too long"));
	    c->err = E_DATA;
	} else {
	    strcpy(c->savename, s);
	    mark_token_done(c->toks[pos]);
	    mark_token_done(c->toks[pos+1]);
	}
    }

    return c->err;
}

#if CDEBUG

static void vstart_line_out (CMD *c)
{
    if (c->ciflags & CI_EXPR) {
	fprintf(stderr, "* expr: '%s'\n", c->vstart);
    } else if (c->ciflags & CI_VARGS) {
	fprintf(stderr, "* varargs: '%s'\n", c->vstart);
    }
}

static void my_printlist (const int *list, const char *s)
{
    int i;

    if (list == NULL) {
	fprintf(stderr, "%s is NULL\n", s);
	return;
    }

    fprintf(stderr, "%s: %d : ", s, list[0]);

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    fputs("; ", stderr);
	} else {
	    fprintf(stderr, "%d ", list[i]);
	}
    }

    fputc('\n', stderr);
}

# if CDEBUG > 1

static char *tokstring (char *s, cmd_token *toks, int i)
{
    cmd_token *tok = &toks[i];

    *s = '\0';

    if (tok->type == TOK_OPT) {
	strcpy(s, "option");
    } else if (tok->type == TOK_SOPT) {
	strcpy(s, "short-option");
    } else if (tok->type == TOK_STRING) {
	strcpy(s, "quoted");
    } else if (tok->type == TOK_NUMBER) {
	strcpy(s, "number");
    } else if (tok->type == TOK_INT) {
	strcpy(s, "integer");
    } else if (tok->type == TOK_OPTDASH) {
	strcpy(s, "option-leader");
    } else if (tok->type == TOK_DASH) {
	strcpy(s, "dash");
    } else if (tok->type == TOK_DOT) {
	strcpy(s, "dot");
    } else if (tok->type == TOK_COMMA) {
	strcpy(s, "comma");
    } else if (tok->type == TOK_SEMIC) {
	strcpy(s, "separator");
    } else if (tok->type == TOK_COLON) {
	strcpy(s, "colon");
    } else if (tok->type == TOK_DDASH) {
	strcpy(s, "double-dash");
    } else if (tok->type == TOK_DPLUS) {
	strcpy(s, "double-plus");
    } else if (tok->type == TOK_DDOT) {
	strcpy(s, "double-dot");
    } else if (tok->type == TOK_EQUALS) {
	strcpy(s, "equals");
    } else if (tok->type == TOK_EQMOD) {
	strcpy(s, "modified-equals");
    } else if (tok->type == TOK_OPTEQ) {
	strcpy(s, "opt-equals");
    } else if (tok->type == TOK_OPTVAL) {
	if (tok->flag & TOK_QUOTED) {
	    strcpy(s, "quoted option-value");
	} else {
	    strcpy(s, "option-value");
	}
    } else if (tok->type == TOK_ASSIGN) {
	strcpy(s, "assign");
    } else if (tok->type == TOK_SYMB) {
	strcpy(s, "symbol/operator");
    } else if (tok->type == TOK_PRSTR) {
	strcpy(s, "paren-delimited");
    } else if (tok->type == TOK_CBSTR) {
	strcpy(s, "brace-delimited");
    } else if (tok->type == TOK_BRSTR) {
	strcpy(s, "bracket-delimited");
    } else if (tok->type == TOK_DOLSTR) {
	strcpy(s, "$-variable");
    } else if (tok->type == TOK_EVAL) {
	strcpy(s, "eval-token");
    } else {
	strcpy(s, "regular");
    }

    if (tok->type != TOK_OPT && tok->type != TOK_SOPT &&
	tok->type != TOK_OPTEQ && tok->type != TOK_OPTVAL &&
	token_joined(tok)) {
	strcat(s, ", joined on left");
    }

    if (token_ignored(tok)) {
	strcat(s, ", ignored!");
    } else if (token_done(tok)) {
	strcat(s, ", handled");
    } else {
	strcat(s, ", not handled");
    }

    return s;
}

# endif

static void print_option_flags (gretlopt opt)
{
    const char *flags = "ABCDEFGHIJKLMNOPQRSTUVWXYZbi";
    gretlopt i, n = 0;

    for (i=OPT_A; i<=OPT_i; i*=2) {
	if (opt & i) {
	    if (n > 0) {
		fputc('|', stderr);
	    }
	    fprintf(stderr, "OPT_%c", *flags);
	    n++;
	}
	flags++;
    }
}

static void print_tokens (CMD *c)
{
    const char *inh = "";
    cmd_token *toks = c->toks;
    int i;

# if CDEBUG > 1
    int nt = c->ntoks;
    char desc[48];

    if (c->ciflags & (CI_EXPR | CI_ADHOC)) {
	fprintf(stderr, "tokens examined so far:\n");
    } else {
	fprintf(stderr, "tokens:\n");
    }

    for (i=0; i<nt; i++) {
	fprintf(stderr, "%3d: '%s' (%s)\n", i, toks[i].s,
		tokstring(desc, toks, i));
    }
# endif

    if (c->context || c->ci == END) {
	inh = " (inherited)";
    };

    if (c->flags & CMD_CATCH) {
	fprintf(stderr, "* catching errors%s\n", inh);
    } else if (*c->savename != '\0') {
	fprintf(stderr, "* assignment to '%s'%s\n",
		c->savename, inh);
    }

    i = c->cstart;

    if (c->ci > 0) {
	if (c->context) {
	    fprintf(stderr, "* command: '%s' (%d), context = '%s' (%d)\n",
		    gretl_command_word(c->ci), c->ci,
		    gretl_command_word(c->context), c->context);
	} else {
	    fprintf(stderr, "* command: '%s' (%d), context = null\n",
		    gretl_command_word(c->ci), c->ci);
	}
    } else {
	fprintf(stderr, "* '%s': don't understand\n", toks[i].s);
    }

    if (c->err == E_BADOPT) {
	fprintf(stderr, "* command contains invalid option\n");
    } else if (c->opt) {
	fprintf(stderr, "* command option(s) valid (");
	print_option_flags(c->opt);
	fprintf(stderr, ")\n");
    }

    if (c->param != NULL) {
	fprintf(stderr, "* param: '%s'\n", c->param);
    }

    if (c->parm2 != NULL) {
	fprintf(stderr, "* parm2: '%s'\n", c->parm2);
    }

    if (c->order != 0) {
	fprintf(stderr, "* order = %d\n", c->order);
    }

    if (c->auxint != 0) {
	fprintf(stderr, "* auxint = %d\n", c->auxint);
    }

    if (c->err) {
	fprintf(stderr, "*** error = %d (%s)\n\n", c->err,
		errmsg_get_with_default(c->err));
	return;
    }

    if (c->ciflags & CI_LIST) {
	my_printlist(c->list, "list");
    } else if (c->vstart != NULL) {
	vstart_line_out(c);
    }

    if (c->auxlist != NULL) {
	my_printlist(c->auxlist, "auxlist");
    }

    fputc('\n', stderr);
}

#endif /* CDEBUG */

static int stray_symbol_error (cmd_token *tok)
{
    GretlType t = user_var_get_type_by_name(tok->s);

    if (t == GRETL_TYPE_NONE) {
	gretl_errmsg_sprintf(_("Parse error at unexpected token '%s'"),
			     tok->s);
	return E_PARSE;
    } else {
	gretl_errmsg_sprintf(_("The variable %s is of type %s, "
			       "not acceptable in context"),
			     tok->s, gretl_type_get_name(t));
	return E_TYPES;
    }
}

static int check_for_stray_tokens (CMD *c)
{
    if (!(c->ciflags & (CI_EXPR | CI_ADHOC))) {
	cmd_token *tok;
	int i;

	for (i=0; i<c->ntoks && !c->err; i++) {
	    tok = &c->toks[i];
	    if (!token_done(tok)) {
		c->err = stray_symbol_error(tok);
	    }
	}
    }

    return c->err;
}

static int check_list_sepcount (int ci, int nsep)
{
    int minsep, maxsep;
    int err = 0;

    maxsep = get_sep_max(ci, &minsep);

    if (nsep < minsep) {
	err = E_ARGS;
    } else if (nsep > maxsep) {
	gretl_errmsg_sprintf(_("The symbol '%c' is not valid in this context\n"),
			     ';');
	err = E_INVARG;
    }

#if CDEBUG
    if (err) {
	fprintf(stderr, "error %d from check_list_sepcount\n", err);
    }
#endif

    return err;
}

/* Determine if the command-line is a "genr"-type expression,
   which will be directed to a separate parser -- this gets
   invoked if we haven't been able to find a recognizable
   command-word.
*/

static int test_for_genr (CMD *c, int i, char cnext, DATASET *dset)
{
    cmd_token *toks = c->toks;
    char *s = toks[i].s;
    int j = c->ntoks - 1;

    if (j > i && (toks[j].type == TOK_EQUALS || toks[j].type == TOK_EQMOD)) {
	c->ci = GENR;
    } else if (dset != NULL && cnext != '(' && current_series_index(dset, s) >= 0) {
	c->ci = GENR;
    } else if (cnext != '(' && gretl_is_user_var(s)) {
	c->ci = GENR;
    } else if (toks[i].type == TOK_NAME && c->ntoks > i + 1) {
	cmd_token *nexttok = &toks[i+1];

	if ((nexttok->type == TOK_DPLUS ||
	     nexttok->type == TOK_DDASH) &&
	    token_joined(nexttok)) {
	    /* increment/decrement */
	    c->ci = GENR;
	} else if (nexttok->type == TOK_EQUALS || nexttok->type == TOK_EQMOD) {
	    /* assignment token in second place */
	    c->ci = GENR;
	} else if (nexttok->type == TOK_BRSTR && token_joined(nexttok)) {
	    /* assignment to array or series element(s) */
	    c->ci = GENR;
	} else if (function_lookup(s) || get_user_function_by_name(s) ||
		   is_function_alias(s)) {
	    /* function call, no assignment */
	    c->ci = GENR;
	    c->opt |= OPT_O;
	}
#ifdef USE_RLIB
	else if (j == 1 && toks[j].type == TOK_DOT &&
		 libset_get_bool(R_FUNCTIONS) &&
		 toks[0].type == TOK_NAME && !strcmp(toks[0].s, "R")) {
	    /* R function call, no assignment */
	    c->ci = GENR;
	    c->opt |= OPT_O;
	}
#endif
    }

    return c->ci;
}

static void set_deprecation (const char *bad, const char *good,
			     int command)
{
    const char *tag = command ? "command" : "construction";

    if (strstr(good, "()")) {
	gretl_warnmsg_sprintf(_("\"%s\": obsolete %s; please use the function %s"),
			      bad, tag, good);
    } else {
	gretl_warnmsg_sprintf(_("\"%s\": obsolete %s; please use \"%s\""),
			      bad, tag, good);
    }
}

static int try_for_command_alias (const char *s, CMD *cmd)
{
    int ci = 0;

    if (!strcmp(s, "exit")) {
	ci = QUIT;
	cmd->opt = OPT_X;
    } else if (!strcmp(s, "ls")) {
	ci = VARLIST;
    } else if (!strcmp(s, "tsplots")) {
	ci = SCATTERS;
	cmd->opt |= OPT_T;
    } else if (!strcmp(s, "pooled")) {
	set_deprecation("pooled", "ols", 1);
	ci = OLS;
    } else if (!strcmp(s, "equations")) {
	/* reached only when compiling loop */
	ci = EQUATION;
	cmd->opt |= OPT_M;
    } else if (*s == '!') {
	ci = SHELL;
    } else if (!strcmp(s, "launch")) {
	ci = SHELL;
	cmd->opt |= OPT_A;
    } else if (!strcmp(s, "fcasterr")) {
	set_deprecation("fcasterr", "fcast", 1);
	ci = FCAST;
    } else if (!strcmp(s, "install")) {
	ci = PKG;
	cmd->opt |= OPT_B; /* back-compat */
    }

    return ci;
}

static char peek_next_char (CMD *cmd, int i)
{
    const char *s;

    s = cmd->toks[i].lp + strlen(cmd->toks[i].s);
    s += strspn(s, " ");
    return *s;
}

static int peek_end_param (CMD *cmd, int i)
{
    const char *s;

    s = cmd->toks[i].lp + strlen(cmd->toks[i].s);
    s += strspn(s, " ");

    if (!strncmp(s, "loop", 4)) {
	return LOOP;
    } else if (!strncmp(s, "function", 8)) {
	return FUNC;
    } else {
	return 0;
    }
}

#if CDEBUG > 1

static void maybe_report_command_index (CMD *cmd, const char *s)
{
    if (cmd->ci > 0) {
	const char *word = gretl_command_word(cmd->ci);

	fprintf(stderr, "try_for_command_index: ci = %d (%s)",
		cmd->ci, word);
	if (strcmp(word, s)) {
	    fprintf(stderr, ", actual word = '%s'", s);
	}
	fputc('\n', stderr);
    }
}

#endif

static int is_plot_keyword (const char *s, int *err)
{
    if (!strcmp(s, "option") ||
	!strcmp(s, "options") ||
        !strcmp(s, "printf") ||
	!strcmp(s, "literal") ||
	*s == '@') {
	return 1;
    } else {
	gretl_errmsg_sprintf(_("'%s': invalid keyword in plot block"), s);
	*err = E_PARSE;
	return 0;
    }
}

/* If we have enough tokens parsed, try to determine the
   current command index.
*/

static int try_for_command_index (CMD *cmd, int i,
				  DATASET *dset,
				  int compmode,
				  int *err)
{
    cmd_token *toks = cmd->toks;
    const char *test = toks[i].s;

    cmd->ci = gretl_command_number(test);

    if (cmd->context && cmd->ci != END) {
	if (cmd->context == FOREIGN || cmd->context == MPI) {
	    /* Do not attempt to parse! Note: we get here only
	       when "compiling" a foreign block into a loop or
	       function.
	    */
	    cmd->ciflags = CI_EXPR;
	    cmd->ci = cmd->context;
	} else {
	    /* We're inside a "native" block of some kind.
	       In that case the line should be passed "as
	       is" to the cumulator for the block, with one
	       exception, namely the "equation" command
	       within a "system" block.
	    */
	    if (cmd->context == SYSTEM && !strcmp(test, "equations")) {
		cmd->ci = EQUATION;
		cmd->opt |= OPT_M;
	    } else if (cmd->context == SYSTEM && cmd->ci == EQUATION) {
		; /* OK */
	    } else if (cmd->context == PLOT && !is_plot_keyword(test, err)) {
		return 0;
	    } else {
		cmd->ci = cmd->context;
	    }
	}
    }

    if (cmd->ci > 0 && has_function_form(cmd->ci)) {
	/* disambiguate command versus function */
	if (peek_next_char(cmd, i) == '(') {
	    /* must be function form, not command proper */
	    cmd->ci = 0;
	    goto gentest;
	}
    } else if (cmd->ci == 0) {
	cmd->ci = try_for_command_alias(test, cmd);
    }

#if CDEBUG > 1
    maybe_report_command_index(cmd, test);
#endif

 gentest:

    if (cmd->ci <= 0 && cmd->ntoks < 5) {
	char cnext = peek_next_char(cmd, i);

	cmd->ci = test_for_genr(cmd, i, cnext, dset);
    }

    if (cmd->ci > 0) {
	mark_token_done(toks[i]);
	if (cmd->ci == cmd->context) {
	    cmd->ciflags = CI_EXPR;
	} else {
	    cmd->ciflags = command_get_flags(cmd->ci);
	    if (cmd->ci == EQUATION && (cmd->opt & OPT_M)) {
		/* the system "equations" keyword */
		cmd->ciflags ^= CI_LIST;
		cmd->ciflags |= CI_PARM1;
		cmd->ciflags |= CI_PARM2;
	    }
	    if (cmd->ci == STORE && (cmd->flags & CMD_PROG)) {
		cmd->ciflags ^= CI_LIST;
		cmd->ciflags ^= CI_DOALL;
		cmd->ciflags |= CI_EXTRA;
	    }
	    if (cmd->ci == GENR) {
		GretlType gtype = 0;

		if (!strcmp(test, "list")) {
		    if (peek_next_char(cmd, i) == '\0') {
			/* just "list" by itself */
			cmd->ci = VARLIST;
			cmd->ciflags = 0;
		    } else {
			/* probably "genr" but might be a special */
			gtype = cmd->gtype = GRETL_TYPE_LIST;
			cmd->ciflags |= CI_LCHK;
		    }
		} else if ((gtype = gretl_get_gen_type(test)) > 0) {
		    cmd->gtype = gtype;
		}
		if (gtype == 0 && get_list_by_name(test)) {
		    cmd->ciflags |= CI_LGEN;
		}
	    }
	    if (compmode && cmd->ci == END) {
		int endci = peek_end_param(cmd, i);

		if (compmode == FUNC && endci == FUNC) {
		    cmd->flags |= CMD_ENDFUN;
		} else if (endci == LOOP) {
		    gretl_errmsg_set(_("'end loop': did you mean 'endloop'"));
		    gretl_abort_compiling_loop();
		    *err = E_PARSE;
		}
	    }
	}
    }

    return cmd->ci;
}

/* Get a count of any tokens not marked as 'done'. */

static int count_remaining_toks (CMD *c)
{
    int i, n = 0;

    for (i=c->cstart+1; i<c->ntoks; i++) {
	if (!(c->toks[i].flag & TOK_DONE)) {
	    n++;
	}
    }

    return n;
}

static int process_auxlist_term (CMD *c, cmd_token *tok,
				 int **ilistptr)
{
    int *ilist = *ilistptr;
    int err = 0;

    if (ilist == NULL) {
	c->auxlist = get_auxlist(tok, &err);
	if (!err) {
	    /* add placeholder to @ilist */
	    gretl_list_append_term(ilistptr, list_max(c->auxlist));
	}
    } else if (c->ci == ARMA && ilist[0] < 3) {
	int *aux2 = get_auxlist(tok, &err);

	if (!err) {
	    int *tmp;

	    tmp = gretl_lists_join_with_separator(c->auxlist, aux2);
	    if (tmp == NULL) {
		c->err = E_ALLOC;
	    } else {
		gretl_list_append_term(ilistptr, list_max(aux2));
		free(c->auxlist);
		c->auxlist = tmp;
	    }
	}
	free(aux2);
    } else {
	/* we got too many auxlist terms */
	err = E_PARSE;
    }

    return err;
}

static int try_auxlist_term (CMD *cmd, cmd_token *tok, int scount)
{
    if (cmd->ci != ARMA && cmd->ci != DPANEL && cmd->ci != VAR) {
	/* only supported by these commands */
	return 0;
    }
    if (scount > 0) {
	/* we must be in first sublist */
	return 0;
    }
    if (tok->type == TOK_CBSTR || (tok->type == TOK_NAME &&
				   get_matrix_by_name(tok->s))) {
	/* not a regular "int list" term */
	return 1;
    } else {
	return 0;
    }
}

static int check_arma_ilist (const int *ilist,
			     const char *vtest)
{
    if (ilist != NULL && ilist[0] == 3) {
	/* got all 3 ARIMA fields: p, d, q */
	if (vtest[1] != 0) {
	    /* we got a vector term for 'd' */
	    return E_PARSE;
	}
    }

    return 0;
}

static int panel_gmm_special (CMD *cmd, const char *s)
{
    if (cmd->ci == DPANEL) {
	if (!strcmp(s, "GMM") || !strcmp(s, "GMMlevel")) {
	    return 1;
	}
    }

    return 0;
}

static int midas_term_special (CMD *cmd, const char *s)
{
    if (cmd->ci == MIDASREG) {
	if (!strcmp(s, "mds") || !strcmp(s, "mdsl")) {
	    return 1;
	}
    }

    return 0;
}

static void rejoin_list_toks (CMD *c, int k1, int *k2,
			      char *lstr, int j)
{
    cmd_token *tok;
    int i;

    if (j > 0) {
	strcat(lstr, " ");
    }

    for (i=k1; i<c->ntoks; i++) {
	tok = &c->toks[i];
	if (i > k1 && (!token_joined(tok) || token_done(tok))) {
	    break;
	}
	*k2 = i;
	if (tok->type == TOK_PRSTR) {
	    /* lag spec, for example */
	    strcat(lstr, "(");
	    strcat(lstr, tok->s);
	    strcat(lstr, ")");
	} else if (tok->type == TOK_BRSTR) {
	    /* element of named list, for example */
	    strcat(lstr, "[");
	    strcat(lstr, tok->s);
	    strcat(lstr, "]");
	} else if (i < c->ntoks - 1 && panel_gmm_special(c, tok->s)) {
	    cmd_token *next = &c->toks[i+1];
	    char *tmp;

	    tmp = gretl_strdup_printf("%s(%s)", tok->s, next->s);
	    c->param = gretl_str_expand(&c->param, tmp, " ");
	    free(tmp);
	    mark_list_token_done(c->toks[i]);
	    *k2 = ++i;
	} else if (i < c->ntoks - 1 && midas_term_special(c, tok->s)) {
	    cmd_token *next = &c->toks[i+1];
	    char *tmp;

	    tmp = gretl_strdup_printf("%s(%s)", tok->s, next->s);
	    c->param = gretl_str_expand(&c->param, tmp, " ");
	    free(tmp);
	    mark_list_token_done(c->toks[i]);
	    *k2 = ++i;
	} else {
	    strcat(lstr, tok->s);
	}
	mark_list_token_done(c->toks[i]);
    }
}

static int validate_list_token (cmd_token *tok)
{
    /* In a "genr" context we might accept a scalar or matrix as
       representing one or more integer series IDs in constructing a
       list, but this is not acceptable in the regular command
       context. So we'll screen out all "uservars" other than named
       lists (this does not include series, which are of course OK).
    */

    if (tok->type == TOK_NAME) {
	GretlType t = user_var_get_type_by_name(tok->s);

	return t == GRETL_TYPE_NONE || t == GRETL_TYPE_LIST;
    } else {
	return 1;
    }
}

/* In case we got something that looks like a list but
   turned out not to be valid as such, remove the TOK_DONE
   flag from the tokens that composed the putative list,
   which will be marked with TOK_LSTR status.
*/

static void rescind_tok_done_status (CMD *c)
{
    cmd_token *tok;
    int i;

    for (i=c->cstart+1; i<c->ntoks; i++) {
	tok = &c->toks[i];
	if (tok->flag & TOK_LSTR) {
	    tok->flag ^= TOK_DONE;
	    tok->flag ^= TOK_LSTR;
	}
    }
}

static int process_command_list (CMD *c, DATASET *dset)
{
    guint8 TOK_PROV = (TOK_DONE | TOK_LSTR);
    char vectest[3] = {0}; /* for arima */
    char lstr[MAXLINE];
    cmd_token *tok;
    int *ilist = NULL;
    int *vlist = NULL;
    int want_ints = 0;
    int scount = 0;
    int i, j, k, ns;

    if (c->ciflags & CI_L1INT) {
	want_ints = 1;
    }

    ns = cmd_get_sepcount(c);

#if CDEBUG
    fprintf(stderr, "process command list, ntoks = %d, sepcount %d\n",
	    c->ntoks, ns);
#endif

    *lstr = '\0';
    j = 0;

    for (i=c->cstart+1; i<c->ntoks && !c->err; i++) {
	tok = &c->toks[i];
	if (tok->type == TOK_SEMIC) {
	    scount++;
	    if (c->ci == LAGS && scount == 1) {
		/* the separator that follows the optional number
		   of lags to create */
		tok->flag |= TOK_DONE;
	    } else if (c->ci == ARMA && ns == 2) {
		if (scount == 1) {
		    gretl_list_append_term(&ilist, LISTSEP);
		    tok->flag |= TOK_DONE;
		} else if (scount == 2) {
		    want_ints = 0;
		    tok->flag |= TOK_DONE;
		}
	    } else if (want_ints) {
		want_ints = 0;
		tok->flag |= TOK_DONE;
	    } else if (c->ci == MPOLS) {
		want_ints = 1;
		tok->flag |= TOK_DONE;
	    }
	}
	if (!token_done(tok)) {
	    if (bracketed_type(tok->type) && c->ci != PRINT) {
		gretl_errmsg_sprintf(_("Parse error at unexpected token '%s'"),
				     tok->type == TOK_PRSTR ? "(" : "[");
		c->err = E_PARSE;
	    } else if (want_ints) {
		if (try_auxlist_term(c, tok, scount)) {
		    /* a vector-style entry */
		    c->err = process_auxlist_term(c, tok, &ilist);
		    if (!c->err) {
			vectest[ilist[0] - 1] = 1;
			tok->flag |= TOK_PROV;
		    }
		} else {
		    k = gretl_int_from_string(tok->s, &c->err);
		    if (!c->err) {
			gretl_list_append_term(&ilist, k);
			tok->flag |= TOK_PROV;
		    }
		}
		if (c->err) {
		    gretl_errmsg_sprintf(_("Parse error at unexpected token '%s'"),
					 tok->s);
		}
	    } else if (next_joined_token(c, i) != NULL) {
		rejoin_list_toks(c, i, &i, lstr, j++);
	    } else if (validate_list_token(tok)) {
		if (j > 0) {
		    strcat(lstr, " ");
		}
		strcat(lstr, tok->s);
		tok->flag |= TOK_PROV;
		j++;
	    }
	}
	if (j == 1 && (c->ciflags & CI_LLEN1)) {
	    break;
	} else if (j == 2 && (c->ciflags & CI_LLEN2)) {
	    break;
	}
    }

    if (!c->err) {
	c->err = check_list_sepcount(c->ci, ns);
    }

    if (!c->err && c->ci == ARMA) {
	c->err = check_arma_ilist(ilist, vectest);
    }

    if (!c->err && *lstr != '\0') {
	tailstrip(lstr);
	if (c->ci == DPANEL || c->ci == MIDASREG) {
	    /* We may have a ';' separator that's not followed
	       by any regular second list, just special terms; so
	       don't error out on a trailing ';' in defining a
	       list.
	    */
	    if (lstr[strlen(lstr)-1] == ';') {
		lstr[strlen(lstr)-1] = '\0';
	    }
	}
    }

#if CDEBUG > 1
    fprintf(stderr, "process_command_list: lstr='%s' (err=%d)\n", lstr, c->err);
#endif

    if ((c->ci == DELEET || c->ci == PRINT) && *lstr == '\0') {
	/* we didn't get a "list string": maybe the terms are
	   names of non-series variables
	*/
	c->ciflags &= ~CI_LIST;
	c->ciflags &= ~CI_DOALL;
	c->ciflags |= CI_ADHOC;
	c->err = 0;
    }

    if (!c->err && dset != NULL && *lstr != '\0') {
	vlist = generate_list(lstr, dset, c->ci, &c->err);
	if (c->err && (c->ci == DELEET || c->ci == PRINT)) {
	    /* we got something that looked like a list string,
	       but list generation failed: again, maybe the
	       the terms are names of non-series variables
	    */
	    rescind_tok_done_status(c);
	    c->ciflags &= ~CI_LIST;
	    c->ciflags &= ~CI_DOALL;
	    c->ciflags |= CI_ADHOC;
	    c->err = 0;
	    goto finish;
	}
    }

    if (!c->err) {
	if (c->ci == MPOLS && vlist != NULL && ilist != NULL) {
	    /* legacy mpols special */
	    c->list = gretl_lists_join_with_separator(vlist, ilist);
	    if (c->list == NULL) {
		c->err = E_ALLOC;
	    }
	} else if (vlist != NULL) {
	    if (ilist != NULL) {
		c->list = gretl_lists_join_with_separator(ilist, vlist);
		if (c->list == NULL) {
		    c->err = E_ALLOC;
		}
	    } else {
		c->list = vlist;
		vlist = NULL;
	    }
	} else if (ilist != NULL) {
	    /* a "pure" ints list */
	    c->list = ilist;
	    ilist = NULL;
	}
    }

 finish:

    free(ilist);
    free(vlist);

#if CDEBUG > 1
    fprintf(stderr, "process_command_list: returning err = %d\n", c->err);
#endif

    return c->err;
}

static int handle_adhoc_string (CMD *c)
{
    PRN *prn;
    cmd_token *tok;
    int i, j = 0;
    int err = 0;

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);
    if (err) {
	return err;
    }

    if (c->ci == FUNC) {
	pputs(prn, "function ");
    }

    for (i=c->cstart+1; i<c->ntoks; i++) {
	tok = &c->toks[i];
	/* 2019-02-27: was !option_type(tok->type) below */
	if (!token_done(tok)) {
	    if (j > 0 && !token_joined(tok)) {
		pputc(prn, ' ');
	    }
	    if (tok->flag & TOK_QUOTED) {
		pprintf(prn, "\"%s\"", tok->s);
	    } else if (tok->type == TOK_PRSTR) {
		pprintf(prn, "(%s)", tok->s);
	    } else if (tok->type == TOK_CBSTR) {
		pprintf(prn, "{%s}", tok->s);
	    } else if (tok->type == TOK_BRSTR) {
		pprintf(prn, "[%s]", tok->s);
	    } else {
		pputs(prn, tok->s);
	    }
	    j++;
	}
    }

    if (c->vstart != NULL) {
	const char *s = c->vstart;

	s += strspn(s, " ");

	if (*s) {
	    if (j > 0) {
		pputc(prn, ' ');
	    }
	    pputs(prn, s);
	} else {
	    c->vstart = NULL;
	}
    }

    /* In general, any ad hoc portion of a command should
       be stuffed into its @param member. However, with
       PRINT, @param is reserved for a string literal so
       we use @parm2 instead.
    */

    if (c->ci == PRINT) {
	c->parm2 = gretl_print_steal_buffer(prn);
    } else {
	c->param = gretl_print_steal_buffer(prn);
    }

    gretl_print_destroy(prn);

    return err;
}

static int n_regular_tokens (CMD *c)
{
    int i, n = c->ntoks;

    if (c->opt) {
	for (i=0; i<c->ntoks; i++) {
	    if (c->toks[i].type == TOK_OPT ||
		c->toks[i].type == TOK_SOPT ||
		c->toks[i].type == TOK_OPTDASH ||
		c->toks[i].type == TOK_OPTEQ ||
		c->toks[i].type == TOK_OPTVAL) {
		n--;
	    }
	}
    }

    return n;
}

static int handle_command_extra (CMD *c)
{
    cmd_token *tok;
    int i;

    if (c->ci == GNUPLOT || c->ci == BXPLOT ||
	c->ci == HFPLOT || c->ci == PANPLOT) {
	/* if present, 'extra' goes into param */
	for (i=c->cstart+1; i<c->ntoks; i++) {
	    tok = &c->toks[i];
	    if (!token_done(tok) && tok->type == TOK_CBSTR) {
		/* catch stuff in braces */
		tok->flag |= TOK_DONE;
		c->param = tok->s;
	    }
	}
    } else if (c->ci == MODPRINT) {
	/* if present, 'extra' gets pushed as an option */
	int regtoks = n_regular_tokens(c);

	if (regtoks - c->cstart - 1 == 3) {
	    /* got three arguments */
	    char *extra;

	    tok = &c->toks[regtoks - 1];
	    if (!token_done(tok) && tok->type == TOK_NAME) {
		tok->flag |= TOK_DONE;
		extra = gretl_strdup(tok->s);
		if (extra == NULL) {
		    c->err = E_ALLOC;
		} else {
		    c->opt |= OPT_A;
		    c->err = push_option_param(c->ci, OPT_A, extra);
		}
	    } else if (!token_done(tok)) {
		gretl_errmsg_set(_("modprint: third argument must be a named matrix"));
		c->err = E_TYPES;
	    }
	}
    } else if (c->ci == STORE) {
	/* progressive loop: 'extra' goes into parm2 */
	for (i=c->cstart+1; i<c->ntoks; i++) {
	    tok = &c->toks[i];
	    if (!token_done(tok) && tok->type == TOK_NAME) {
		tok->flag |= TOK_DONE;
		c->parm2 = gretl_str_expand(&c->parm2, tok->s, " ");
	    }
	}
    } else if (c->ci == DATAMOD) {
	/* "datamod compact": 'extra' into parm2 */
	for (i=c->cstart+1; i<c->ntoks; i++) {
	    tok = &c->toks[i];
	    if (!token_done(tok) &&
		(tok->type == TOK_NAME || tok->type == TOK_INT)) {
		tok->flag |= TOK_DONE;
		c->parm2 = gretl_str_expand(&c->parm2, tok->s, " ");
	    }
	}
    } else if (c->ci == JOIN) {
	/* join: allow for multiple import names */
	for (i=c->cstart+1; i<c->ntoks; i++) {
	    tok = &c->toks[i];
	    if (!token_done(tok) && tok->type == TOK_NAME) {
		tok->flag |= TOK_DONE;
		c->parm2 = gretl_str_expand(&c->parm2, tok->s, " ");
	    } else if (!token_done(tok) && wildsym(tok)) {
		tok->flag |= TOK_DONE;
		c->parm2 = gretl_str_expand(&c->parm2, tok->s, "");
	    }
	}
    }

    return c->err;
}

/* @vstart is a const pointer into the incoming command
   line, holding a "genr"-type expression, a string to
   be passed to the shell, or a varargs expression.
*/

static int set_command_vstart (CMD *cmd, ExecState *state,
			       PRN *prn)
{
    cmd_token *tok;
    const char *s = NULL;

    if (cmd->ciflags & CI_EXPR) {
	tok = &cmd->toks[cmd->cstart];
	s = tok->lp;
	if (!cmd->context && expr_keep_cmdword(cmd->ci)) {
	    ; /* leave it alone */
	} else if (cmd->ci == EVAL) {
	    cmd->gtype = GRETL_TYPE_NONE;
	    s += 4;
	} else if (cmd->ci == GENR && cmd->gtype != GRETL_TYPE_ANY) {
	    /* skip initial command word */
	    s += strlen(tok->s);
	} else if (cmd->ci != GENR && cmd->ci != cmd->context) {
	    /* skip initial command word */
	    s += strlen(tok->s);
	}
    } else if (cmd->ciflags & CI_VARGS) {
	/* vstart should point beyond the last token */
	tok = &cmd->toks[cmd->ntoks-1];
	s = tok->lp + real_toklen(tok);
    }

    if (s != NULL) {
	s += strspn(s, " \t");
	if (*s == '\0') {
	    s = NULL;
	}
    }

    if (cmd->gtype == GRETL_TYPE_LIST || (cmd->ciflags & CI_LGEN)) {
	/* we won't accept ';' as list separator outside of an
	   appropriate command context. This breaks the old version
	   of the system "equations" mechanism but I'm afraid that's
	   just too bad.
	*/
	if (strchr(s, ';')) {
	    if (prn != NULL) {
		pputc(prn, '\n');
		pputs(prn, "If you are trying to assemble a compound list for use in a "
		      "\"system\" block,\nplease see section 31.2 of the Gretl User's "
		      "Guide for the current method.\n");
		pputc(prn, '\n');
	    }
	    gretl_errmsg_sprintf(_("The symbol '%c' is not valid in this context\n"),
				 ';');
	    return E_INVARG;
	}
    }

    cmd_set_vstart(cmd, s);

    return 0;
}

/* For a command that ends with varargs, do we have the required
   leading non-vararg parameter(s)? (This means either one or two
   parameters, followed by a comma.) If so, we can stop parsing
   and designate the remainder of the command line as the
   varargs portion.
*/

static int got_param_tokens (CMD *cmd)
{
    int mintoks = (cmd->ciflags & CI_PARM2)? 3 : 2;
    int i, n = cmd->ntoks - cmd->cstart - 1;

    for (i=cmd->cstart+1; i<cmd->ntoks; i++) {
	if (cmd->toks[i].flag & TOK_JOINED) {
	    /* avoid over-counting */
	    n--;
	}
    }

    if (cmd->ci == PRINTF) {
	i = cmd->cstart + 2;
	if (i < cmd->ntoks && cmd->toks[i].type == TOK_COMMA) {
	    if (cmd->ci == PRINTF) {
		mark_token_done(cmd->toks[i]);
	    } else {
		/* sprintf: redundant comma after varname */
		mark_token_ignored(cmd->toks[i]);
		mintoks++;
	    }
	}
    }

    return n == mintoks;
}

static int check_end_command (CMD *cmd)
{
    int endci = gretl_command_number(cmd->param);

    if (endci == OUTFILE || endci == GPBUILD) {
	/* special cases where no "context" is required */
	return 0;
    }

    if (endci == 0 || endci != cmd->context) {
	gretl_errmsg_sprintf(_("end: invalid parameter '%s'"), cmd->param);
	cmd->err = E_DATA;
    }

    /* on "end", scrub the context */
    cmd->context = 0;

    return cmd->err;
}

/* For a command that (usually) requires a list, check
   that we got one, and if so, check that it doesn't
   contain duplicates.
*/

static int check_for_list (CMD *cmd)
{
    if (cmd->list == NULL) {
	if (cmd->ciflags & CI_DOALL) {
	    ; /* list defaults to all series, OK */
	} else if (cmd->ci == OMIT && (cmd->opt & OPT_A)) {
	    ; /* the auto-omit option, OK */
	} else if ((cmd->ci == FREQ || cmd->ci == BDS) && (cmd->opt & OPT_X)) {
	    ; /* using a matrix: may be OK */
	} else {
	    cmd->err = E_ARGS;
	}
    } else {
	/* check for duplicated variables */
	int dupv = gretl_list_duplicates(cmd->list, cmd->ci);

	if (dupv >= 0) {
	    printlist(cmd->list, "command with duplicate(s)");
	    cmd->err = E_DATA;
	    gretl_errmsg_sprintf(_("variable %d duplicated in the "
				   "command list."), dupv);
	}
    }

    return cmd->err;
}

/* @cmd has the CI_LCHK (ambiguity) flag set: see if
   we're able to disambiguate by this point
*/

static int scrub_list_check (CMD *cmd)
{
    int maxtoks = cmd->toks[0].type == TOK_CATCH ? 4 : 3;
    int ret = 0;

    if (cmd->ntoks == maxtoks) {
	const char *s = cmd->toks[maxtoks-1].s;
	int ci = gretl_command_number(s);

	if (ci == DELEET || ci == PRINT) {
	    cmd->toks[maxtoks-1].flag |= TOK_DONE;
	    cmd->ci = ci;
	    cmd->opt = OPT_L;
	    cmd->ciflags = CI_PARM1;
	    ret = 1;
	} else {
	    cmd->ciflags ^= CI_LCHK;
	    ret = 1;
	}
    }

    return ret;
}

static int unexpected_symbol_error (char c)
{
    if (c == '\'') {
	gretl_errmsg_sprintf(_("Unexpected symbol %c"), c);
    } else {
	gretl_errmsg_sprintf(_("Unexpected symbol '%c'"), c);
    }
    return E_PARSE;
}

static int utf8_fail (char *s)
{
    if (!g_utf8_validate(s, -1, NULL)) {
	gretl_errmsg_set(_("Command line is not valid UTF-8"));
	return E_DATA;
    } else {
	return 0;
    }
}

static int print_do_eval (CMD *cmd, int imin, char *s)
{
    char c = peek_next_char(cmd, imin);

    if (c == '$') {
	/* display value of accessor */
	cmd->ci = EVAL;
	cmd->gtype = GRETL_TYPE_NONE;
	cmd_set_vstart(cmd, strchr(s, '$'));
	return 1;
    } else if (c != '"') {
	cmd->ciflags |= CI_LIST;
    }

    return 0;
}

/* @c is A-Z or a-z */

static int is_ascii_alpha (int c)
{
    return (c >= 0x41 && c <= 0x5A) ||
	(c >= 0x61 && c <= 0x7A);
}

#define MAY_START_NUMBER(c) (c == '.' || c == '-' || c == '+')

/* tokenize_line: parse @line into a set of tokens on a
   lexical basis. In some cases constitution of command
   arguments will require compositing tokens. We get a
   little semantic help from the CI_FNAME flag: if this
   is present for a given command, that tells us to
   consider a filename containing directory separators
   as a unitary token.
*/

static int tokenize_line (ExecState *state, DATASET *dset,
			  int compmode)
{
    char tok[FN_NAMELEN];
    char *s = state->line;
    CMD *cmd = state->cmd;
    char *vtok;
    int n, m, pos = 0;
    int wild_ok = 0;
    int at_ok = compmode;
    int want_fname = 0;
    int err = 0;

#if CDEBUG || TDEBUG
    fprintf(stderr, "*** %s: line = '%s'\n",
	    compmode ? "get_command_index" : "parse_command_line", s);
#endif

    if (utf8_fail(s)) {
	return E_DATA;
    }

    gretl_push_c_numeric_locale();

    if (!gretl_in_batch_mode() && *s == '=') {
	/* treat as a bare expression to be evaluated */
	cmd->ci = EVAL;
	cmd->gtype = GRETL_TYPE_NONE;
	cmd_set_vstart(cmd, s + 1);
	goto skipit;
    }

    while (!err && *s) {
	int skipped = 0;

	*tok = '\0';

	if (*s == '-') {
	    want_fname = 0;
	}

	if (*s == '#') {
	    break;
	} else if (want_fname && *s != '"' && !isspace(*s)) {
	    n = strcspn(s, " \t");
	    if (n < FN_NAMELEN) {
		strncat(tok, s, n);
		err = push_string_token(cmd, tok, s, pos);
	    } else {
		vtok = gretl_strndup(s, n);
		if (vtok == NULL) {
		    err = E_ALLOC;
		} else {
		    err = push_string_token(cmd, vtok, s, pos);
		    free(vtok);
		}
	    }
	} else if (wild_ok && (is_ascii_alpha((int) *s) || *s == '*')) {
	    n = 1 + wild_spn(s+1);
	    m = (n < FN_NAMELEN)? n : FN_NAMELEN - 1;
	    strncat(tok, s, m);
	    err = push_string_token(cmd, tok, s, pos);
	} else if (is_ascii_alpha((int) *s) || *s == '$' || (at_ok && *s == '@')) {
	    /* regular or accessor identifier */
	    if (*s == '@' && !compmode) {
		fprintf(stderr, "tokenize: found '@':\n '%s'\n", state->line);
	    }
	    n = 1 + namechar_spn(s+1);
	    m = (n < FN_NAMELEN)? n : FN_NAMELEN - 1;
	    strncat(tok, s, m);
	    err = push_string_token(cmd, tok, s, pos);
	} else if (is_greek_letter(s)) {
	    n = 2;
	    strncat(tok, s, n);
	    err = push_string_token(cmd, tok, s, pos);
	} else if (ldelim(*s)) {
	    /* left-hand delimiter that needs to be paired */
	    n = closing_delimiter_pos(s);
	    if (n < 0) {
		gretl_errmsg_sprintf(_("Unmatched '%c'\n"), *s);
		err = E_PARSE;
	    } else if (n < FN_NAMELEN) {
		strncat(tok, s+1, n);
		err = push_delimited_token(cmd, tok, s, pos);
	    } else {
		vtok = gretl_strndup(s+1, n);
		if (vtok == NULL) {
		    err = E_ALLOC;
		} else {
		    err = push_delimited_token(cmd, vtok, s, pos);
		    free(vtok);
		}
	    }
	    n += 2;
	} else if (*s == '"') {
	    int esc = 0;

	    n = closing_quote_pos(s, cmd->ci, &esc);
	    if (n < 0) {
		gretl_errmsg_sprintf(_("Unmatched '%c'\n"), '"');
		err = E_PARSE;
	    } else {
		err = push_quoted_token(cmd, s, n, pos, esc);
	    }
	    n += 2;
	} else if ((n = symbol_spn(s)) > 0) {
	    if (n == 1 && MAY_START_NUMBER(*s) && isdigit(*(s+1))) {
		n = numeric_spn(s, 0);
		if (n == 0) {
		    err = unexpected_symbol_error(*s);
		} else {
		    m = (n < FN_NAMELEN)? n : FN_NAMELEN - 1;
		    strncat(tok, s, m);
		    err = push_numeric_token(cmd, tok, s, pos);
		}
	    } else {
		/* operator / symbol */
		m = (n < FN_NAMELEN)? n : FN_NAMELEN - 1;
		strncat(tok, s, m);
		err = push_symbol_token(cmd, tok, s, pos);
	    }
	} else if (isdigit(*s)) {
	    /* numeric string */
	    n = numeric_spn(s, 1);
	    m = (n < FN_NAMELEN)? n : FN_NAMELEN - 1;
	    strncat(tok, s, m);
	    err = push_numeric_token(cmd, tok, s, pos);
	} else if (isspace(*s) || *s == (char) 0xA0) {
	    /* handle stupid "non-breaking space" here too */
	    n = 1;
	    skipped = 1;
	} else if (*s == '@' && (compmode || gretl_if_state_false())) {
	    /* string substitution not yet done */
	    n = 1;
	    skipped = 1;
	} else if (*s == '_') {
	    /* unassigned call to a hidden function? */
	    n = 1 + namechar_spn(s+1);
	    m = (n < FN_NAMELEN)? n : FN_NAMELEN - 1;
	    strncat(tok, s, m);
	    err = push_string_token(cmd, tok, s, pos);
	} else {
	    err = unexpected_symbol_error(*s);
	}

	if (err) {
	    break;
	}

	if (!skipped && want_fname) {
	    want_fname = 0;
	}

	if (!skipped && cmd->ci == 0 && cmd->ntoks > 0) {
	    /* use current info to determine command index? */
	    int imin = min_token_index(cmd, compmode);

	    if (cmd->ntoks > imin) {
		try_for_command_index(cmd, imin, dset, compmode, &err);
#if TDEBUG
		if (cmd->ci > 0) {
		    fprintf(stderr, "ntoks=%d, imin=%d, ci=%d (%s)\n",
			    cmd->ntoks, imin, cmd->ci, gretl_command_word(cmd->ci));
		} else {
		    fprintf(stderr, "ntoks=%d, imin=%d, ci not yet known\n",
			    cmd->ntoks, imin);
		}
#endif
		if (cmd->ci == PRINT && print_do_eval(cmd, imin, s)) {
		    goto skipit;
		}
		if (cmd->ciflags & CI_FNAME) {
		    want_fname = 1;
		}
		if (cmd->ci == LOOP && compmode == LOOP) {
		    /* we need to pick up the specifics of the loop */
		    compmode = 0;
		}
	    }
	}

	/* when we're just looking for a command index (compmode),
	   we may be able to get out early
	*/
	if (compmode && simple_flow_control(cmd)) {
	    if (string_is_blank(s)) {
		break;
	    }
	} else if (compmode && (cmd->ci > 0 || cmd->ntoks == 3)) {
	    /* either we've got the command index or it seems
	       we're not going to get it */
	    break;
	}

	if (cmd->ciflags & CI_LCHK) {
	    /* handle ambiguity of "list ..." */
	    wild_ok = !scrub_list_check(cmd);
	}

	if (cmd->ci == DELEET || (cmd->ciflags & CI_LIST)) {
	    /* flag acceptance of wildcard expressions */
	    wild_ok = 1;
	} else if ((cmd->ciflags & CI_EXPR) && !(cmd->ciflags & CI_LCHK)) {
	    /* the remainder of line will be parsed elsewhere */
	    break;
	} else if ((cmd->ciflags & CI_ADHOC) && (cmd->ciflags & CI_NOOPT)) {
	    /* ditto */
	    cmd_set_vstart(cmd, s + n);
	    break;
	}

	if ((cmd->ciflags & CI_VARGS) && got_param_tokens(cmd)) {
	    /* remaining args will be parsed elsewhere */
	    break;
	}

	s += n;
	pos += n;
    }

 skipit:

    gretl_pop_c_numeric_locale();

#if CDEBUG
    if (err) {
	fprintf(stderr, "tokenize_line: err = %d\n", err);
    }
#endif

    return err;
}

/* for use with spreadsheet option params */

static int small_positive_int (const char *s)
{
    if (integer_string(s)) {
	int k = atoi(s);

	if (k > 0 && k <= 10) {
	    return 1;
	}
    }

    return 0;
}

/* The following is kind of a curiosity -- translation from
   option parameters to a special list, for handling
   spreadsheet-specific options for "open" or "append".
   For for the moment I'm going to leave it close to what
   was in the "old" (pre-tokenize) interact.c. It can be
   revisited later. AC, 2014-08-30
*/

static int post_process_spreadsheet_options (CMD *cmd)
{
    int err = 0;

    if (cmd->opt & OPT_O) {
	/* odbc: spreadsheet-specific options not acceptable */
	err = incompatible_options(cmd->opt, OPT_O | OPT_C |
				   OPT_R | OPT_S);
    } else if (cmd->opt & OPT_W) {
	/* web database: ditto */
	err = incompatible_options(cmd->opt, OPT_W | OPT_C |
				   OPT_R | OPT_S);
    }

    if (!err) {
	err = incompatible_options(cmd->opt, OPT_O | OPT_W);
    }

    if (!err && (cmd->opt & (OPT_R | OPT_C | OPT_S))) {
	/* row offset, column offset, sheet name/number */
	const char *s = NULL;
	int r0 = 0, c0 = 0;

	if (cmd->opt & OPT_R) {
	    /* --rowoffset */
	    r0 = get_optval_int(cmd->ci, OPT_R, &err);
	}

	if (!err && (cmd->opt & OPT_C)) {
	    /* --coloffset */
	    c0 = get_optval_int(cmd->ci, OPT_C, &err);
	}

	if (!err && (cmd->opt & OPT_S)) {
	    /* --sheet */
	    s = get_optval_string(cmd->ci, OPT_S);
	    if (s == NULL) {
		err = E_DATA;
	    }
	}

	if (!err) {
	    int slist[4] = {3, 0, c0, r0};

	    free(cmd->list);
	    cmd->list = gretl_list_copy(slist);
	    if (cmd->list == NULL) {
		err = E_ALLOC;
	    } else {
		/* note: dodgy heuristic here? */
		if (small_positive_int(s)) {
		    /* take the --sheet spec as giving a sheet
		       number (1-based) */
		    cmd->list[1] = atoi(s);
		} else if (s != NULL) {
		    /* take it as giving a sheet name */
		    free(cmd->parm2);
		    cmd->parm2 = gretl_strdup(s);
		    if (cmd->parm2 == NULL) {
			err = E_ALLOC;
		    }
		}
	    }
	}
    }

    return err;
}

static int post_process_rename_param (CMD *cmd,
				      DATASET *dset)
{
    int err = 0;

    if (integer_string(cmd->param)) {
	cmd->auxint = atoi(cmd->param);
	if (cmd->auxint < 1 || cmd->auxint >= dset->v) {
	    err = E_DATA;
	}
    } else {
	cmd->auxint = current_series_index(dset, cmd->param);
	if (cmd->auxint < 0) {
	    err = E_UNKVAR;
	}
    }

    return err;
}

/* check the commands that have the CI_INFL flag:
   the precise line-up of required arguments may
   depend on the option(s) specified
*/

static void handle_option_inflections (CMD *cmd)
{
    if (cmd->ci == BXPLOT) {
	if (cmd->opt & OPT_Z) {
	    /* factorized: two variables wanted */
	    cmd->ciflags |= CI_LLEN2;
	} else if (cmd->opt & OPT_X) {
	    /* with --matrix, default to all columns */
	    cmd->ciflags |= CI_DOALL;
	}
    } else if (cmd->ci == SMPL) {
	if (cmd->opt & (OPT_M | OPT_A | OPT_C)) {
	    /* no-missing, no-all-missing or contiguous */
	    cmd->ciflags = CI_LIST | CI_DOALL;
	} else if (cmd->opt & OPT_R) {
	    /* restrict */
	    cmd->ciflags = CI_ADHOC;
	} else if (cmd->opt & OPT_F) {
	    /* full: no args */
	    cmd->ciflags = 0;
	} else if (cmd->opt & OPT_O) {
	    /* using dummy variable */
	    cmd->ciflags &= ~CI_PARM2;
	} else if (cmd->opt & OPT_N) {
	    /* random sample */
	    cmd->ciflags = CI_PARM1;
	}
    } else if (cmd->ci == SET) {
	if (cmd->opt & (OPT_F | OPT_T)) {
	    /* from file, to file */
	    cmd->ciflags = 0;
	}
    } else if (cmd->ci == GNUPLOT) {
	if (cmd->opt & (OPT_I | OPT_i)) {
	    /* we got --input or --inbuf, so no args wanted,
	       but still accept "{...}" stuff?
	    */
	    cmd->ciflags = CI_EXTRA;
	} else if (cmd->opt & OPT_X) {
	    /* with --matrix, default to all columns */
	    cmd->ciflags |= CI_DOALL;
	}
    } else if (cmd->ci == OPEN) {
	if (cmd->opt & OPT_O) {
	    /* --odbc */
	    cmd->ciflags = CI_ADHOC;
	}
    } else if (cmd->ci == DELEET) {
	if (cmd->opt == OPT_T) {
	    /* --type=... */
	    cmd->ciflags = 0;
	} else if (cmd->opt == OPT_D) {
	    /* --db */
	    cmd->ciflags = CI_ADHOC;
	}
    } else if (cmd->ci == PRINT) {
	if (cmd->opt == OPT_L) {
	    /* --list */
	    cmd->ciflags = CI_PARM1;
	}
    } else if (cmd->ci == MODELTAB) {
	if (cmd->opt == OPT_O) {
	    /* --output: no arg needed */
	    cmd->ciflags &= ~CI_PARM1;
	}
    } else if (cmd->ci == XTAB) {
	if (cmd->opt & OPT_X) {
	    /* --matrix: no list wanted */
	    cmd->ciflags &= ~CI_LIST;
	}
    } else if (cmd->ci == SETINFO) {
	if (cmd->opt & (OPT_M | OPT_C | OPT_D)) {
	    /* midas, continuous or discrete */
	    if (!(cmd->opt & (OPT_G | OPT_I))) {
		/* but not graph-name or description */
		cmd->ciflags &= ~CI_LLEN1;
	    }
	}
    }
}

/* Check whether a "print" command includes a function call
   taking a pointer argument: we don't allow this.
*/

static int command_includes_fncall (CMD *c)
{
    const char *s;
    int i, ptr;
    int err = 0;

    for (i=c->cstart+1; i<c->ntoks && !err; i++) {
	if (i > 1 && c->toks[i].type == TOK_PRSTR) {
	    ptr = strchr(c->toks[i].s, '&') != NULL;
	    s = c->toks[i-1].s;
	    if (ptr && (function_lookup(s) || is_user_function(s))) {
		err = 1;
	    }
	}
    }

    if (err) {
	gretl_errmsg_set(_("\"print\": pointer arguments not allowed: "
                           "please use \"eval\" or \"printf\""));
    }

    return err;
}

static int assemble_command (CMD *cmd, DATASET *dset,
			     ExecState *s, char *line)
{
    /* defer handling option(s) till param is known? */
    int options_later = cmd->ci == SETOPT;
    PRN *prn = NULL;

    if (cmd->ntoks == 0) {
	return cmd->err;
    }

    if (s != NULL) {
	prn = s->prn;
    }

#if CDEBUG > 1
    fprintf(stderr, "doing assemble_command...\n");
#endif

    if (!never_takes_options(cmd) && !options_later) {
	cmd->err = check_command_options(cmd);
    }

    if (!cmd->err) {
	handle_command_preamble(cmd);
    }

    if (!cmd->err && cmd->ci == PRINT) {
	cmd->err = command_includes_fncall(cmd);
    }

    if (cmd->err) {
	goto bailout;
    }

    if (matrix_data_option(cmd->ci, cmd->opt)) {
	/* using matrix argument, plain ints ok in list */
	cmd->ciflags |= CI_L1INT;
    }

    if (cmd->ci == PRINT && !(cmd->opt & OPT_L)) {
	if (first_arg_quoted(cmd)) {
	    /* printing a string literal */
	    cmd->ciflags = CI_PARM1;
	} else if (dset != NULL && dset->v > 0) {
	    /* assume for now that we're printing series */
	    cmd->ciflags |= (CI_LIST | CI_DOALL);
	}
    } else if (option_inflected(cmd)) {
	handle_option_inflections(cmd);
    }

    /* legacy stuff */

    if (cmd->ci == TABPRINT || cmd->ci == TABPRINT) {
	legacy_get_filename(cmd);
    } else if (cmd->ci == SETINFO) {
	get_quoted_dash_fields(cmd, "dn");
    }

    if (cmd->err) {
	goto bailout;
    }

    /* main command assembly begins */

    if (cmd->ciflags & CI_PARM1) {
	get_param(cmd, dset);
    } else if (cmd->ciflags & CI_ORD1) {
	get_command_order(cmd);
    }

    if (!cmd->err && cmd->ci == VECM) {
	get_vecm_rank(cmd);
    }

    if (!cmd->err && (cmd->ciflags & CI_EXTRA)) {
	handle_command_extra(cmd);
    }

    if (!cmd->err && options_later) {
	cmd->err = check_command_options(cmd);
    }

    if (!cmd->err && (cmd->ciflags & CI_PARM2)) {
	get_parm2(cmd, options_later);
    }

    if (!cmd->err && (cmd->ciflags & CI_LIST)) {
	if (count_remaining_toks(cmd) > 0) {
	    process_command_list(cmd, dset);
	}
    }

    if (!cmd->err && (cmd->ciflags & CI_ORD2)) {
	get_optional_order(cmd);
    }

    if (!cmd->err && (cmd->ciflags & CI_LIST)) {
	check_for_list(cmd);
    }

    if (!cmd->err) {
	if (cmd->ciflags & CI_ADHOC) {
	    handle_adhoc_string(cmd);
	} else if (cmd->ciflags & (CI_EXPR | CI_VARGS)) {
	    cmd->err = set_command_vstart(cmd, s, prn);
	}
    }

    if (!cmd->err && cmd->ci == END) {
	check_end_command(cmd);
    }

    if (!cmd->err) {
	check_for_stray_tokens(cmd);
    }

 bailout:

#if CDEBUG
    print_tokens(cmd);
#endif

    if (!cmd->err) {
	if (cmd->opt != OPT_NONE &&
	    (cmd->ci == OPEN || cmd->ci == APPEND)) {
	    cmd->err = post_process_spreadsheet_options(cmd);
	} else if (cmd->ci == RENAME) {
	    cmd->err = post_process_rename_param(cmd, dset);
	}
    }

    return cmd->err;
}

static void maybe_init_shadow (void)
{
    static int shadow_initted;

    if (!shadow_initted) {
	check_for_shadowed_commands();
	shadow_initted = 1;
    }
}

static char *get_or_set_errline (const char *s, int set)
{
    static char *errline;

    if (set) {
	free(errline);
	errline = gretl_strdup(s);
    }

    return errline;
}

const char *get_parser_errline (void)
{
    return get_or_set_errline(NULL, 0);
}

/* When the current "if-state" is FALSE, and we're scanning
   the current command/statement, it should be sufficient to
   determine if we have a command which potentially modifies
   the if-state.
*/

static int get_flow_control_ci (ExecState *state)
{
    char word[6], *s = state->line;
    int ci = 0;

    if (sscanf(s, "%5s", word) == 1) {
	int n = 0;

	if (!strcmp(word, "if")) {
	    ci = IF;
	    n = 2;
	} else if (!strcmp(word, "else")) {
	    ci = ELSE;
	    n = 4;
	} else if (!strcmp(word, "elif")) {
	    ci = ELIF;
	    n = 4;
	} else if (!strcmp(word, "endif")) {
	    ci = ENDIF;
	    n = 5;
	}
	if (ci > 0 && s[n] != '\0' && !isspace(s[n])) {
	    ci = 0;
	}
	if (s[n] != '\0' && (ci == IF || ci == ELIF)) {
	    /* set pointer to the condition */
	    cmd_set_vstart(state->cmd, s + n + 1);
	}
    }

    return ci;
}

static int real_parse_command (ExecState *s,
			       DATASET *dset,
			       int compmode,
			       void *ptr)
{
    char *line = s->line;
    CMD *cmd = s->cmd;
    int err = 0;

#if CDEBUG
    fprintf(stderr, "real_parse_command: '%s', compmode = %d (%s)\n",
	    line, compmode, gretl_command_word(compmode));
#endif

    maybe_init_shadow();

    if (*line != '\0') {
	if (!compmode && gretl_if_state_false()) {
	    /* take a short-cut */
	    cmd->ci = get_flow_control_ci(s);
	} else {
	    /* not compiling or not blocked */
	    err = tokenize_line(s, dset, compmode);
	}

	if (!err && simple_flow_control(cmd)) {
	    /* These don't go to assemble_command(), so check
	       them here for extraneous junk.
	    */
	    err = check_for_stray_tokens(cmd);
	}

	if (err) {
	    goto parse_exit;
	}

	if (compmode) {
	    /* Are we doing get_command_index(), for compilation?
	       In that case we shouldn't do any further processing
	       unless we got a nested loop command (in which case we
	       want to extract the options), or we got a ci that
	       ought to be unitary.
	    */
	    if (compmode == LOOP && cmd->ci == LOOP) {
		err = assemble_command(cmd, dset, s, line);
		compmode = 0;
	    } else if (cmd->ci == IF || cmd->ci == ELIF) {
		err = set_command_vstart(cmd, s, s->prn);
	    }
	    goto parse_exit;
	}

	/* cmd->vstart must be set for the benefit of flow_control():
	   it will hold the condition attached to IF or ELIF
	*/
	if ((cmd->ci == IF || cmd->ci == ELIF) && cmd->vstart == NULL) {
	    err = set_command_vstart(cmd, s, s->prn);
	}

	/* If we haven't already hit an error, then we need to consult
	   and perhaps modify the flow control state -- and if we're
	   blocked, return.
	*/
	if (ptr != NULL && (cmd->flags & CMD_SUBST)) {
	    /* don't attempt to compile if string substitution
	       is going on
	    */
	    ptr = NULL;
	}
	if (!err && flow_control(s, dset, ptr)) {
	    if (cmd->err) {
		/* we hit an error evaluating the if state */
		err = cmd->err;
	    } else {
		cmd->ci = CMD_MASKED;
	    }
	    goto parse_exit;
	}

	/* Otherwise proceed to "assemble" the parsed command */
	if (!err && !simple_flow_control(cmd)) {
	    err = assemble_command(cmd, dset, s, line);
	}
    }

 parse_exit:

#if CDEBUG
    if (cmd->ci == CMD_MASKED) {
	fprintf(stderr, "breaking on flow control, current state = %s\n\n",
		gretl_if_state_false() ? "false" : "true");
    } else if (compmode) {
	fputc('\n', stderr);
    }
#endif

    if (err) {
#if CDEBUG
	fprintf(stderr, "+++ tokenizer: err=%d on '%s'\n", err, line);
#endif
	get_or_set_errline(line, 1);
    } else {
	get_or_set_errline(NULL, 1);
    }

    return err;
}

/* Here we're parsing a command line that was assembled via the gretl
   GUI (menus and dialogs). We can take some shortcuts in this case,
   since we don't have to worry about filtering comments, carrying out
   string substitution, or "if-state" conditionality.
*/

int parse_gui_command (char *line, CMD *cmd, DATASET *dset)
{
    ExecState s = {0};
    int err = 0;

    maybe_init_shadow();

    s.line = line;
    s.cmd = cmd;

    gretl_cmd_clear(cmd);
    gretl_error_clear();

    if (*line != '\0') {
	err = tokenize_line(&s, dset, 0);
	if (!err) {
	    err = assemble_command(cmd, dset, NULL, NULL);
	}
    }

    if (err) {
	fprintf(stderr, "+++ parse_gui_command: err=%d on '%s'\n",
		err, line);
    }

    return err;
}
