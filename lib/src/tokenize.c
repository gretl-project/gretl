#include "libgretl.h"
#include "gretl_func.h"
#include "uservar.h"

/* start of what should be in commands.h */

typedef enum {
    CI_LIST  = 1 << 0,  /* list may be present */
    CI_LLEN1 = 1 << 1,  /* list must contain exactly 1 member */
    CI_LLEN2 = 1 << 2,  /* list must contain exactly 2 members */
    CI_ORD1  = 1 << 3,  /* args start with order specifier */
    CI_ORD2  = 1 << 4,  /* args (may) end with order */
    CI_L1INT = 1 << 5,  /* first portion of list contains ints */
    CI_PARM1 = 1 << 6,  /* args start with (non-list) parameter */
    CI_PARM2 = 1 << 7,  /* may use second param */
    CI_EXPR  = 1 << 8,  /* uses "genr"-type expression */
    CI_VARGS = 1 << 9,  /* ends with varargs field */
    CI_EXTRA = 1 << 10, /* uses some special "extra" feature */
    CI_ADHOC = 1 << 11, /* needs special-purpose parser */
    CI_DOALL = 1 << 12, /* operates on all series if no list given */
    CI_NOOPT = 1 << 13, /* command never takes any options */
    CI_BLOCK = 1 << 14  /* command starts a block */
} CIFlags;

int gretl_command_get_flags (int ci);

/* end commands.h */

/* start of what should be in commands.c */

struct gretl_cmd_new {
    int cnum;
    const char *cword;
    CIFlags flags;
};

static struct gretl_cmd_new gretl_cmds_new[] = {
    { SEMIC,    ";",        0 },     
    { ADD,      "add",      CI_LIST },
    { ADF,      "adf",      CI_ORD1 | CI_LIST }, 
    { ANOVA,    "anova",    CI_LIST }, 
    { APPEND,   "append",   CI_PARM1 },
    { AR,       "ar",       CI_LIST | CI_L1INT },  
    { AR1,      "ar1",      CI_LIST },
    { ARBOND,   "arbond",   CI_LIST | CI_L1INT },
    { ARCH,     "arch",     CI_ORD1 | CI_LIST },
    { ARMA,     "arima",    CI_LIST | CI_L1INT },
    { BIPROBIT, "biprobit", CI_LIST },
    { BREAK,    "break",    0 },
    { BXPLOT,   "boxplot",  CI_LIST | CI_EXTRA },
    { CHOW,     "chow",     CI_PARM1 },
    { CLEAR,    "clear",    0 },
    { COEFFSUM, "coeffsum", CI_LIST  },
    { COINT,    "coint",    CI_ORD1 | CI_LIST },
    { COINT2,   "coint2",   CI_ORD1 | CI_LIST },
    { CORR,     "corr",     CI_LIST | CI_DOALL },     
    { CORRGM,   "corrgm",   CI_LIST | CI_LLEN1 | CI_ORD2 },   
    { CUSUM,    "cusum",    0 },
    { DATA,     "data",     CI_ADHOC }, /* special: needs whole line */
    { DATAMOD,  "dataset",  CI_PARM1 | CI_LIST | CI_PARM2 }, /* ? */
    { DELEET,   "delete",   CI_LIST },
    { DIFF,     "diff",     CI_LIST },
    { DIFFTEST, "difftest", CI_LIST | CI_LLEN2 },
    { DISCRETE, "discrete", CI_LIST },
    { DPANEL  , "dpanel",   CI_LIST | CI_L1INT },
    { DUMMIFY,  "dummify",  CI_LIST },
    { DURATION, "duration", CI_LIST },
    { ELIF,     "elif",     CI_EXPR },
    { ELSE,     "else",     0 },
    { END,      "end",      CI_PARM1 },
    { ENDIF,    "endif",    0 },
    { ENDLOOP,  "endloop",  0 },
    { EQNPRINT, "eqnprint", 0 }, /* special, handled later */
    { EQUATION, "equation", CI_LIST }, /* FIXME is CI_LIST right? */
    { ESTIMATE, "estimate", CI_PARM1 | CI_PARM2 }, /* params optional */
    { FCAST,    "fcast",    CI_ADHOC },
    { FLUSH,    "flush",    0 },
    { FOREIGN,  "foreign",  CI_PARM1 | CI_BLOCK }, 
    { FRACTINT, "fractint", CI_LIST | CI_LLEN1 | CI_ORD2 }, 
    { FREQ,     "freq",     CI_LIST | CI_LLEN1 }, 
    { FUNC,     "function", CI_ADHOC | CI_NOOPT | CI_BLOCK },
    { FUNCERR,  "funcerr",  CI_PARM1 },
    { GARCH,    "garch",    CI_LIST | CI_L1INT },
    { GENR,     "genr",     CI_EXPR },  
    { GMM,      "gmm",      CI_EXPR | CI_BLOCK },
    { GNUPLOT,  "gnuplot",  CI_LIST | CI_EXTRA },
    { GRAPHPG,  "graphpg",  CI_PARM1 }, 
    { HAUSMAN,  "hausman",  0 },
    { HECKIT,   "heckit",   CI_LIST },
    { HELP,     "help",     CI_PARM1 },    
    { HSK,      "hsk",      CI_LIST }, 
    { HURST,    "hurst",    CI_LIST | CI_LLEN1 },
    { IF,       "if",       CI_EXPR },
    { INCLUDE,  "include",  CI_PARM1 },
    { INFO,     "info",     0 }, 
    { INTREG,   "intreg",   CI_LIST },
    { JOIN,     "join",     CI_PARM1 | CI_PARM2 },
    { KALMAN,   "kalman",   CI_BLOCK },
    { KPSS,     "kpss",     CI_ORD1 | CI_LIST },
    { LABELS,   "labels",   CI_LIST | CI_DOALL },
    { LAD,      "lad",      CI_LIST },
    { LAGS,     "lags",     CI_ORD1 | CI_LIST },
    { LDIFF,    "ldiff",    CI_LIST },
    { LEVERAGE, "leverage", 0 },
    { LEVINLIN, "levinlin", CI_PARM1 | CI_LIST | CI_LLEN1 },
    { LOGISTIC, "logistic", CI_LIST }, /* + special */
    { LOGIT,    "logit",    CI_LIST },
    { LOGS,     "logs",     CI_LIST },
    { LOOP,     "loop",     CI_ADHOC }, /* ? */
    { MAHAL,    "mahal",    CI_LIST },
    { MAKEPKG,  "makepkg",  CI_PARM1 },
    { MARKERS,  "markers",  0 },
    { MEANTEST, "meantest", CI_LIST | CI_LLEN2 },
    { MLE,      "mle",      CI_EXPR | CI_BLOCK },
    { MODELTAB, "modeltab", CI_PARM1 },
    { MODPRINT, "modprint", CI_PARM1 | CI_PARM2 | CI_EXTRA },
    { MODTEST,  "modtest",  CI_ORD1 },
    { MPI,      "mpi",      CI_BLOCK },
    { MPOLS,    "mpols",    CI_LIST },
    { NEGBIN,   "negbin",   CI_LIST },
    { NLS,      "nls",      CI_EXPR | CI_BLOCK },
    { NORMTEST, "normtest", CI_LIST | CI_LLEN1 },
    { NULLDATA, "nulldata", CI_PARM1 }, /* maybe convert to ORD1? */
    { OLS,      "ols",      CI_LIST },     
    { OMIT,     "omit",     CI_LIST },
    { OPEN,     "open",     CI_PARM1 },
    { ORTHDEV,  "orthdev",  CI_LIST },
    { OUTFILE,  "outfile",  CI_PARM1 },
    { PANEL,    "panel",    CI_LIST },
    { PCA,      "pca",      CI_LIST | CI_DOALL },
    { PERGM,    "pergm",    CI_LIST | CI_LLEN1 | CI_ORD2 },
    { PLOT,     "textplot", CI_LIST },    
    { POISSON,  "poisson",  CI_LIST },
    { PRINT,    "print",    0 }, /* special: handled later */
    { PRINTF,   "printf",   CI_PARM1 | CI_VARGS },
    { PROBIT,   "probit",   CI_LIST },
    { PVAL,     "pvalue",   CI_ADHOC | CI_NOOPT }, 
    { QUANTREG, "quantreg", CI_PARM1 | CI_LIST },
    { QLRTEST,  "qlrtest",  0 }, 
    { QQPLOT,   "qqplot",   CI_LIST },
    { QUIT,     "quit",     0 }, 
    { RENAME,   "rename",   CI_PARM1 | CI_PARM2 }, /* code needs fixing */
    { RESET,    "reset",    0 },
    { RESTRICT, "restrict", CI_PARM1 | CI_BLOCK },
    { RMPLOT,   "rmplot",   CI_LIST | CI_LLEN1 },
    { RUN,      "run",      CI_PARM1 },
    { RUNS,     "runs",     CI_LIST | CI_LLEN1 },
    { SCATTERS, "scatters", CI_LIST },
    { SDIFF,    "sdiff",    CI_LIST },
    { SET,      "set",      CI_PARM1 | CI_PARM2 },
    { SETINFO,  "setinfo",  CI_LIST | CI_LLEN1 }, /* + special: handled later */
    { SETOBS,   "setobs",   CI_PARM1 | CI_PARM2 },
    { SETOPT,   "setopt",   CI_PARM1 | CI_PARM2 },
    { SETMISS,  "setmiss",  CI_PARM1 | CI_LIST | CI_DOALL },
    { SHELL,    "shell",    CI_ADHOC | CI_NOOPT },
    { SMPL,     "smpl",     CI_PARM1 | CI_PARM2 }, /* + alternate forms */
    { SPEARMAN, "spearman", CI_LIST | CI_LLEN2 },
    { SPRINTF,  "sprintf",  CI_PARM1 | CI_PARM2 | CI_VARGS },
    { SQUARE,   "square",   CI_LIST },
    { SSCANF,   "sscanf",   CI_PARM1 | CI_PARM2 | CI_VARGS },
    { STORE,    "store",    CI_PARM1 | CI_LIST | CI_DOALL },   
    { SUMMARY,  "summary",  CI_LIST | CI_DOALL },
    { SYSTEM,   "system",   CI_PARM1 | CI_BLOCK }, /* "method" param optional */
    { TABPRINT, "tabprint", 0 }, /* special, handled later */
    { TOBIT,    "tobit",    CI_LIST },
    { IVREG,    "tsls",     CI_LIST },
    { VAR,      "var",      CI_ORD1 | CI_LIST },
    { VARLIST,  "varlist",  0 },
    { VARTEST,  "vartest",  CI_LIST | CI_LLEN2 },
    { VECM,     "vecm",     CI_ORD1 | CI_LIST },
    { VIF,      "vif",      0 },
    { WLS,      "wls",      CI_LIST },
    { XCORRGM,  "xcorrgm",  CI_LIST | CI_LLEN2 | CI_ORD2 },
    { XTAB,     "xtab",     CI_LIST },
    { FUNDEBUG, "debug",    CI_PARM1 },
    { FUNCRET,  "return",   CI_EXPR },
    { CATCH,    "catch",    0 },
    { NC,       NULL,       0 }
}; 

int gretl_command_get_flags (int ci)
{
    if (ci >= 0 && ci < NC) {
	return gretl_cmds_new[ci].flags;
    } else {
	return 0;
    }
}

/* end commands.c */

/* prototype tokenizer for gretl commands, April 2010 
   partially updated October 2010; updated again August
   2013; and again August 2014.
*/

#define CDEBUG 2

enum {
    TOK_JOINED = 1 << 0, /* token is joined on the left (no space) */
    TOK_DONE   = 1 << 1  /* token has been handled */
} TokenFlags;

enum {
    C_CATCH   = 1 << 0
};

enum {
    TOK_NAME,    /* potentially valid identifier */
    TOK_ATSTR,   /* '@' plus potentially valid identifier */
    TOK_DOLSTR,  /* '$' plus potentially valid identifier */
    TOK_OPT,     /* long-form option flag */
    TOK_SOPT,    /* short-form option flag */ 
    TOK_CATCH,   /* "catch" keyword */
    TOK_QUOTED,  /* string in double quotes */
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
    TOK_OPTEQ,   /* '=' that joins option to value */
    TOK_OPTVAL,  /* value attached to option flag */
    TOK_AST,     /* single asterisk */
    TOK_SYMB,    /* symbols, not otherwise handled */
    TOK_SAVNAME  /* token holds command 'savename' */
} TokenTypes;    

typedef struct cmd_info_ cmd_info;
typedef struct cmd_token_ cmd_token;

struct cmd_info_ {
    int ci;          /* current command index */
    int err;         /* error code */
    int context;     /* for block commands, index of current context */
    int ciflags;     /* status flags pertaining to @ci */
    gretlopt opt;    /* option(s) for command */
    int flags;       /* status flags for command invocation (only C_CATCH so far)  */
    int order;       /* lag order, where appropriate */
    int aux;         /* auxiliary int (e.g. VECM rank) */
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

struct cmd_token_ {
    char *s;         /* allocated token string */
    const char *lp;  /* pointer to line position */
    char type;       /* one of TokenTypes */
    char flag;       /* zero or more of TokenFlags */
}; 

#define param_optional(c) (c == SET || c == HELP || c == RESTRICT || c == SYSTEM)
#define parm2_optional(c) (c == SET || c == GNUPLOT || c == BXPLOT || \
			   c == SETOPT || c == ESTIMATE)

#define not_catchable(c) (c == IF || c == ENDIF || c == ELIF)

#define token_joined(t) (t->flag & TOK_JOINED)
#define token_done(t)   (t->flag & TOK_DONE)

#define mark_token_done(t) (t.flag |= TOK_DONE)

#define option_type(t) (t == TOK_OPT || t == TOK_SOPT || \
			t == TOK_OPTDASH || t == TOK_OPTEQ || \
			t == TOK_OPTVAL)

#define delimited_type(t) (t == TOK_QUOTED || \
			   t == TOK_PRSTR ||  \
			   t == TOK_CBSTR ||  \
			   t == TOK_BRSTR)

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

static void cmd_info_destroy (cmd_info *c)
{
    int i;

    free(c->param);
    free(c->parm2);
    free(c->list);
    free(c->auxlist);

    for (i=0; i<c->nt_alloced; i++) {
	free(c->toks[i].s);
    }

    free(c->toks);
} 

static int cmd_info_init (cmd_info *c)
{
    int i, n = 16;
    int err = 0;

    c->ci = 0;
    c->context = 0;
    c->ciflags = 0;
    c->opt = 0;
    c->flags = 0;
    c->order = 0;
    c->aux = 0;
    c->err = 0;
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

    c->toks = malloc(n * sizeof *c->toks);
    if (c->toks == NULL) {
	return E_ALLOC;
    } 

    for (i=0; i<n; i++) {
	cmd_token_init(&c->toks[i]);
    }

    if (err) {
	cmd_info_destroy(c);
    } else {
	c->nt_alloced = n;
    }

    return err;
} 

static void cinfo_drop_context (cmd_info *c)
{
    c->opt = 0;
    c->flags = 0; 
    c->context = 0;
    *c->savename = '\0';
}

static void cmd_info_clear (cmd_info *c)
{
    cmd_token *tok;
    int i;

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

    /* note: c->context may persist */

    c->ci = 0;
    c->ciflags = 0;
    c->opt = 0;
    c->order = 0;
    c->aux = 0;
    c->err = 0;
    c->cstart = 0;
    c->ntoks = 0;
    c->vstart = NULL;

    if (!c->context || c->err) {
	/* is the conditionality right here? */
	c->opt = 0;
	c->flags = 0; 
	*c->savename = '\0';
    }

    if (c->param != NULL) {
	free(c->param);
	c->param = NULL;
    }

    if (c->parm2 != NULL) {
	free(c->parm2);
	c->parm2 = NULL;
    }   

    if (c->list != NULL) {
	free(c->list);
	c->list = NULL;
    } 

    if (c->auxlist != NULL) {
	free(c->auxlist);
	c->auxlist = NULL;
    }    
} 

static int real_add_token (cmd_info *c, const char *tok,
			   const char *lp, char type, char flag)
{
    int n = c->ntoks;
    int err = 0;

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

static int push_token (cmd_info *c, const char *tok, const char *s, 
		       int pos, char type, char flag)
{
    if (pos > 0 && !isspace(*(s-1))) {
	flag |= TOK_JOINED;
    }

    return real_add_token(c, tok, s, type, flag);
}

static int push_string_token (cmd_info *c, const char *tok, 
			      const char *s, int pos)
{
    char type = TOK_NAME;

    if (c->ntoks == 0 && !strcmp(tok, "catch")) {
	type = TOK_CATCH; 
    } else if (*tok == '$') {
	type = TOK_DOLSTR;
    } else if (*tok == '@') {
	/* handle string substitution if possible */
	const char *svar = get_string_by_name(tok + 1);

	if (svar != NULL) {
	    tok = svar;
	} else {
	    type = TOK_ATSTR; 
	}
    }

    return push_token(c, tok, s, pos, type, 0);
}

static int push_symbol_token (cmd_info *c, const char *tok, 
			      const char *s, int pos)
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

    return push_token(c, tok, s, pos, type, 0);
}

static unsigned int digit_spn (const char *s)
{
    const char *digits = "0123456789";

    return strspn(s, digits);
}

static int push_numeric_token (cmd_info *c, const char *tok, 
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

#define ldelim(c) (c == '(' || c == '{' || c == '"' || c == '[')

static int push_delimited_token (cmd_info *c, const char *tok, 
				 const char *s, int pos)
{
    char type = TOK_PRSTR;

    if (*s == '{') {
	type = TOK_CBSTR;
    } else if (*s == '[') {
	type = TOK_BRSTR;
    } else if (*s == '"') {
	type = TOK_QUOTED;
    }

    return push_token(c, tok, s, pos, type, 0);
}

static int symbol_spn (const char *s)
{
    const char *ok = "=+-/*<>?|~^!%&.,:;\\";

    if (*s == '=' && *(s+1) != '=') {
	return 1;
    }

    return strspn(s, ok);
}

/* we'll treat observation identifiers such as "1995:04"
   as numeric in this context, provided the string
   starts with a digit
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

static int matching_delim (int ltype)
{
    if (ltype == '(') {
	return ')';
    } else if (ltype == '{') {
	return '}';
    } else if (ltype == '"') {
	return '"';
    } else {
	return ']';
    }
}

static int closing_delimiter_pos (const char *s)
{
    int targ = matching_delim(*s);
    int n = 0;

    s++;

    while (*s) {
	if (*s == targ) {
	    if (targ == '"') {
		if (*(s-1) != '\\') {
		    return n;
		}
	    } else {
		return n;
	    }
	}
	s++;
	n++;
    }

    return -1;
}

/* determine the index of the first 'real command' token, beyond 
   "catch" or assignment to an object */

static int min_token_index (cmd_info *c)
{
    int pos = 0;

    if (c->toks[0].type == TOK_CATCH) {
	pos = 1;
    }

    if (c->ntoks > pos + 1 && c->toks[pos+1].type == TOK_ASSIGN) {
	pos += 2;
    }

    c->cstart = pos;

    return c->cstart;
}

/* old-style "-%c param" settings */

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
    }

    return 0;
}

/* build a string composed of tokens k1 to k2 */

static char *fuse_tokens (cmd_info *c, int k1, int k2, int n)
{
    char *ret = malloc(n);
    int i;

    if (ret == NULL) {
	c->err = E_ALLOC;
    } else {
	*ret = '\0';
	for (i=k1; i<=k2; i++) {
	    strcat(ret, c->toks[i].s);
	    mark_token_done(c->toks[i]);
	}
    }

    return ret;
}

/* merge tokens from position k1 rightward, stopping at the
   first token that is not left-joined
*/

static char *merge_toks_l_to_r (cmd_info *c, int k1)
{
    cmd_token *tok;
    int n = strlen(c->toks[k1].s) + 1;
    int i, k2 = k1;

    for (i=k1+1; i<c->ntoks; i++) {
	tok = &c->toks[i];
	if (token_joined(tok) && !token_done(tok)) {
	    n += strlen(tok->s);
	    k2 = i;
	} else {
	    break;
	}
    }

    return fuse_tokens(c, k1, k2, n);
}

/* merge tokens from position k2 leftward, stopping at the
   first token that is not left-joined
*/

static char *merge_toks_r_to_l (cmd_info *c, int k2)
{
    cmd_token *tok, *prevtok;
    int n = strlen(c->toks[k2].s) + 1;
    int i, k1 = k2;

    for (i=k2; i>c->cstart+1; i--) {
	tok = &c->toks[i];
	prevtok = &c->toks[i-1];
	if (token_joined(tok) && !token_done(prevtok)) {
	    n += strlen(prevtok->s);
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

static void mark_option_tokens (cmd_info *c)
{
    cmd_token *tok, *prevtok;
    int i;

    for (i=1; i<c->ntoks; i++) {
	tok = &c->toks[i];
	if (token_joined(tok)) {
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
		/* previous token is 'joined' */
		if (tok->type == TOK_EQUALS) {
		    tok->type = TOK_OPTEQ;
		} else if (tok->type == TOK_DASH ||
			   tok->type == TOK_NAME ||
			   tok->type == TOK_INT) {
		    tok->type = TOK_OPT; /* continuation of option */
		}
	    }
	}
    }
}

static cmd_token *next_joined_token (cmd_info *c, int i)
{
    if (i < c->ntoks - 1) {
	if (c->toks[i+1].flag & TOK_JOINED) {
	    return &c->toks[i+1];
	}
    }

    return NULL;
}

static int handle_option_value (cmd_info *c, int i, gretlopt opt,
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
	/* need but didn't get a value */
	gretl_errmsg_sprintf(_("The option '--%s' requires a parameter"), 
			     tok->s);
	c->err = E_BADOPT;
    }

    if (getval) {
	char *val = NULL;

	nexttok = next_joined_token(c, i + 1);
	if (nexttok == NULL) {
	    c->err = E_PARSE;
	} else if (delimited_type(nexttok->type)) {
	    val = gretl_strdup(nexttok->s);
	    nexttok->flag |= TOK_DONE;
	} else {
	    val = merge_toks_l_to_r(c, i + 1);
	}

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

static int assemble_option_flag (cmd_info *c, cmd_token *tok,
				 char *flag, int *pk)
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
		tok->flag |= TOK_DONE;
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

static int check_command_options (cmd_info *c)
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
    } 

    for (i=1; i<c->ntoks && !err; i++) {
	tok = &c->toks[i];
	if (token_done(tok)) {
	    continue;
	}
	if (tok->type == TOK_OPT) {
	    /* long-form option, possibly with attached value */
	    err = assemble_option_flag(c, tok, optflag, &i);
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
		    gretl_errmsg_sprintf(_("Invalid option '-%c'"), tok->s[j]);
		    err = E_BADOPT;
		} else {
		    c->opt |= opt;
		}
	    }
	}
    }

    c->ci = save_ci;

    return err;
}

/* get the 0-based token position of the 'real' command argument
   in 1-based position @k (i.e. k = 1 gets first real arg),
   or -1 if none
*/

static int real_arg_index (cmd_info *c, int k)
{
    int i, i0 = c->cstart + k;

    for (i=i0; i<c->ntoks; i++) {
	if (!(c->toks[i].flag & TOK_DONE)) {
	    return i;
	}
    }

    return -1;
}

/* get the 0-based token position of the last command token
   that is not already handled (or -1 if none)
*/

static int last_arg_index (cmd_info *c)
{
    int i;

    for (i=c->ntoks-1; i>c->cstart; i--) {
	if (!(c->toks[i].flag & TOK_DONE)) {
	    return i;
	}
    }

    return -1;
}

/* record the content of the token at @pos
   as either param or parm2, depending on @i
*/

static int token_to_param (cmd_info *c, int pos, int i)
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

/* look for, e.g., "-f filename" */

static int dash_char_index (cmd_info *c, const char *s)
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
	    if (!token_done(tok) && 
		(tok->type == TOK_NAME || tok->type == TOK_QUOTED) &&
		!token_joined(tok)) {
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

static int legacy_get_filename (cmd_info *c)
{
    int pos = dash_char_index(c, "f");

    if (pos > 0) {
	if (c->toks[pos].type == TOK_QUOTED) {
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

static int get_logistic_ymax (cmd_info *c)
{
    cmd_token *tok;
    int i;

    for (i=c->cstart+1; i<c->ntoks; i++) {
	tok = &c->toks[i];
	if (tok->type == TOK_EQUALS && token_joined(tok)) {
	    if (!strcmp(c->toks[i-1].s, "ymax")) {
		c->param = merge_toks_l_to_r(c, i-1);
		break;
	    }
	}
    }
		
    return c->err;
}

/* get command parameter in first position; may involve
   compositing tokens */

static int get_param (cmd_info *c)
{
    int pos = real_arg_index(c, 1);

    if (pos < 0) {
	if (!param_optional(c->ci)) {
	    c->err = E_ARGS;
	    fprintf(stderr, "%s: required param is missing\n", 
		    c->toks[c->cstart].s);
	}
	return c->err;
    }

    if (delimited_type(c->toks[pos].type)) {
	c->param = c->toks[pos].s;
	mark_token_done(c->toks[pos]);
    } else if (next_joined_token(c, pos) == NULL) {
	c->param = c->toks[pos].s;
	mark_token_done(c->toks[pos]);
    } else {
	c->param = merge_toks_l_to_r(c, pos);
    }

    if (c->ci == DATAMOD) {
	/* "dataset" command: remove the list-wanted flag if
	   the parameter doesn't require a list */
	if (strcmp(c->param, "sortby") && 
	    strcmp(c->param, "dsortby")) {
	    c->ciflags ^= CI_LIST;
	}
    } else if (c->ci == SMPL) {
	/* "smpl" command: drop the requirement for a second param
	   if the first is "full" */
	if (!strcmp(c->param, "full")) {
	    c->ciflags ^= CI_PARM2;
	}
    }	
		
    return c->err;
}

/* get command parameter in last position; may involve
   compositing tokens */

static int get_parm2 (cmd_info *c)
{
    int pos = last_arg_index(c);
    cmd_token *tok;

    if (pos < 0) {
	if (!parm2_optional(c->ci)) {
	    c->err = E_ARGS;
	    fprintf(stderr, "%s: required parm2 is missing\n", 
		    c->toks[c->cstart].s);
	}
	return c->err;
    }

    tok = &c->toks[pos];

    if (c->ciflags & CI_VARGS) {
	/* check for trailing comma: the token we want will
	   precede this */
	if (tok->type != TOK_COMMA) {
	    c->err = E_PARSE;
	    return c->err;
	} else {
	    tok->flag |= TOK_DONE;
	    pos--;
	}
    }

    tok = &c->toks[pos];

    if (delimited_type(tok->type)) {
	c->parm2 = tok->s;
	tok->flag |= TOK_DONE;
    } else if (!token_joined(tok)) {
	c->parm2 = tok->s;
	tok->flag |= TOK_DONE;
    } else {
	c->parm2 = merge_toks_r_to_l(c, pos);
    }

    return c->err;
}

/* legacy: handle SETINFO fields */

static int get_quoted_dash_fields (cmd_info *c, const char *s)
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
	    if (c->toks[pos].type == TOK_QUOTED) {
		token_to_param(c, pos, i+1);
		mark_token_done(c->toks[pos]);
	    } else {
		c->err = E_PARSE;
	    }
	}
    }

    return c->err;
}

static int first_arg_quoted (cmd_info *c)
{
    int pos = real_arg_index(c, 1);

    if (pos < 0) {
	return 0;
    } else {
	return c->toks[pos].type == TOK_QUOTED;
    }
}

/* count instances of list separator, ';', in the
   command line */

static int cinfo_get_sepcount (cmd_info *c)
{
    int i, n = 0;

    for (i=c->cstart+1; i<c->ntoks; i++) {
	if (c->toks[i].type == TOK_SEMIC) {
	    n++;
	}
    }

    return n;
}

/* see if the line contains a parenthesized token */

static int got_parenthetical_token (cmd_info *c)
{
    int i;

    for (i=c->cstart+1; i<c->ntoks; i++) {
	if (c->toks[i].type == TOK_PRSTR) {
	    return 1;
	}
    }

    return 0;
}

/* convert token @k to integer */

static int token_to_int (cmd_info *c, int k)
{
    cmd_token *tok = &c->toks[k];
    int ret = -1;

    if (tok->type == TOK_INT) {
	ret = atoi(tok->s);
	tok->flag |= TOK_DONE;
    } else if (tok->type == TOK_NAME) {
	double x = gretl_scalar_get_value(tok->s, &c->err);

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

/* some commands require an integer order as the first
   argument, but we also have the cases where an order
   in first argument position is optional, viz:

   modtest [ order ]
   lags [ order ; ] list 
*/

static int get_command_order (cmd_info *c)
{
    int pos = real_arg_index(c, 1);

    if (c->ci == MODTEST && pos < 0) {
	/* order is optional, not present, OK */
	return 0;
    }

    if (c->ci == LAGS && cinfo_get_sepcount(c) == 0) {
	/* order is optional, not present, OK */
	return 0;
    }

    if (pos < 0) {
	c->err = E_ARGS;
    } else {
	c->order = token_to_int(c, pos);
    }

    return c->err;
}

/* special for VECM only, so far */

static int get_vecm_rank (cmd_info *c)
{
    int pos = real_arg_index(c, 2);

     if (pos < 0) {
	c->err = E_ARGS;
    } else {
	c->aux = token_to_int(c, pos);
    }

    return c->err;
}

/* get an optional, trailing "order" field */

static int get_optional_order (cmd_info *c)
{
    int pos = last_arg_index(c);

    if (pos > 0) {
	c->order = token_to_int(c, pos);
    }

    return c->err;
}

/* stuff that can come before a command proper: the "catch" 
   keyword or assignment to a named object
*/

static int handle_command_preamble (cmd_info *c)
{
    int pos = 0;

    if (c->toks[0].type == TOK_CATCH) {
	if (not_catchable(c->ci)) {
	    printf("catch: cannot be applied to this command\n");
	    c->err = E_DATA;
	    return c->err;
	} else {
	    c->flags |= C_CATCH;
	    mark_token_done(c->toks[0]);
	    pos = 1;
	}
    }

    if (c->ntoks > pos + 1 && c->toks[pos+1].type == TOK_ASSIGN) {
	char *s = c->toks[pos].s;
	int n = strlen(s);

	if (n >= MAXSAVENAME) {
	    printf("savename is too long\n");
	    c->err = E_DATA;
	} else {
	    strcpy(c->savename, s);
	    c->toks[pos].type = TOK_SAVNAME;
	    mark_token_done(c->toks[pos]);
	    mark_token_done(c->toks[pos+1]);
	}
    }

    return c->err;
}

/* check for an undigested would-be string substitution */

static int got_atstr (cmd_info *c)
{
    int i;

    for (i=0; i<c->ntoks; i++) {
	if (c->toks[i].type == TOK_ATSTR) {
	    return 1;
	}
    }

    return 0;
}

#if CDEBUG

static void expr_line_out (cmd_info *c)
{
    c->vstart = c->toks[c->cstart].lp;

    if (c->ci != GENR && c->ci != c->context) {
	c->vstart += strlen(c->toks[c->cstart].s);
    }

    if (c->ciflags & CI_ADHOC) {
	printf("ad hoc line: '%s'\n", c->vstart);
    } else {
	printf("expr: '%s'\n", c->vstart);
    }
}

/* print line with preamble, command word, and any option flags
   omitted */

static void adhoc_line_out (cmd_info *c)
{
    cmd_token *tok;
    int i, j = 0;

    printf("ad hoc line: '");

    for (i=c->cstart+1; i<c->ntoks; i++) {
	tok = &c->toks[i];
	if (!option_type(tok->type)) {
	    if (j > 0 && !token_joined(tok)) {
		putchar(' ');
	    }
	    if (tok->type == TOK_QUOTED) {
		printf("\"%s\"", tok->s);
	    } else if (tok->type == TOK_PRSTR) {
		printf("(%s)", tok->s);
	    } else if (tok->type == TOK_CBSTR) {
		printf("{%s}", tok->s);
	    } else if (tok->type == TOK_BRSTR) {
		printf("[%s]", tok->s);
	    } else {
		printf("%s", tok->s);
	    }
	    j++;
	}
    }

    fputs("'\n", stdout);
}

static void varargs_line_out (cmd_info *c)
{
    int n = strlen(c->toks[c->ntoks-1].s);
    const char *v = c->toks[c->ntoks-1].lp;
    
    c->vstart = v + n;
    printf("varargs: '%s'\n", c->vstart);
}

static void my_printlist (const int *list, const char *s)
{
    int i;

    if (list == NULL) {
	printf("%s is NULL\n", s);
	return;
    }

    printf("%s: %d : ", s, list[0]);

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    fputs("; ", stdout);
	} else {
	    printf("%d ", list[i]);
	}
    }

    putchar('\n');
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
    } else if (tok->type == TOK_QUOTED) {
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
	strcpy(s, "option-value");
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
    } else if (tok->type == TOK_ATSTR) {
	strcpy(s, "@-variable");
    } else if (tok->type == TOK_DOLSTR) {
	strcpy(s, "$-variable");
    } else {
	strcpy(s, "regular");
    }

    if (tok->type != TOK_OPT && tok->type != TOK_SOPT &&
	tok->type != TOK_OPTEQ && tok->type != TOK_OPTVAL &&
	token_joined(tok)) {
	strcat(s, ", joined on left");
    }

    if (token_done(tok)) {
	strcat(s, ", handled");
    } else {
	strcat(s, ", not handled");
    }

    return s;
}

# endif

static void print_option_flags (gretlopt opt)
{
    const char *flags = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    gretlopt i, n = 0;

    for (i=OPT_A; i<=OPT_Y; i*=2) {
	if (opt & i) {
	    if (n > 0) {
		putchar('|');
	    }
	    printf("OPT_%c", *flags);
	    n++;
	}
	flags++;
    }
}

static void print_tokens (cmd_info *c)
{
    cmd_token *toks = c->toks;
    int i;

# if CDEBUG > 1
    int nt = c->ntoks;
    char desc[48];

    if (c->ciflags & (CI_EXPR | CI_ADHOC)) {
	printf("tokens examined so far:\n");
    } else {
	printf("tokens:\n");
    }

    for (i=0; i<nt; i++) {
	printf("%3d: '%s' (%s)\n", i, toks[i].s, tokstring(desc, toks, i));
    }
# endif

    if (c->ci != c->context) {
	if (c->flags & C_CATCH) {
	    printf("* catching errors\n");
	} else if (*c->savename != '\0') {
	    printf("* assignment to '%s'\n", c->savename);
	}
    } 

    i = c->cstart;

    if (c->ci == GENR) {
	printf("* command: genr (%d)\n", c->ci);
    } else if (c->ci == c->context) {
	printf("* context: '%s' (%d)\n", 
	       gretl_command_word(c->ci), c->ci);
    } else if (c->ci > 0) {
	printf("* command: '%s' (%d)\n", 
	       gretl_command_word(c->ci), c->ci);
    } else {
	printf("* '%s': don't understand\n", toks[i].s);
    }

    if (c->err == E_BADOPT) {
	printf("* command contains invalid option\n");
    } else if (c->opt) {
	printf("* command option(s) valid (");
	print_option_flags(c->opt);
	printf(")\n");
    }

    if (got_atstr(c)) {
	printf("* warning: unhandled @-variable\n");
    }

    if (c->param != NULL) {
	printf("* param: '%s'\n", c->param);
    }

    if (c->parm2 != NULL) {
	printf("* parm2: '%s'\n", c->parm2);
    }    

    if (c->order != 0) {
	printf("* order = %d\n", c->order);
    }

    if (c->err) {
	printf("*** error = %d (%s)\n\n", c->err, errmsg_get_with_default(c->err));
	return;
    }    

    if (c->ciflags & CI_LIST) {
	my_printlist(c->list, "list");
    } else if ((c->ciflags & CI_ADHOC) && !(c->ciflags & CI_NOOPT)) {
	/* ad hoc line with options trimmed out */
	adhoc_line_out(c);
    } else if (c->ciflags & CI_VARGS) {
	/* trailing varargs portion of command */
	varargs_line_out(c);
    } else if (c->ciflags & (CI_EXPR | CI_ADHOC)) {
	/* "genr" expression or non-trimmed ad hoc line */
	expr_line_out(c);
    }

    if (c->auxlist != NULL) {
	my_printlist(c->auxlist, "auxlist");
    }

    putchar('\n');
}

#endif /* CDEBUG */

static int check_for_stray_tokens (cmd_info *c)
{
    if (!(c->ciflags & (CI_EXPR | CI_ADHOC))) {
	cmd_token *tok;
	int i;

	for (i=0; i<c->ntoks && !c->err; i++) {
	    tok = &c->toks[i];
	    if (!token_done(tok)) {
		printf("parse error at unexpected token '%s'\n",
		       tok->s);
		c->err = E_PARSE;
	    }
	}
    }

    return c->err;
}

#if 0

/* for commands that take a list: how many separators are
   required (minsep) and how many are permitted (maxsep)?
*/

static void get_separator_range (int ci, int *minsep, int *maxsep)
{
    /* default: no separators */
    *minsep = *maxsep = 0;

    switch (ci) {
    case AR:
    case GARCH:
    case HECKIT:
    case SCATTERS:
	*minsep = *maxsep = 1;
	break;
    case ARBOND:
    case ARMA:
    case IVREG:
	*minsep = 1;
	*maxsep = 2;
	break;
    case COINT2:
    case VECM:
	*maxsep = 2;
	break;
    case BIPROBIT:
    case DURATION:
    case EQUATION:
    case MPOLS:
    case NEGBIN:
    case POISSON:
    case VAR:
    case XTAB:
	*maxsep = 1;
	break;
    default:
	break;
    }
}

#endif

static int is_varname (const char *s, DATASET *dset)
{
    if (current_series_index(dset, s) >= 0) {
	return 1;
    } else if (gretl_is_user_var(s)) {
	return 1;
    } else {
	return 0;
    }
}

/* determine if the command-line is a "genr"-type expression,
   which will be directed to a separate parser -- this gets
   invoked if we haven't been able to find a recognizable
   comand-word
*/

static int test_for_genr (cmd_info *c, int i, DATASET *dset)
{
    cmd_token *toks = c->toks;
    char *s = toks[i].s;
    
    if (is_varname(s, dset)) {
	c->ci = GENR;
    } else if (toks[i].type == TOK_NAME && c->ntoks > i + 1) {
	cmd_token *nexttok = &toks[i+1];

	if ((nexttok->type == TOK_DPLUS ||
	     nexttok->type == TOK_DDASH) &&
	    token_joined(nexttok)) {
	    /* increment/decrement */
	    c->ci = GENR;
	} else if (nexttok->type == TOK_EQUALS || nexttok->type == TOK_EQMOD) {
	    /* equals sign in second place */
	    c->ci = GENR;
	} else if (nexttok->type == TOK_BRSTR && token_joined(nexttok)) {
	    /* assignment to array or series element(s) */
	    c->ci = GENR;
	} else if (!strcmp(s, "list") && c->ntoks > i + 1) {
	    if (c->ntoks == i + 2) {
		/* bare declaration, as in "list L" */
		c->ci = GENR;
	    } else {
		/* may be "list L = ..." or similar */
		nexttok = &toks[i+2];
		if (nexttok->type == TOK_EQUALS || nexttok->type == TOK_EQMOD) {
		    c->ci = GENR;
		}
	    }
	} else if (function_lookup(s) || get_user_function_by_name(s)) {
	    /* function call, no assignment */
	    c->ci = GENR;
	    c->opt |= OPT_O;
	}
    }

    return c->ci;
}

/* bodge for testing */
static void deprecate_alias (const char *bad, const char *good, int cmd);

static int try_for_command_alias (const char *s, cmd_info *cinfo)
{
    int ci = 0;

    if (!strcmp(s, "exit")) {
	ci = QUIT;
	cinfo->opt = OPT_X;
    } else if (!strcmp(s, "ls")) {
	ci = VARLIST;
    } else if (!strcmp(s, "pooled")) {
	deprecate_alias("pooled", "ols", 1);
	ci = OLS;
    } else if (!strcmp(s, "equations")) {
	ci = EQUATION;
	cinfo->opt |= OPT_M;
    } else if (*s == '!' || !strcmp(s, "launch")) {
	ci = SHELL;
    } else if (!strcmp(s, "fcasterr")) {
	deprecate_alias("fcasterr", "fcast", 1);
	ci = FCAST; 	 
    } else if (!strcmp(s, "continue")) {
	ci = FUNDEBUG;
	cinfo->opt |= OPT_C;
    } else if (!strcmp(s, "next")) {
	ci = FUNDEBUG;
	cinfo->opt |= OPT_N;
    } else if (!strcmp(s, "undebug")) {
	ci = FUNDEBUG;
	cinfo->opt |= OPT_Q;
    }

    return ci;
}

/* if we have enough tokens ready, try to determine the
   current command index */

static int try_for_command_index (cmd_info *cinfo,
				  const char *s,
				  DATASET *dset)
{
    cmd_token *toks = cinfo->toks;
    int i = min_token_index(cinfo);

    if (i >= cinfo->ntoks) {
	/* not ready to test */
	return 0;
    }

    if (*s == '(') {
	; /* a function call? */
    } else if (cinfo->ntoks > i+1 && (toks[i+1].flag & TOK_JOINED)) {
	; /* a function call? */
    } else {
	cinfo->ci = gretl_command_number(toks[i].s);
	if (cinfo->ci == 0) {
	    cinfo->ci = try_for_command_alias(toks[i].s, cinfo);
	}
    }

#if CDEBUG > 1
    if (cinfo->ci > 0) {
	printf("try_for_command_index: ci = %d (%s)\n", cinfo->ci,
	       gretl_command_word(cinfo->ci));
    }
#endif

    if (cinfo->context == FOREIGN && cinfo->ci != END) {
	/* do not attempt to parse! */
	cinfo->ci = FOREIGN;
    }

    if (cinfo->ci <= 0 && cinfo->context) {
	/* FIXME is this placement right? */
	cinfo->ci = cinfo->context;
    }

    if (cinfo->ci <= 0 && cinfo->ntoks < 4) {
	cinfo->ci = test_for_genr(cinfo, i, dset);
    }

    if (cinfo->ci > 0) {
	mark_token_done(toks[i]);
	if (cinfo->ci == cinfo->context) {
	    cinfo->ciflags = CI_EXPR;
	} else {
	    cinfo->ciflags = gretl_command_get_flags(cinfo->ci);
	    if (cinfo->ci == EQUATION && (cinfo->opt & OPT_M)) {
		cinfo->ciflags ^= CI_LIST;
		cinfo->ciflags |= CI_ADHOC;
	    }
	}
    }

    return cinfo->ci;
}

/* tricky, because "list" can be a type-word, but by itself
   it's a venerable alias for "varlist"; and in addition it
   can start "list <listname> delete"
*/

static int recognize_list_command (cmd_info *c)
{
    int k = c->cstart;

    if (c->ntoks == k + 1) {
	if (!strcmp(c->toks[k].s, "list")) {
	    mark_token_done(c->toks[k]);
	    c->ci = VARLIST;
	    return 1;
	}
    } else if (c->ntoks == k + 3) {
	const char *cword = c->toks[k+2].s;

	if (!strcmp(c->toks[k].s, "list") && 
	    (!strcmp(cword, "delete") ||
	     !strcmp(cword, "print"))) {
	    mark_token_done(c->toks[k]);
	    mark_token_done(c->toks[k+2]);
	    c->ci = strcmp(cword, "print") ? DELEET : PRINT;
	    c->ciflags = CI_PARM1;
	    c->opt = OPT_L;
	    return 1;
	}
    }	

    return 0;
}

/* get a count of any tokens not marked as 'done' */

static int count_remaining_toks (cmd_info *c)
{
    int i, n = 0;

    for (i=c->cstart+1; i<c->ntoks; i++) {
	if (!(c->toks[i].flag & TOK_DONE)) {
	    n++;
	}
    }

    return n;
}

static int process_auxlist_term (cmd_info *c, cmd_token *tok, int j)
{
    int err = 0;

    if (c->auxlist == NULL) {
	if (tok->type == TOK_CBSTR) {
	    c->auxlist = gretl_list_from_string(tok->s, &err);
	} else {
	    gretl_matrix *m = get_matrix_by_name(tok->s);

	    c->auxlist = gretl_list_from_vector(m, &err);
	}
    } else if (j == 1 && c->ci == ARMA) {
	int *aux2 = NULL;

	if (tok->type == TOK_CBSTR) {
	    aux2 = gretl_list_from_string(tok->s, &err);
	} else {
	    gretl_matrix *m = get_matrix_by_name(tok->s);

	    aux2 = gretl_list_from_vector(m, &err);
	}
	if (!err) {
	    int *tmp;

	    tmp = gretl_lists_join_with_separator(c->auxlist, aux2);
	    if (tmp == NULL) {
		c->err = E_ALLOC;
	    } else {
		free(c->auxlist);
		c->auxlist = tmp;
	    }
	}
	free(aux2);
    } else {
	/* can't have more than two such terms */
	err = E_PARSE;
    }

    return err;
}

static int try_auxlist_term (cmd_info *c, cmd_token *tok, int scount)
{
    if (c->ci != ARMA && c->ci != DPANEL && c->ci != VAR) {
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

static int panel_gmm_special (cmd_info *c, const char *s)
{
    if (c->ci == ARBOND || c->ci == DPANEL) {
	if (!strcmp(s, "GMM") || !strcmp(s, "GMMlevel")) {
	    return 1;
	}
    }

    return 0;
}

static void rejoin_list_toks (cmd_info *c, int k1, int *k2,
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
	    strcat(lstr, "(");
	    strcat(lstr, tok->s);
	    strcat(lstr, ")");
	} else if (i < c->ntoks - 1 && panel_gmm_special(c, tok->s)) {
	    cmd_token *next = &c->toks[i+1];
	    char *tmp;

	    tmp = gretl_strdup_printf("%s(%s)", tok->s, next->s);
	    c->param = gretl_str_expand(&c->param, tmp, " ");
	    free(tmp);
	    mark_token_done(c->toks[i]);
	    *k2 = ++i;
	} else {
	    strcat(lstr, tok->s);
	}
	mark_token_done(c->toks[i]);
    }
}

static int cinfo_process_command_list (cmd_info *c, DATASET *dset)
{
    char lstr[MAXLINE];
    cmd_token *tok;
    int *ilist = NULL;
    int *vlist = NULL;
    int want_ints = 0;
    int scount = 0;
    int mcount = 0;
    int i, j, k, ns;

    if (c->ciflags & CI_L1INT) {
	want_ints = 1;
    }

    ns = cinfo_get_sepcount(c);

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
	    if (want_ints) {
		if (try_auxlist_term(c, tok, scount)) {
		    /* a matrix-style entry */
		    c->err = process_auxlist_term(c, tok, mcount++);
		    if (!c->err) {
			gretl_list_append_term(&ilist, 0); /* placeholder */
			tok->flag |= TOK_DONE;
		    }
		} else {
		    k = gretl_int_from_string(tok->s, &c->err);
		    if (!c->err) {
			gretl_list_append_term(&ilist, k);
			tok->flag |= TOK_DONE;
		    }
		}
	    } else if (next_joined_token(c, i) != NULL) {
		rejoin_list_toks(c, i, &i, lstr, j++);
	    } else {
		if (j > 0) {
		    strcat(lstr, " ");
		}
		strcat(lstr, tok->s);
		tok->flag |= TOK_DONE;
		j++;
	    }
	}
	if (j == 1 && (c->ciflags & CI_LLEN1)) {
	    break;
	} else if (j == 2 && (c->ciflags & CI_LLEN2)) {
	    break;
	}
    }

    printf("cinfo_process_command_list (err=%d): lstr='%s'\n", c->err, lstr);

    if (dset != NULL && *lstr != '\0') {
	vlist = generate_list(lstr, dset, &c->err);
	if (c->err && (c->ci == PRINT || c->ci == DELEET)) {
	    /* the terms may be names of non-series variables */
	    c->ciflags ^= CI_LIST;
	    c->ciflags &= ~CI_DOALL;
	    c->ciflags |= CI_ADHOC;
	    c->err = 0;
	}
    }

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

    free(ilist);
    free(vlist);

    printf("cinfo_process_command_list: returning err = %d\n", c->err);

    return c->err;
}

static int handle_command_extra (cmd_info *c)
{
    cmd_token *tok;
    int i;

    if (c->ci == GNUPLOT || c->ci == BXPLOT) {
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
	if (c->ntoks - c->cstart - 1 == 3) {
	    /* got three arguments */
	    char *extra;

	    tok = &c->toks[c->ntoks - 1];
	    if (!token_done(tok) && tok->type == TOK_NAME) {
		tok->flag |= TOK_DONE;
		extra = gretl_strdup(tok->s);
		if (extra == NULL) {
		    c->err = E_ALLOC;
		} else {
		    c->opt |= OPT_A;
		    c->err = push_option_param(c->ci, OPT_A, extra);
		}
	    }
	}
    }

    return c->err;
}

/* For a command that ends with varargs, do we have the required
   leading non-vararg tokens? This means either one or two 
   regular tokens, followed by a comma.
*/

static int got_param_tokens (cmd_info *c)
{
    int mintoks = (c->ciflags & CI_PARM2)? 3 : 2;
    int n = c->ntoks - c->cstart - 1;

    if (c->ci == PRINTF) {
	/* pre-handle the trailing comma */
	int i = c->cstart + 2;

	if (i < c->ntoks && c->toks[i].type == TOK_COMMA) {
	    c->toks[i].flag |= TOK_DONE;
	}
    } else if (c->ci == SSCANF) {
	/* backward compatibility for comma after varname */
	int i = c->cstart + 2;

	if (i < c->ntoks && c->toks[i].type == TOK_COMMA) {
	    c->toks[i].flag |= TOK_DONE;
	    mintoks = 4;
	}
    }

    return n == mintoks;
}

static int cinfo_check_end_command (cmd_info *c)
{
    int endci = gretl_command_number(c->param);

    if (endci != c->context) {
	printf("end: invalid parameter '%s'\n", c->param);
	c->err = E_DATA;
    }

    return c->err;
}

static int check_for_list (cmd_info *c)
{
    if (c->ciflags & CI_DOALL) {
	return 0; /* no list is OK */
    }

    if (c->list == NULL) {
	if (c->ci == OMIT && (c->opt & OPT_A)) {
	    ; /* auto-omit, OK */
	} else {
	    c->err = E_ARGS;
	}
    }

    if (!c->err && c->list != NULL) {
	/* check for duplicated variables */
	int dupv = gretl_list_duplicates(c->list, c->ci);

	if (dupv >= 0) {
	    printlist(c->list, "command with duplicate(s)");
	    c->err = E_DATA;
	    gretl_errmsg_sprintf(_("variable %d duplicated in the "
				   "command list."), dupv);
	} 
    }
    
    return c->err;
}

#define MAY_START_NUMBER(c) (c == '.' || c == '-' || c == '+')

static int tokenize_line (cmd_info *cinfo, const char *line,
			  DATASET *dset)
{
    char tok[FN_NAMELEN];
    const char *s = line;
    char *vtok;
    int n, m, pos = 0;
    int err = 0;

#if CDEBUG
    printf("*** line = '%s'\n", line);
#endif

    while (!err && *s) {
	int skipped = 0;

	*tok = '\0';

	if (*s == '#') {
	    break;
	} else if (isalpha(*s) || *s == '$' || *s == '@') {
	    /* regular, dollar or string identifier */
	    n = 1 + namechar_spn(s+1);
	    m = (n < FN_NAMELEN)? n : FN_NAMELEN - 1;
	    strncat(tok, s, m);
	    err = push_string_token(cinfo, tok, s, pos);
	} else if (ldelim(*s)) {
	    /* left-hand delimiter that needs to be paired */
	    n = closing_delimiter_pos(s);
	    if (n < 0) {
		fprintf(stderr, "missing closing delimiter\n");
		err = E_PARSE;
	    } else if (n < FN_NAMELEN) {
		strncat(tok, s+1, n);
		err = push_delimited_token(cinfo, tok, s, pos);
	    } else {
		vtok = gretl_strndup(s+1, n);
		if (vtok == NULL) {
		    err = E_ALLOC;
		} else {
		    err = push_delimited_token(cinfo, vtok, s, pos);
		    free(vtok);
		}
	    }
	    n += 2;
	} else if ((n = symbol_spn(s)) > 0) {
	    if (n == 1 && MAY_START_NUMBER(*s) && isdigit(*(s+1))) {
		n = numeric_spn(s, 0);
		m = (n < FN_NAMELEN)? n : FN_NAMELEN - 1;
		strncat(tok, s, m);
		err = push_numeric_token(cinfo, tok, s, pos);
	    } else {
		/* operator / symbol */
		m = (n < FN_NAMELEN)? n : FN_NAMELEN - 1;
		strncat(tok, s, m);
		err = push_symbol_token(cinfo, tok, s, pos);
	    }
	} else if (isdigit(*s)) {
	    /* numeric string */
	    n = numeric_spn(s, 1);
	    m = (n < FN_NAMELEN)? n : FN_NAMELEN - 1;
	    strncat(tok, s, m);
	    err = push_numeric_token(cinfo, tok, s, pos);
	} else if (isspace(*s)) {	    
	    n = 1;
	    skipped = 1;
	} else {
	    printf("Unexpected symbol '%c'\n", *s);
	    err = E_PARSE;
	}

	if (err) {
	    break;
	}

	if (!skipped && cinfo->ci == 0 && cinfo->ntoks > 0) {
	    /* use current info to determine command index? */
	    try_for_command_index(cinfo, s + n, dset);
	}

	if (cinfo->ciflags & CI_EXPR) {
	    /* the remainder of line will be parsed elsewhere */
	    break;
	} else if ((cinfo->ciflags & CI_ADHOC) && (cinfo->ciflags & CI_NOOPT)) {
	    /* ditto */
	    break;
	}

	if ((cinfo->ciflags & CI_VARGS) && got_param_tokens(cinfo)) {
	    /* remaining args will be parsed elsewhere: save a
	       reference to where they start */
	    break;
	} 

	s += n;
	pos += n;
    }

    if (err) {
	cmd_info_clear(cinfo);
    } 

    return err;
}

static int assemble_command (cmd_info *cinfo, DATASET *dset)
{
    if (cinfo->ntoks == 0) {
	return cinfo->err;
    }

#if CDEBUG
    printf("doing assemble_command...\n");
#endif

    if (!(cinfo->ciflags & (CI_EXPR | CI_NOOPT))) {
	cinfo->err = check_command_options(cinfo);
    }

    if (!cinfo->err) {
	handle_command_preamble(cinfo);
    }

    if (!cinfo->err && cinfo->ci == 0) {
	recognize_list_command(cinfo);
    }

    if (cinfo->err) {
	goto bailout;
    }

    /* option-related inflections of various commands */

    if (cinfo->ci == PRINT) {
	if (first_arg_quoted(cinfo)) {
	    /* printing string literal */
	    cinfo->ciflags |= CI_PARM1;
	} else {
	    /* assume for now that we're printing series */
	    cinfo->ciflags |= (CI_LIST | CI_DOALL);
	}
    } else if (cinfo->ci == BXPLOT) {
	if (got_parenthetical_token(cinfo)) {
	    /* using specials such as "foo(x=1)" */
	    cinfo->ciflags = CI_ADHOC;
	} 
    } else if (cinfo->ci == XTAB) {
	if (cinfo->opt & OPT_M) {
	    /* using named matrix, ints ok in list */
	    cinfo->ciflags |= CI_L1INT;
	}
    } else if (matrix_data_option(cinfo->ci, cinfo->opt)) {
	/* using matrix argument, plain ints ok in list */
	cinfo->ciflags |= CI_L1INT;
    } else if (cinfo->ci == SMPL) {
	if (cinfo->opt & (OPT_M | OPT_C)) {
	    /* no-missing or contiguous */
	    cinfo->ciflags = CI_LIST | CI_DOALL;
	} else if (cinfo->opt & OPT_R) {
	    /* restrict */
	    cinfo->ciflags = CI_ADHOC;
	} else if (cinfo->opt & OPT_F) {
	    /* full: no args */
	    cinfo->ciflags = 0;
	} else if (cinfo->opt & OPT_O) {
	    /* using dummy variable */
	    cinfo->ciflags = CI_LIST | CI_LLEN1;
	} else if (cinfo->opt & OPT_N) {
	    /* random sample */
	    cinfo->ciflags = CI_PARM1; /* or CI_ORD1 ? */
	}
    } else if (cinfo->ci == SET) {
	if (cinfo->opt & (OPT_F | OPT_T)) {
	    /* from file, to file */
	    cinfo->ciflags = 0;
	}
    }

    /* legacy stuff */

    if (cinfo->ci == TABPRINT || cinfo->ci == TABPRINT) {
	legacy_get_filename(cinfo);
    } else if (cinfo->ci == SETINFO) {
	get_quoted_dash_fields(cinfo, "dn");
    } else if (cinfo->ci == LOGISTIC) {
	get_logistic_ymax(cinfo);
    } 

    if (cinfo->err) {
	goto bailout;
    }

    /* main command assembly begins */

    if (cinfo->ciflags & CI_PARM1) {
	get_param(cinfo);
    } else if (cinfo->ciflags & CI_ORD1) {
	get_command_order(cinfo);
    }

    if (!cinfo->err && cinfo->ci == VECM) {
	get_vecm_rank(cinfo);
    }

    if (!cinfo->err && (cinfo->ciflags & CI_EXTRA)) {
	handle_command_extra(cinfo);
    }	

    if (!cinfo->err && (cinfo->ciflags & CI_PARM2)) {
	get_parm2(cinfo);
    }

    if (!cinfo->err && (cinfo->ciflags & CI_LIST)) {
	if (count_remaining_toks(cinfo) > 0) {
	    cinfo_process_command_list(cinfo, dset);
	}
    }

    if (!cinfo->err && (cinfo->ciflags & CI_ORD2)) {
	get_optional_order(cinfo);
    }

    if (!cinfo->err && (cinfo->ciflags & CI_LIST)) {
	printf("assemble: doing check_for_list\n");
	check_for_list(cinfo);
    }

    if (!cinfo->err && cinfo->ci == END) {
	cinfo_check_end_command(cinfo);
    }

    if (!cinfo->err) {
	check_for_stray_tokens(cinfo);
    }

 bailout:

    print_tokens(cinfo);

    if (cinfo->ci == END) {
	/* FIXME placement */
	cinfo_drop_context(cinfo);
    } else if (cinfo->ciflags & CI_BLOCK) {
	cinfo->context = cinfo->ci;
    }

    return cinfo->err;
}

int test_tokenize (char *line, CMD *cmd, DATASET *dset, void *ptr)
{
    static cmd_info cinfo;
    static int initted;
    int err = 0;

    if (line == NULL) {
	/* quit signal */
	if (initted) {
	    cmd_info_destroy(&cinfo);
	    initted = 0;
	}
	return 0;
    }    

    if (!initted) {
	err = cmd_info_init(&cinfo);
	initted = 1;
    }

    if (!err && *line != '\0') {
	err = tokenize_line(&cinfo, line, dset);
	if (!err) {
	    err = assemble_command(&cinfo, dset);
	}
    }

    if (err) {
	printf("+++ shadow: err = %d on '%s'\n", err, line);
    }

    cmd_info_clear(&cinfo);

    return err;
}
