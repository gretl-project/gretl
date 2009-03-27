/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

/* gretl_commands.c */

#include <stdlib.h>
#include <string.h>
#include <glib.h>

#include "gretl_commands.h"

struct gretl_cmd {
    int cnum;
    const char *cword;
};

static struct gretl_cmd gretl_cmds[] = {
    { SEMIC,    ";" },     
    { ADD,      "add" },
    { ADF,      "adf" }, 
    { APPEND,   "append" },
    { AR,       "ar" },  
    { AR1,      "ar1" },
    { ARBOND,   "arbond" },
    { ARCH,     "arch" },
    { ARMA,     "arima" },
    { BREAK,    "break" },
    { BXPLOT,   "boxplot" },
    { CHOW,     "chow" },     
    { COEFFSUM, "coeffsum" },
    { COINT,    "coint" },
    { COINT2,   "coint2" },
    { CORR,     "corr" },     
    { CORRGM,   "corrgm" },   
    { CRITERIA, "criteria" },
    { CUSUM,    "cusum" },
    { DATA,     "data" },
    { DATAMOD,  "dataset" },
    { DELEET,   "delete" },
    { DIFF,     "diff" },
    { DIFFTEST, "difftest" },
    { DISCRETE, "discrete" },
    { DUMMIFY,  "dummify" },
    { ELIF,     "elif" },
    { ELSE,     "else" },
    { END,      "end" },
    { ENDIF,    "endif" },
    { ENDLOOP,  "endloop" },
    { EQNPRINT, "eqnprint" },
    { EQUATION, "equation" },
    { ESTIMATE, "estimate" },
    { FCAST,    "fcast" }, 
    { FOREIGN,  "foreign" }, 
    { FREQ,     "freq" }, 
    { FUNC,     "function" },
    { FUNCERR,  "funcerr" },
    { GARCH,    "garch" },
    { GENR,     "genr" },  
    { GMM,      "gmm" },
    { GNUPLOT,  "gnuplot" },  
    { GRAPH,    "graph" }, 
    { HAUSMAN,  "hausman" },
    { HECKIT,   "heckit" },
    { HELP,     "help" },    
    { HSK,      "hsk" }, 
    { HURST,    "hurst" },
    { IF,       "if" },
    { INCLUDE,  "include" },
    { INFO,     "info" }, 
    { INTREG,   "intreg" }, 
    { KALMAN,   "kalman" },
    { KPSS,     "kpss" },
    { LABELS,   "labels" },
    { LAD,      "lad" },
    { LAGS,     "lags" },    
    { LDIFF,    "ldiff" },
    { LEVERAGE, "leverage" },
    { LMTEST,   "lmtest" }, 
    { LOGISTIC, "logistic" },
    { LOGIT,    "logit" },
    { LOGS,     "logs" },
    { LOOP,     "loop" },
    { MAHAL,    "mahal" },
    { MEANTEST, "meantest" },
    { MLE,      "mle" },
    { MODELTAB, "modeltab" },
    { MODPRINT, "modprint" },
    { MPOLS,    "mpols" },
    { NLS,      "nls" },
    { NORMTEST, "normtest" },
    { NULLDATA, "nulldata" },
    { OLS,      "ols" },     
    { OMIT,     "omit" },
    { OPEN,     "open" },
    { ORTHDEV,  "orthdev" },
    { OUTFILE,  "outfile" },
    { PANEL,    "panel" },
    { PCA,      "pca" },
    { PERGM,    "pergm" },
    { PLOT,     "plot" },    
    { POISSON,  "poisson" },
    { PRINT,    "print" }, 
    { PRINTF,   "printf" },
    { PROBIT,   "probit" },
    { PVALUE,   "pvalue" }, 
    { QUANTREG, "quantreg" },
    { QLRTEST,  "qlrtest" }, 
    { QUIT,     "quit" }, 
    { RENAME,   "rename" },
    { RESET,    "reset" },
    { RESTRICT, "restrict" },
    { RMPLOT,   "rmplot" },
    { RUN,      "run" },
    { RUNS,     "runs" },
    { SCATTERS, "scatters" },
    { SDIFF,    "sdiff" },
    { SET,      "set" },
    { SETINFO,  "setinfo" },
    { SETOBS,   "setobs" },
    { SETMISS,  "setmiss" },
    { SHELL,    "shell" },
    { SMPL,     "smpl" },
    { SPEARMAN, "spearman" },
    { SPRINTF,  "sprintf" },
    { SQUARE,   "square" },
    { SSCANF,   "sscanf" },
    { STORE,    "store" },   
    { SUMMARY,  "summary" },
    { SYSTEM,   "system" },
    { TABPRINT, "tabprint" },
    { TESTUHAT, "testuhat" },
    { TOBIT,    "tobit" },
    { IVREG,    "tsls", },
    { VAR,      "var" },
    { VARLIST,  "varlist" },
    { VARTEST,  "vartest" },
    { VECM,     "vecm" },
    { VIF,      "vif" },
    { WLS,      "wls" },
    { XCORRGM,  "xcorrgm" },
    { XTAB,     "xtab" },
    { NC,       NULL}
}; 

static struct gretl_cmd gretl_cmd_aliases[] = {
    { GENR, "eval" },   
    { GENR, "series" },   
    { GENR, "scalar" },  
    { GENR, "matrix" },  
    { GENR, "string" },
    { LOOP, "for" },   
    { LOOP, "foreach" },
    { LOOP, "while" },
    { ARMA, "arma" },
    { NC,   NULL }
};

const char *gretl_command_word (int i)
{
    if (i >= 0 && i < NC) {
	return gretl_cmds[i].cword;
    } else {
	return "";
    }
}

static GHashTable *ht;

static void gretl_command_hash_init (void)
{
    int i;

    ht = g_hash_table_new(g_str_hash, g_str_equal);

    for (i=0; gretl_cmds[i].cword != NULL; i++) {
	g_hash_table_insert(ht, (gpointer) gretl_cmds[i].cword, 
			    GINT_TO_POINTER(gretl_cmds[i].cnum));
    }

    for (i=0; gretl_cmd_aliases[i].cword != NULL; i++) {
	g_hash_table_insert(ht, (gpointer) gretl_cmd_aliases[i].cword, 
			    GINT_TO_POINTER(gretl_cmd_aliases[i].cnum));
    }
}

int gretl_command_number (const char *s)
{    
    gpointer p;
    int ret = 0;

    if (ht == NULL) {
	gretl_command_hash_init();
    }

    p = g_hash_table_lookup(ht, s);
    if (p != NULL) {
	ret = GPOINTER_TO_INT(p);
    }

    return ret;
}

void gretl_command_hash_cleanup (void)
{
    if (ht != NULL) {
	g_hash_table_destroy(ht);
    }
}

const char *gretl_command_complete_next (const char *s,
					 int idx)
{
    size_t n = strlen(s);
    int i;

    for (i=idx; i<NC; i++) {
	if (!strncmp(s, gretl_cmds[i].cword, n)) {
	    return gretl_cmds[i].cword;
	}
    }
    return NULL;
}  

const char *gretl_command_complete (const char *s)
{
    return gretl_command_complete_next(s, 0);
}  
