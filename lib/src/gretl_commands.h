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

#ifndef COMMANDS_H
#define COMMANDS_H

typedef enum {
    SEMIC = 0,    
    ADD,
    ADDTO,
    ADF,
    APPEND,
    AR,  
    ARBOND,
    ARCH,
    ARMA,
    BREAK,
    BXPLOT,
    CHOW,
    COEFFSUM,
    COINT,
    COINT2,
    CORC,     
    CORR,     
    CORRGM,   
    CRITERIA,
    CUSUM,
    DATA,
    DATAMOD,
    DELEET,
    DIFF,
    DIFFTEST,
    DISCRETE,
    DUMMIFY,
    ELSE,
    END,
    ENDIF,
    ENDLOOP,
    EQNPRINT, 
    EQUATION,
    ESTIMATE,
    FCAST,   
    FCASTERR, 
    FIT,  
    FREQ, 
    FUNC,
    FUNCERR,
    GARCH,
    GENR,  
    GMM,
    GNUPLOT, 
    GRAPH,
    HAUSMAN,
    HCCM,    
    HECKIT,
    HELP,    
    HILU,    
    HSK,
    HURST,
    IF,
    INCLUDE,
    INFO,
    KPSS,
    LABELS, 
    LAD,
    LAGS,    
    LDIFF,
    LEVERAGE,
    LMTEST, 
    LOGISTIC,
    LOGIT,
    LOGS,
    LOOP,
    MAHAL,
    MEANTEST,
    MLE,
    MODELTAB,
    MPOLS,
    MULTIPLY,
    NLS,
    NULLDATA,
    OLS,     
    OMIT,
    OMITFROM,
    OPEN,
    OUTFILE,
    PANEL,
    PCA,
    PERGM,
    PLOT,    
    POISSON,
    PRINT, 
    PRINTF,
    PROBIT,
    PVALUE, 
    PWE,
    QLRTEST,
    QUIT,
    REMEMBER,
    RENAME,
    RESET,
    RESTRICT,
    RHODIFF, 
    RMPLOT,
    RUN,
    RUNS,
    SCATTERS,
    SDIFF,
    SET,
    SETINFO,
    SETOBS,
    SETMISS,
    SHELL,  
    SMPL,
    SPEARMAN,
    SPRINTF,
    SQUARE,  
    STORE, 
    STRING,
    SUMMARY,
    SYSTEM,
    TABPRINT,
    TESTUHAT,
    TOBIT,
    TSLS,    
    VAR,
    VARLIST,
    VARTEST,
    VECM,
    VIF,
    WLS,
    XCORRGM,
    XTAB,
    NC
} GretlCmdIndex;

#define NEEDS_TWO_VARS(c)  ((c) == AR || (c) == ARCH || (c) == COINT || \
                            (c) == CORC || (c) == CORR || (c) == HCCM || \
                            (c) == HILU || (c) == HSK || (c) == LOGIT || \
                            (c) == SPEARMAN || (c) == OLS || (c) == PROBIT || \
                            (c) == TSLS || (c) == VAR || (c) == WLS || \
			    (c) == PWE || (c) == XTAB)

#define TEXTSAVE_OK(c) (c == ADD || \
                        c == ADDTO || \
                        c == ADF || \
                        c == CHOW || \
                        c == COEFFSUM || \
                        c == COINT || \
                        c == COINT2 || \
                        c == CORR || \
                        c == CORRGM || \
                        c == CRITERIA || \
                        c == CUSUM || \
                        c == FCAST || \
                        c == FCASTERR || \
                        c == FIT || \
                        c == FREQ || \
                        c == GRAPH || \
                        c == HAUSMAN || \
                        c == KPSS || \
                        c == LEVERAGE || \
                        c == LMTEST || \
                        c == MEANTEST || \
                        c == MPOLS || \
                        c == OMIT || \
                        c == OMITFROM || \
                        c == PCA || \
                        c == PERGM || \
                        c == PLOT || \
                        c == PRINT || \
                        c == PRINTF || \
                        c == PVALUE || \
                        c == QLRTEST || \
                        c == RESET || \
                        c == RMPLOT || \
                        c == RUNS || \
                        c == SPEARMAN || \
                        c == SUMMARY || \
                        c == TESTUHAT || \
                        c == VIF || \
                        c == VARTEST || \
                        c == XCORRGM || \
                        c == XTAB)

#define NEEDS_MODEL_CHECK(c) (c == ADD || \
                              c == OMIT || \
                              c == COEFFSUM || \
                              c == CUSUM || \
                              c == RESET || \
                              c == CHOW || \
                              c == QLRTEST || \
                              c == VIF || \
                              c == TABPRINT || \
                              c == EQNPRINT || \
                              c == FCAST || \
                              c == FIT || \
                              c == FCASTERR || \
                              c == HAUSMAN || \
                              c == LEVERAGE || \
                              c == LMTEST || \
                              c == TESTUHAT)

#define USES_BFGS(c) ((c) == ARMA || \
                      (c) == GARCH || \
                      (c) == GMM || \
		      (c) == HECKIT || \
                      (c) == MLE || \
                      (c) == TOBIT)

int gretl_command_number (const char *s);

const char *gretl_command_word (int i);

const char *gretl_command_complete_next (const char *s, int idx);

const char *gretl_command_complete (const char *s);

void gretl_command_hash_cleanup (void);

#endif /* COMMANDS_H */
