/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef COMMANDS_H
#define COMMANDS_H

typedef enum {
    SEMIC,    
    ADD,
    ADDOBS,
    ADDTO,
    ADF,
    APPEND,
    AR,      
    ARCH,
    ARIMA,
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
    CRITICAL,
    CUSUM,
    DATA,
    DELEET,
    DIFF,
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
    GNUPLOT, 
    GRAPH,
    HAUSMAN,
    HCCM,    
    HELP,    
    HILU,    
    HSK,
    HURST,
    IF,
    IMPORT,
    INCLUDE,
    INFO,
    KPSS,
    LABEL,
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
    MATRIX,
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
    POOLED,
    PRINT, 
    PRINTF,
    PROBIT,
    PVALUE, 
    PWE,
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
    SETOBS,
    SETMISS,
    SHELL,   
    SMPL,
    SPEARMAN,
    SQUARE,  
    STORE,   
    SUMMARY,
    SYSTEM,
    TABPRINT,
    TESTUHAT,
    TOBIT,
    TRANSPOSE,
    TSLS,    
    VAR,
    VARLIST,
    VARTEST,
    VECM,
    VIF,
    WLS,
    NC
} GretlCmdIndex;

#define NEEDS_TWO_VARS(c)  ((c) == AR || (c) == ARCH || (c) == COINT || \
                           (c) == CORC || (c) == HCCM || (c) == HILU || \
                           (c) == HSK || (c) == LOGIT || (c) == SPEARMAN || \
                           (c) == OLS || (c) == POOLED || (c) == PROBIT || \
                           (c) == TSLS || (c) == VAR || (c) == WLS || (c) == PWE)

#define TEXTSAVE_OK(c) (c == ADD || \
                        c == ADDTO || \
                        c == ADF || \
                        c == ARCH || \
                        c == CHOW || \
                        c == COEFFSUM || \
                        c == COINT || \
                        c == COINT2 || \
                        c == CORR || \
                        c == CORRGM || \
                        c == CRITERIA || \
                        c == CRITICAL || \
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
                        c == RESET || \
                        c == RMPLOT || \
                        c == RUNS || \
                        c == SPEARMAN || \
                        c == SUMMARY || \
                        c == TESTUHAT || \
                        c == VIF || \
                        c == VARTEST)

int gretl_command_number (const char *s);

const char *gretl_command_word (int i);

const char *gretl_command_complete_next (const char *s, int idx);

const char *gretl_command_complete (const char *s);

void gretl_command_hash_cleanup (void);

#endif /* COMMANDS_H */
