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
    ADDTO,
    ADF,
    AR,      
    ARCH,
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
    FCAST,   
    FCASTERR, 
    FIT,     
    FREQ,     
    GENR,     
    GNUPLOT, 
    GRAPH,
    HAUSMAN,
    HCCM,    
    HELP,    
    HILU,    
    HSK,
    IF,
    IMPORT,
    INFO,
    LABEL,
    LABELS, 
    LAD,
    LAGS,    
    LDIFF,
    LEVERAGE,
    LMTEST, 
    LOGIT,
    LOGS,
    LOOP,
    MEANTEST,
    MPOLS,
    MULTIPLY,
    MVAVG,
    NLS,
    NOECHO,
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
    POOLED,
    PRINT, 
    PRINTF,
    PROBIT,
    PVALUE,  
    QUIT, 
    RENAME,
    RESET,
    RHODIFF, 
    RMPLOT,
    RUN,
    RUNS,
    SCATTERS,
    SEED,
    SETOBS,
    SETMISS,
    SHELL,   
    SIM,     
    SMPL,
    SPEARMAN,
    SQUARE,  
    STORE,   
    SUMMARY,
    SYSTEM,
    TABPRINT,
    TESTUHAT,
    TSLS,    
    VAR,
    VARLIST,
    VARTEST,
    WLS,
    NC
} cmd_codes;

#define NEEDS_TWO_VARS(c)  ((c) == AR || (c) == ARCH || (c) == COINT || \
                           (c) == CORC || (c) == HCCM || (c) == HILU || \
                           (c) == HSK || (c) == LOGIT || (c) == MVAVG || \
                           (c) == OLS || (c) == POOLED || (c) == PROBIT || \
                           (c) == TSLS || (c) == VAR || (c) == WLS || \
                           (c) == SPEARMAN)

extern char *gretl_commands[];  /* see gretl_cmdlist.h */

#endif /* COMMANDS_H */
