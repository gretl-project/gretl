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

enum cmd_codes {
    SEMIC,    
    ADD,
    ADDTO,
    ADF,
    AR,      
    ARCH,
    BXPLOT,
    CHOW,
    COINT,
    CORC,     
    CORR,     
    CORRGM,   
    CRITERIA, 
    CUSUM,
    DELEET,
    DIFF,
    ENDLOOP,
    EQNPRINT, 
    FCAST,   
    FCASTERR, 
    FIT,     
    FREQ,     
    GENR,     
    GNUPLOT, 
    GRAPH,   
    HCCM,    
    HELP,    
    HILU,    
    HSK,   
    IMPORT,
    INFO,
    LABELS,  
    LAGS,    
    LDIFF,
    LIST,    
    LMTEST, 
    LOGIT,
    LOGS,
    LOOP,
    MEANTEST,
    MULTIPLY,
    MVAVG,
    NULLDATA,
    OLS,     
    OMIT,
    OMITFROM,
    OPEN,
    PERGM,
    PLOT,    
    POOLED,
    PRINT, 
    PROBIT,
    PVALUE,  
    QUIT,    
    RHODIFF, 
    RUN,
    RUNS,
    SCATTERS,
    SEED,
    SETOBS,
    SHELL,   
    SIM,     
    SMPL,
    SPEARMAN,
    SQUARE,  
    STORE,   
    SUMMARY,
    TABPRINT,
    TESTUHAT,
    TSLS,    
    VAR,
    VARTEST,
    WLS,
    NC
};

extern char *commands[];  /* see cmdlist.h */

#endif /* COMMANDS_H */
