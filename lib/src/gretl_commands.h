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
    ADF,
    ANOVA,
    APPEND,
    AR,  
    AR1,
    ARBOND,
    ARCH,
    ARMA,
    BIPROBIT,
    BREAK,
    BXPLOT,
    CHOW,
    CLEAR,
    COEFFSUM,
    COINT,
    COINT2,
    CORR,     
    CORRGM,   
    CUSUM,
    DATA,
    DATAMOD,
    DELEET,
    DIFF,
    DIFFTEST,
    DISCRETE,
    DPANEL,
    DUMMIFY,
    DURATION,
    ELIF,
    ELSE,
    END,
    ENDIF,
    ENDLOOP,
    EQNPRINT, 
    EQUATION,
    ESTIMATE,
    FCAST,
    FLUSH,
    FOREIGN,
    FRACTINT,
    FREQ, 
    FUNC,
    FUNCERR,
    GARCH,
    GENR,  
    GMM,
    GNUPLOT, 
    GRAPHPG,
    HAUSMAN,
    HECKIT,
    HELP,    
    HSK,
    HURST,
    IF,
    INCLUDE,
    INFO,
    INTREG,
    JOIN,
    KALMAN,
    KPSS,
    LABELS, 
    LAD,
    LAGS,    
    LDIFF,
    LEVERAGE,
    LEVINLIN,
    LOGISTIC,
    LOGIT,
    LOGS,
    LOOP,
    MAHAL,
    MAKEPKG,
    MARKERS,
    MEANTEST,
    MLE,
    MODELTAB,
    MODPRINT,
    MODTEST,
    MPI,
    MPOLS,
    NEGBIN,
    NLS,
    NORMTEST,
    NULLDATA,
    OLS,     
    OMIT,
    OPEN,
    ORTHDEV,
    OUTFILE,
    PANEL,
    PCA,
    PERGM,
    PLOT,    
    POISSON,
    PRINT, 
    PRINTF,
    PROBIT,
    PVAL, 
    QUANTREG,
    QLRTEST,
    QQPLOT,
    QUIT,
    RENAME,
    RESET,
    RESTRICT,
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
    SSCANF,
    STORE, 
    SUMMARY,
    SYSTEM,
    TABPRINT,
    TOBIT,
    IVREG,
    VAR,
    VARLIST,
    VARTEST,
    VECM,
    VIF,
    WLS,
    XCORRGM,
    XTAB,
    FUNDEBUG,
    FUNCRET,
    CATCH,
    NC
} GretlCmdIndex;

#define NEEDS_TWO_VARS(c)  ((c) == AR || (c) == ARCH || (c) == COINT || \
                            (c) == AR1 || (c) == CORR || (c) == HSK || \
                            (c) == LOGIT || (c) == PROBIT || \
                            (c) == SPEARMAN || (c) == OLS || \
                            (c) == IVREG || (c) == WLS || \
			    (c) == XTAB)

#define TEXTSAVE_OK(c) (c == ADF || \
                        c == CHOW || \
                        c == COEFFSUM || \
                        c == COINT || \
                        c == COINT2 || \
                        c == CORR || \
                        c == CORRGM || \
                        c == CUSUM || \
                        c == FCAST || \
                        c == FREQ || \
                        c == HAUSMAN || \
                        c == KPSS || \
                        c == LEVERAGE || \
                        c == LEVINLIN || \
                        c == MODTEST || \
                        c == MEANTEST || \
                        c == MPOLS || \
			c == NORMTEST || \
                        c == PCA || \
                        c == PERGM || \
                        c == PLOT || \
                        c == PRINT || \
                        c == PRINTF || \
                        c == PVAL || \
                        c == QLRTEST || \
                        c == RESET || \
                        c == RMPLOT || \
                        c == RUNS || \
                        c == SPEARMAN || \
                        c == SUMMARY || \
                        c == VIF || \
                        c == VARTEST || \
                        c == XCORRGM || \
                        c == XTAB)

#define NEEDS_MODEL_CHECK(c) (c == ADD || \
                              c == CHOW || \
                              c == COEFFSUM || \
                              c == CUSUM || \
                              c == EQNPRINT || \
                              c == FCAST || \
                              c == HAUSMAN || \
                              c == LEVERAGE || \
                              c == MODTEST || \
                              c == OMIT || \
                              c == QLRTEST || \
                              c == RESET || \
                              c == TABPRINT || \
                              c == VIF)

#define GRAPHING_COMMAND(c) (c == GNUPLOT || \
			     c == BXPLOT ||  \
			     c == SCATTERS || \
			     c == CORRGM || \
			     c == XCORRGM || \
			     c == PERGM || \
			     c == RMPLOT || \
			     c == HURST || \
			     c == LEVERAGE || \
			     c == QQPLOT)
	
int gretl_command_number (const char *s);

const char *gretl_command_word (int i);

int word_is_genr_alias (const char *s);

const char *gretl_command_complete_next (const char *s, int idx);

const char *gretl_command_complete (const char *s);

void gretl_command_hash_cleanup (void);

#endif /* COMMANDS_H */
