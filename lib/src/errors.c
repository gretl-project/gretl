/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* errors.c - error messages for gretl */

#include "libgretl.h"

int gretl_errno;
char gretl_errmsg[ERRLEN];

const char *gretl_error_messages[] = {
    NULL,
    NULL,
    "Data error.",                                            /* E_DATA = 2 */
    "Exact or near collinearity encountered.",                /* E_SINGULAR */
    "Insufficient degrees of freedom for regression.",        /* E_DF */
    "Y-prime * Y equals zero.",                               /* E_YPY */
    "Dependent variable is all zeros, aborting regression.",  /* E_ZERO */
    "Total sum of squares was not positive.",                 /* E_TSS */
    "Sum of squared residuals negative!",                     /* E_ESS */
    "First character of new name is not a letter.",           /* E_NOTALPH */
    "You can't redefine the constant in genr.",               /* E_CONST */
    "Unbalanced parentheses in genr command.",                /* E_UNBAL */
    "Too many nested parentheses.  genr not done.",           /* E_NEST */
    "Argument must be an integer.",                           /* E_NOTINTG */
    "Sorry, command not available for this estimator.",       /* E_NOTIMP */
    "Evaluation of genr expression failed.",                  /* E_IGNONZERO */
    "Unrecognized term in genr formula.",                     /* E_CASEU */
    "Unspecified error -- FIXME",                             /* E_UNSPEC */
    "Syntax error in genr formula.",                          /* E_SYNTAX */
    "\"getxvec\" failed in generating new variable.",         /* E_ER */
    "Unrecognized operator in genr formula.",                 /* E_BADOP */
    "Variable already exists.",                               /* E_VAREXISTS */
    "This command won't work with the current periodicity.",  /* E_PDWRONG */
    "The \"-o\" flag is not implemented for this command.",   /* E_OFLAG */
    "Error attempting to open file.",                         /* E_FOPEN */
    "Out of memory error.",                                   /* E_ALLOC */
    "Missing equals sign in genr.",                           /* E_NOEQ */
    "No formula supplied in genr.",                           /* E_EQN */
    "Unknown variable name in command.",                      /* E_UNKVAR */
    "The observations specified for the regression "
    "exceed those in the data set.",                          /* E_NODATA */
    "Command has insufficient arguments.",                    /* E_ARGS */
    "This command is implemented only for OLS models.",       /* E_OLSONLY */
    "Couldn't open data file.",                               /* E_DATA_DAT */
    "Invalid argument for coeff, corr, stderr, or rho.",      /* E_INVARG */
    "Invalid lag order for adf command.",                     /* E_ADF */
    "Invalid sample split for Chow test.",                    /* E_SPLIT */
    "Syntax error in command line.",                          /* E_PARSE */
    "No independent variables left after omissions.",         /* E_NOVARS */
    "No independent variables were omitted.",                 /* E_NOOMIT */
    "Can't do this: some vars in original model "
    "have been redefined.",                                   /* E_VARCHANGE */
    "No new independent variables were added.",               /* E_NOADD */
    "One or more \"added\" vars were already present.",       /* E_ADDDUP */
    "Error generating logarithms.",                           /* E_LOGS */
    "Error generating squares.",                              /* E_SQUARES */
    "Error generating lagged variables.",                     /* E_LAGS */
    "Error in auxiliary regression for rho.",                 /* E_RHO */
    "Attempting to take square root of negative number.",     /* E_SQRT */
    "Excessive exponent in genr formula.",                    /* E_HIGH */
    "Weight variable is all zeros, aborting regression.",     /* E_WTZERO */
    "Need valid starting and ending observations.",           /* E_OBS */
    "New variable name was not supplied.",                    /* E_NOVAR */
    "You must include a constant in this sort of model.",     /* E_NOCONST */
    "There were missing observations for the added "
    "variable(s).\nReset the sample and rerun the original "
    "regression first.",                                      /* E_MISS */
    "The statistic you requested is not meaningful "
    "for this model.",                                        /* E_BADSTAT */
    "Missing sub-sample information; can't merge data."       /* E_NOMERGE */
};

/**
 * get_errmsg:
 * @errcode: gretl error code (see #error_codes).
 * @errtext: pre-allocated string or NULL.
 * @prn: gretl printing struct.
 *
 * Print an error message, given an error code number.  The message
 * is printed to the string variable errtext, if it is non-NULL,
 * or otherwise to the printing struct @prn.
 * 
 * Returns: the error text string, or NULL if @errtext is NULL.
 */

char *get_errmsg (const int errcode, char *errtext, print_t *prn)
{
    if (errcode < E_MAX && gretl_error_messages[errcode]) {
	if (errtext == NULL) {
	    pprintf(prn, "%s\n", gretl_error_messages[errcode]);
	    return NULL;
	} else {
	    strcpy(errtext, gretl_error_messages[errcode]);
	    return errtext;
	}
    } else 
	return NULL;
}

/**
 * errmsg:
 * @errcode: gretl error code (see #error_codes).
 * @prn: gretl printing struct.
 *
 * Print an error message looked up from a given an error code number, 
 * or a more specific error message if available.  
 * 
 */

void errmsg (const int errcode, print_t *prn)
{
    if (gretl_errmsg[0] == '\0') 
	get_errmsg(errcode, NULL, prn);
    else pprintf(prn, "%s\n", gretl_errmsg);
}
