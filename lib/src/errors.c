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
char gretl_msg[ERRLEN];

const char *gretl_error_messages[] = {
    NULL,
    NULL,
    N_("Data error"),                                            /* E_DATA = 2 */
    N_("Exact or near collinearity encountered"),                /* E_SINGULAR */
    N_("Insufficient degrees of freedom for regression"),        /* E_DF */
    N_("Y-prime * Y equals zero"),                               /* E_YPY */
    N_("Dependent variable is all zeros, aborting regression"),  /* E_ZERO */
    N_("Total sum of squares was not positive"),                 /* E_TSS */
    N_("Sum of squared residuals negative!"),                    /* E_ESS */
    N_("First character of new name is not a letter"),           /* E_NOTALPH */
    N_("You can't redefine the constant in genr"),               /* E_CONST */
    N_("Unbalanced parentheses in genr command"),                /* E_UNBAL */
    N_("Too many nested parentheses.  genr not done"),           /* E_NEST */
    N_("Argument must be an integer"),                           /* E_NOTINTG */
    N_("Sorry, command not available for this estimator"),       /* E_NOTIMP */
    N_("Evaluation of genr expression failed"),                  /* E_IGNONZERO */
    N_("Unrecognized term in genr formula"),                     /* E_CASEU */
    N_("Unspecified error -- FIXME"),                            /* E_UNSPEC */
    N_("Syntax error in genr formula"),                          /* E_SYNTAX */
    N_("\"getxvec\" failed in generating new variable"),         /* E_ER */
    N_("Unrecognized operator in genr formula"),                 /* E_BADOP */
    N_("This command won't work with the current periodicity"),  /* E_PDWRONG */
    N_("The \"-o\" flag is not implemented for this command"),   /* E_OFLAG */
    N_("Error attempting to open file"),                         /* E_FOPEN */
    N_("Out of memory error"),                                   /* E_ALLOC */
    N_("Missing equals sign in genr"),                           /* E_NOEQ */
    N_("No formula supplied in genr"),                           /* E_EQN */
    N_("Unknown variable name in command"),                      /* E_UNKVAR */
    N_("The observations specified for the regression "
       "exceed those in the data set"),                          /* E_NODATA */
    N_("Command has insufficient arguments"),                    /* E_ARGS */
    N_("This command is implemented only for OLS models"),       /* E_OLSONLY */
    N_("Invalid argument for coeff, corr, stderr, rho, "
       "pvalue or mpow"),                                        /* E_INVARG */
    N_("Invalid lag order for adf command"),                     /* E_ADF */
    N_("Invalid sample split for Chow test"),                    /* E_SPLIT */
    N_("Syntax error in command line"),                          /* E_PARSE */
    N_("No independent variables left after omissions"),         /* E_NOVARS */
    N_("No independent variables were omitted"),                 /* E_NOOMIT */
    N_("Can't do this: some vars in original model "
       "have been redefined"),                                   /* E_VARCHANGE */
    N_("No new independent variables were added"),               /* E_NOADD */
    N_("One or more \"added\" vars were already present"),       /* E_ADDDUP */
    N_("Error generating logarithms"),                           /* E_LOGS */
    N_("Error generating squares"),                              /* E_SQUARES */
    N_("Error generating lagged variables"),                     /* E_LAGS */
    N_("Error in auxiliary regression for rho"),                 /* E_RHO */
    N_("Attempting to take square root of negative number"),     /* E_SQRT */
    N_("Excessive exponent in genr formula"),                    /* E_HIGH */
    N_("Weight variable is all zeros, aborting regression"),     /* E_WTZERO */
    N_("Need valid starting and ending observations"),           /* E_OBS */
    N_("New variable name was not supplied"),                    /* E_NOVAR */
    N_("You must include a constant in this sort of model"),     /* E_NOCONST */
    N_("There were missing observations for the added "
       "variable(s).\nReset the sample and rerun the original "
       "regression first"),                                      /* E_MISS */
    N_("The statistic you requested is not meaningful "
       "for this model"),                                        /* E_BADSTAT */
    N_("Missing sub-sample information; can't merge data")       /* E_NOMERGE */
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

char *get_errmsg (const int errcode, char *errtext, PRN *prn)
{
    if (errcode < E_MAX && gretl_error_messages[errcode]) {
	if (errtext == NULL) {
	    pprintf(prn, "%s\n", _(gretl_error_messages[errcode]));
	    return NULL;
	} else {
	    strcpy(errtext, _(gretl_error_messages[errcode]));
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

void errmsg (const int errcode, PRN *prn)
{
    if (gretl_errmsg[0] == '\0') 
	get_errmsg(errcode, NULL, prn);
    else pprintf(prn, "%s\n", gretl_errmsg);
}

int get_gretl_errno (void)
{
    return gretl_errno;
}

char *get_gretl_errmsg (void)
{
    return gretl_errmsg;
}

char *get_gretl_msg (void)
{
    if (gretl_msg[0] == 0) return NULL;
    return gretl_msg;
}
