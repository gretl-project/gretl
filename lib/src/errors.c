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
    char tmpstr[MAXLEN];

    switch (errcode) {

    case E_DATA:
	strcpy(tmpstr, "Data error.");
	break;
    case E_DF:
	strcpy(tmpstr, "Insufficient degrees of freedom for regression.");
	break;
    case E_NODATA:
        strcpy(tmpstr, "The observations specified for the regression "
	       "exceed those in the data set.");
	break;
    case E_EQN:
	strcpy(tmpstr, "No formula supplied in genr.");
	break;
    case E_NOVAR:
	strcpy(tmpstr, "New variable name was not supplied.");
	break;
    case E_NOEQ:
	strcpy(tmpstr, "Missing equals sign in genr.");
	break;
    case E_UNKVAR:
	strcpy(tmpstr, "Unknown variable name in command.");
	break;	
    case E_SINGULAR:
	strcpy(tmpstr, "Exact or near collinearity encountered.");
	break;
    case E_ESS:
	strcpy(tmpstr, "Sum of squared residuals negative!");
	break;
    case E_NOTALPH:
	strcpy(tmpstr, "First character of new name is not a letter.");
	break;
    case E_CONST:
	strcpy(tmpstr, "You can't redefine the constant in genr.");
	break;
    case E_INEG:
	strcpy(tmpstr, "New variable name not supplied.");
	break;
    case E_UNBAL:
	strcpy(tmpstr, "Unbalanced parentheses in genr command.");
	break;
    case E_NEST:
	strcpy(tmpstr, "Too many nested parentheses.  genr not done.");
	break;
    case E_NOTINTG:
	strcpy(tmpstr, "Argument must be an integer.");
	break;
    case E_NOTIMP:
	strcpy(tmpstr, "Sorry, command not available for this estimator.");
	break;
    case E_IGNONZERO:
	strcpy(tmpstr, "Evaluation of genr expression failed.");
	break;
    case E_ZERO:
	strcpy(tmpstr, "Dependent variable is all zeros, aborting regression.");
	break;
    case E_WTZERO:
	strcpy(tmpstr, "Weight variable is all zeros, aborting regression.");
	break;
    case E_CASEU:
	strcpy(tmpstr, "Unrecognized term in genr formula.");
	break;
    case E_UNSPEC:
	strcpy(tmpstr, "Unspecified error -- FIXME");
	break;
    case E_SYNTAX:
	strcpy(tmpstr, "Syntax error in genr formula.");
	break;
    case E_ER:
	strcpy(tmpstr, "\"getxvec\" failed in generating new variable.");
	break;
    case E_BADOP:
	strcpy(tmpstr, "Unrecognized operator in genr formula.");
	break;
    case E_VAREXISTS:
	strcpy(tmpstr, "Variable already exists.");
	break;
    case E_PDWRONG:
	strcpy(tmpstr, "This command won't work with the current periodicity.");
	break;
    case E_OFLAG:
	strcpy(tmpstr, "The \"-o\" flag is not implemented for this command.");
	break;
    case E_FOPEN:
	strcpy(tmpstr, "Error attempting to open file.");
	break;
    case E_ALLOC:
	strcpy(tmpstr, "Out of memory error.");
	break;
    case E_ARGS:
	strcpy(tmpstr, "Command has insufficient arguments.");
	break;
    case E_OLSONLY:
	strcpy(tmpstr, "This command is implemented only for OLS models.");
	break;
    case E_DATA_DAT:
	strcpy(tmpstr, "Couldn't open data file.");
	break;
    case E_INVARG:
	strcpy(tmpstr, "Invalid argument for coeff, corr, stderr, or rho.");
	break;
    case E_ADF:
	strcpy(tmpstr, "Invalid lag order for adf command.");
	break;
    case E_SPLIT:
	strcpy(tmpstr, "Invalid sample split for Chow test.");
	break;
    case E_PARSE:
	strcpy(tmpstr, "Syntax error in command line.");
	break;
    case E_NOVARS:
	strcpy(tmpstr, "No independent variables left after omissions.");
	break;
    case E_NOOMIT:
	strcpy(tmpstr, "No independent variables were omitted.");
	break;
    case E_NOADD:
	strcpy(tmpstr, "No new independent variables were added.");
	break;
    case E_ADDDUP:
	strcpy(tmpstr, "One or more \"added\" vars were already present.");
	break;
    case E_VARCHANGE:
	strcpy(tmpstr, "Can't do this: some vars in original model "
	       "have been redefined.");
	break;
    case E_LOGS:
	strcpy(tmpstr, "Error generating logarithms.");
	break;
    case E_SQUARES:
	strcpy(tmpstr, "Error generating squares.");
	break;
    case E_LAGS:
	strcpy(tmpstr, "Error generating lagged variables.");
	break;
    case E_RHO:
	strcpy(tmpstr, "Error in auxiliary regression for rho.");
	break;
    case E_HIGH:
	strcpy(tmpstr, "Excessive exponent in genr formula.");
	break;
    case E_SQRT:
	strcpy(tmpstr, "Attempting to take square root of negative number.");
	break;
    case E_TSS:
	strcpy(tmpstr, "Total sum of squares was not positive.");
	break;
    case E_OBS:
	strcpy(tmpstr, "Need valid starting and ending observations.");
	break;
    case E_NOCONST:
	strcpy(tmpstr, "You must include a constant in this sort of model.");
	break;
    case E_MISS:
	strcpy(tmpstr, "There were missing observations for the added "
	       "variable(s).\nReset the sample and rerun the original "
	       "regression first.");
	break;
    case E_BADSTAT:
	strcpy(tmpstr, "The statistic you requested is not meaningful "
	       "for this model.");
	break;
    default:
	strcpy(tmpstr, "Unclassified error");
	fprintf(stderr, "Numeric error code = %d\n", errcode); 
	break;
    }
    if (errtext == NULL) {
	pprintf(prn, "%s\n", tmpstr);
	return NULL;
    } else {
	strcpy(errtext, tmpstr);
	return errtext;
    }
}

/**
 * errmsg:
 * @errcode: gretl error code (see #error_codes).
 * @msg: pre-allocated string, or NULL.
 * @prn: gretl printing struct.
 *
 * Print an error message looked up from a given an error code number, 
 * or a supplied error message.  
 * 
 */

void errmsg (const int errcode, const char *msg, print_t *prn)
{
    if (msg == NULL || strlen(msg) == 0) 
	get_errmsg(errcode, NULL, prn);
    else pprintf(prn, "%s\n", msg);
}
