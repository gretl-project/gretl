/*
 *  Copyright (c) by Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* libset.c for gretl */

#include "libgretl.h"
#include "internal.h"

static int use_qr;

enum {
    AUTO_LAG_STOCK_WATSON,
    AUTO_LAG_WOOLDRIDGE
};

static struct {
    int auto_lag;
    int user_lag;
    int hc_version;
} robust_opts;

int get_hc_version (void)
{
    return robust_opts.hc_version;
}

int get_hac_lag (int m)
{
    /* Variants of Newey-West */

    if (robust_opts.user_lag != 0 && robust_opts.user_lag < m - 2) {
	/* FIXME upper limit? */
	return robust_opts.user_lag;
    }

    if (robust_opts.auto_lag == AUTO_LAG_STOCK_WATSON) {
	return 0.75 * pow(m, 1.0 / 3.0);
    } else if (robust_opts.auto_lag == AUTO_LAG_WOOLDRIDGE) {
	return 4.0 * pow(m / 100.0, 2.0 / 9.0);
    }

    /* fallback -- should not be reached */
    return 0.75 * pow(m, 1.0 / 3.0);
}

int parse_set_line (const char *line, int *echo_off)
{
    char setobj[16], setarg[16];
    int nw, err = E_PARSE;

    *setobj = *setarg = '\0';

    nw = sscanf(line, "%*s %15s %15s", setobj, setarg);
    
    if (nw == 1) {
	if (!strcmp(setobj, "echo")) {
	    *echo_off = 0;
	    err = 0;
	}
    }
	    
    else if (nw == 2) {
	if (!strcmp(setobj, "echo")) {
	    if (!strcmp(setarg, "off")) {
		*echo_off = 1;
		err = 0;
	    } 
	    else if (!strcmp(setarg, "on")) {
		*echo_off = 0;
		err = 0;
	    } 
	}
	else if (!strcmp(setobj, "hac_lag")) {
	    if (!strcmp(setarg, "nw1")) {
		robust_opts.auto_lag = AUTO_LAG_STOCK_WATSON;
		robust_opts.user_lag = 0;
		err = 0;
	    }
	    else if (!strcmp(setarg, "nw2")) {
		robust_opts.auto_lag = AUTO_LAG_WOOLDRIDGE;
		robust_opts.user_lag = 0;
		err = 0;
	    }
	    else if (isdigit(*setarg)) {
		int p = atoi(setarg);

		if (p >= 0) {
		    robust_opts.user_lag = p;
		    err = 0;
		} 
	    }
	}
	else if (!strcmp(setobj, "hc_version")) {
	    if (!strcmp(setarg, "0")) {
		robust_opts.hc_version = 0;
		err = 0;
	    }
	    else if (!strcmp(setarg, "1")) {
		robust_opts.hc_version = 1;
		err = 0;
	    }
	    else if (!strcmp(setarg, "2")) {
		robust_opts.hc_version = 2;
		err = 0;
	    }
	}
	else if (!strcmp(setobj, "qr")) {
	    if (!strcmp(setarg, "on")) {
		use_qr = 1;
		err = 0;
	    }
	    else if (!strcmp(setarg, "off")) {
		use_qr = 0;
		err = 0;
	    }
	}
    }
		    
    return err;
}

/* use Cholesky or QR for regression? */

void set_use_qr (int set)
{
    use_qr = set;
}

int get_use_qr (void)
{
    if (getenv("GRETL_USE_QR")) {
	char *s = getenv("GRETL_USE_QR");

	if (*s && *s != '0') use_qr = 1;
    } 

    return use_qr;
}
