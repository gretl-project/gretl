/*
 * Copyright (c) 1989, 1993, 1994
 *	The Regents of the University of California.  All rights reserved.
 *
 * This code is derived from software contributed to Berkeley by
 * Kim Letkeman.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *	This product includes software developed by the University of
 *	California, Berkeley and its contributors.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

/* The following calendar code is based on the source code for the
   unix "cal" program, the provenance of which is given in the above
   notice.  Allin Cottrell, March 2002.
*/

#include "libgretl.h"

/* leap year -- account for gregorian reformation in 1752 */
#define	leap_year(yr) \
	((yr) <= 1752 ? !((yr) % 4) : \
	(!((yr) % 4) && ((yr) % 100)) || !((yr) % 400))

/* number of centuries since 1700, not inclusive */
#define	centuries_since_1700(yr) \
	((yr) > 1700 ? (yr) / 100 - 17 : 0)

/* number of centuries since 1700 whose modulo of 400 is 0 */
#define	quad_centuries_since_1700(yr) \
	((yr) > 1600 ? ((yr) - 1600) / 400 : 0)

/* number of leap years between year 1 and this year, not inclusive */
#define	leap_years_since_year_1(yr) \
	((yr) / 4 - centuries_since_1700(yr) + quad_centuries_since_1700(yr))

static int days_in_month[2][13] = {
    {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
};

int week1stday = 0; /* 1 for Monday, 0 for Sunday */

#define	SATURDAY 		6		/* 1 Jan 1 was a Saturday */
#define	NUMBER_MISSING_DAYS 	11		/* 11 day correction */

static int day_in_year (int day, int month, int year)
{
    int i, leap;

    leap = leap_year(year);
    for (i=1; i<month; i++)
	day += days_in_month[leap][i];
    return day;
}

static int day_in_week (int day, int month, int year)
{
    long temp;

    temp = (long)(year - 1) * 365 + leap_years_since_year_1(year - 1)
	+ day_in_year(day, month, year);
    return (((temp - 1 + SATURDAY) - NUMBER_MISSING_DAYS) % 7);
}

long get_epoch_day (const char *date)
{
    long temp;
    int year, month, day;

    if (sscanf(date, "%d/%d/%d", &year, &month, &day) != 3)
	return -1;
    if (year < 0 || month < 0 || day < 0) return -1;
    if (year > 9999 || month > 12 || day > 31) return -1;
    year += (year < 50)? 2000 : 1900;
    temp = (long)(year - 1) * 365 + leap_years_since_year_1(year - 1)
	+ day_in_year(day, month, year);
    return temp;
}

static int t_to_epoch_day (int t, long start)
{
    int startday = (((start - 1 + SATURDAY) - NUMBER_MISSING_DAYS) % 7);
    int wkends = (t + startday - 1) / 5;

    return start + t + (2 * wkends);
}

void daily_date_string (char *str, int t, const DATAINFO *pdinfo)
{
    int rem, yr;
    int add, day, mo = 0, modays = 0;
    long yrstart, dfind;

    if (pdinfo->pd == 7) 
	dfind = (long) pdinfo->sd0 + t;
    else
	dfind = t_to_epoch_day(t, (long) pdinfo->sd0);

    yr = 1 + (double) dfind / 365.248;
    
    yrstart = (long)(yr - 1) * 365 + leap_years_since_year_1(yr - 1);
    rem = dfind - yrstart;
    if (rem == 0) {
	mo = 12;
	day = 31;
    } else {
	while (modays < rem) {
	    mo++;
	    add = days_in_month[leap_year(yr)][mo];
	    if (modays + add < rem) modays += add;
	    else break;
	}
	day = rem - modays;
    }
    sprintf(str, "%02d/%02d/%02d", yr-1900, mo, day);    
}

void get_date_from_x (double x)
{
    char s[8], yrstr[3], mostr[3], daystr[3];
    int yr, mo, day;

    sprintf(s, "%.4f", x);
    sscanf(s, "%[^.].%2s%2s", yrstr, mostr, daystr);
    yr = atoi(yrstr);
    mo = atoi(mostr);
    day = atoi(daystr);
    printf("%02d/%02d/%02d\n", yr, mo, day);
}

/* The following functions have nothing to do with the "cal" program,
   but are specific to the handling of daily data so I put them in
   here */

char *missobs_vector (double **Z, const DATAINFO *pdinfo, int *misscount)
{
    int i, t;
    char *missvec = malloc(pdinfo->t2 - pdinfo->t1 + 1);
    
    if (missvec == NULL) return NULL;

    *misscount = 0;
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	missvec[t] = 0;
	for (i=1; i<pdinfo->v; i++) {
	    if (!pdinfo->vector[i]) continue;
	    if (na(Z[i][t])) {
		missvec[t] = 1;
		*misscount += 1;
		break;
	    }
	}
    }
    return missvec;
}

int undo_repack_missing (double **Z, const DATAINFO *pdinfo, 
			 const char *missvec, int misscount)
{
    int i, t, m, g;
    double *tmpmiss, *tmpgood;

    tmpmiss = malloc(misscount * sizeof *tmpmiss);
    if (tmpmiss == NULL) return 1;
    tmpgood = malloc((pdinfo->t2 - pdinfo->t1 + 1 - misscount) 
		     * sizeof *tmpgood);
    if (tmpgood == NULL) {
	free(tmpmiss);
	return 1;
    }

    for (i=1; i<pdinfo->v; i++) {
	if (!pdinfo->vector[i]) continue;
	g = 0;
	for (t=pdinfo->t1; t<=pdinfo->t2 - misscount; t++)
	     tmpgood[g++] = Z[i][t];
	m = 0;
	for (t=pdinfo->t2 + 1 - misscount; t<=pdinfo->t2; t++)
	    tmpmiss[m++] = Z[i][t];
	m = 0;
	g = 0;
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (missvec[t]) Z[i][t] = tmpmiss[m++];
	    else Z[i][t] = tmpgood[g++];
	}
    }
    free(tmpmiss);
    free(tmpgood);
    return 0;
}

int repack_missing (double **Z, const DATAINFO *pdinfo, 
		    const char *missvec, int misscount)
{
    int i, t, m, g;
    double *tmpmiss, *tmpgood;

    tmpmiss = malloc(misscount * sizeof *tmpmiss);
    if (tmpmiss == NULL) return 1;
    tmpgood = malloc((pdinfo->t2 - pdinfo->t1 + 1 - misscount) 
		     * sizeof *tmpgood);
    if (tmpgood == NULL) {
	free(tmpmiss);
	return 1;
    }

    for (i=1; i<pdinfo->v; i++) {
	if (!pdinfo->vector[i]) continue;
	m = 0;
	g = 0;
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (missvec[t]) tmpmiss[m++] = Z[i][t];
	    else tmpgood[g++] = Z[i][t];
	}
	g = 0;
	for (t=pdinfo->t1; t<=pdinfo->t2 - misscount; t++)
	     Z[i][t] = tmpgood[g++];
	m = 0;
	for (t=pdinfo->t2 + 1 - misscount; t<=pdinfo->t2; t++)
	    Z[i][t] = tmpmiss[m++];
    }
    free(tmpmiss);
    free(tmpgood);
    return 0;
}


