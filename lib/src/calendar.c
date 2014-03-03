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

/* The following calendar code is partly based on the source code for
   the unix "cal" program, the provenance of which is given in the above
   notice.  Allin Cottrell, March 2002.
*/

#include "libgretl.h"

/**
 * SECTION:calendar
 * @short_description: functions for working with dates
 * @title: Calendar
 * @include: libgretl.h
 *
 * Here we have various functions dealing with calendar dates;
 * for the most part these are designed for handling daily
 * time-series data. Most of the basic code is derived from 
 * the unix "cal" program which is released under the BSD 
 * license (reproduced in the libgretl source file calendar.c).
 */

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
    int i, leap = leap_year(year);

    for (i=1; i<month; i++) {
	day += days_in_month[leap][i];
    }

    return day;
}

static int day_of_week_from_ymd (int yr, int mo, int day)
{
    int c, d;

    /* Uspensky and Heaslet, Elementary Number Theory (1939) */

    if (mo < 3) {
	yr--;
	mo += 10;
    } else {
	mo -= 2;
    }

    c = yr / 100;
    d = yr % 100;

    return ((day % 7) + ((int) floor(2.6 * mo - 0.2) % 7) + 
	    (d % 7) + ((int) floor(d / 4.0) % 7) + ((int) floor(c / 4.0) % 7)
	    - ((2 * c) % 7)) % 7; 
}

/**
 * epoch_day_from_ymd:
 * @y: year (1 <= y <= 9999).
 * @m: month (1 <= m <= 12).
 * @d: day of month (1 <= d <= 31).
 * 
 * Returns: the epoch day number, which equals 1 for the first of
 * January in the year 1 AD, or -1 on error.
 */

long epoch_day_from_ymd (int y, int m, int d)
{
    long ret;

    if (y < 1 || m < 1 || d < 1) {
	return -1;
    }

    if (y > 9999 || m > 12 || d > 31) {
	return -1;
    }

    ret = (long)(y - 1) * 365 + leap_years_since_year_1(y - 1)
	+ day_in_year(d, m, y);

    return ret;
}

static int real_ymd_from_epoch_day (long ed, int *py, int *pm, int *pd)
{
    int y, m, d, delta;
    long days = 0L;

    if (ed < 1) {
	return E_INVARG;
    }

    for (y=1; ; y++) {
	delta = leap_year(y) ? 366 : 365;
	if (days + delta <= ed) {
	    days += delta;
	} else {
	    break;
	}
    }

    if (days == ed) {
	if (y > 1) y--;
	m = 12;
	d = 31;
    } else {
	int i = leap_year(y);

	for (m=1; ; m++) {
	    delta = days_in_month[i][m];
	    if (days + delta <= ed) {
		days += delta;
	    } else {
		break;
	    }
	}
	if (days == ed) {
	    m--;
	    d = days_in_month[i][m];
	} else {
	    d = ed - days;
	}
    }

    *py = y; *pm = m; *pd = d;

    return 0;
}

/**
 * ymd_extended_from_epoch_day:
 * @ed: epoch day (ed >= 1).
 * @err: location to receive error code.
 * 
 * Returns: a string on the pattern YYYY-MM-DD (ISO 8601 extended 
 * date format) given the epoch day number, which equals 1 for the
 * first of January in the year 1 AD, or NULL on error.
 */

char *ymd_extended_from_epoch_day (long ed, int *err)
{
    char *ret = NULL;
    int y, m, d;
    int myerr;

    myerr = real_ymd_from_epoch_day(ed, &y, &m, &d);

    if (!err) {
	ret = calloc(12, 1);
	if (ret == NULL) {
	    myerr = E_ALLOC;
	} else {
	    sprintf(ret, "%04d-%02d-%02d", y, m, d);
	}
    }

    if (err != NULL) {
	*err = myerr;
    }

    return ret;
}

/**
 * ymd_basic_from_epoch_day:
 * @ed: epoch day (ed >= 1).
 * @err: location to receive error code.
 * 
 * Returns: an 8-digit number on the pattern YYYYMMDD (ISO 8601 basic 
 * date format) given the epoch day number, which equals 1 for the
 * first of January in the year 1 AD, or #NADBL on error.
 */

double ymd_basic_from_epoch_day (long ed, int *err)
{
    double x = NADBL;
    int y, m, d;

    *err = real_ymd_from_epoch_day(ed, &y, &m, &d);

    if (!*err) {
	x = 10000*y + 100*m + d;
    }

    return x;
}

/**
 * weekday_from_epoch_day:
 * @ed: epoch day (ed >= 1).
 * 
 * Returns: the weekday (Sunday = 0) corrsponding to @ed.
 */

int weekday_from_epoch_day (long ed)
{
    int y, m, d, delta;
    long days = 0L;

    for (y=1; ; y++) {
	delta = leap_year(y) ? 366 : 365;
	if (days + delta <= ed) {
	    days += delta;
	} else {
	    break;
	}
    }

    if (days == ed) {
	if (y > 1) y--;
	m = 12;
	d = 31;
    } else {
	int i = leap_year(y);

	for (m=1; ; m++) {
	    delta = days_in_month[i][m];
	    if (days + delta <= ed) {
		days += delta;
	    } else {
		break;
	    }
	}
	if (days == ed) {
	    m--;
	    d = days_in_month[i][m];
	} else {
	    d = ed - days;
	}
    }

    return day_of_week_from_ymd(y, m, d);
}

/**
 * get_epoch_day:
 * @date: string representation of calendar date, in form
 * YY[YY]-MM-DD.
 * 
 * Returns: the epoch day number, or -1 on failure.
 */

long get_epoch_day (const char *date)
{
    int nf, year, month, day;
    long ret;

    nf = sscanf(date, YMD_READ_FMT, &year, &month, &day);

    if (nf != 3) {
	/* backward compatibility: try slashes instead */
	nf = sscanf(date, "%d/%d/%d", &year, &month, &day);
    }

    if (nf != 3 || year < 0 || month < 0 || day < 0) {
	return -1;
    }

    if (year > 9999 || month > 12 || day > 31) {
	return -1;
    }

    if (year < 100) {
	year = FOUR_DIGIT_YEAR(year);
    }

    ret = (long)(year - 1) * 365 + leap_years_since_year_1(year - 1)
	+ day_in_year(day, month, year);

    return ret;
}

/**
 * calendar_obs_number:
 * @date: string representation of calendar date, in form
 * YY[YY]/MM/DD.
 * @dset: pointer to dataset information.
 * 
 * Returns: The zero-based observation number for the given
 * date within the current data set.
 */

int calendar_obs_number (const char *date, const DATASET *dset)
{
    long ed0 = (long) dset->sd0;
    long t = get_epoch_day(date);

    if (t == -1) return -1;
    
    /* subtract starting day for dataset */
    t -= ed0;

    if (dset->pd == 52) {
	/* weekly data */
	t /= 7;
    } else if (dset->pd == 5 || dset->pd == 6) { 
	/* daily, 5- or 6-day week: subtract number of irrelevant days */
	int startday = (((ed0 - 1 + SATURDAY) - NUMBER_MISSING_DAYS) % 7);
	int wkends = (t + startday - 1) / 7;

#ifdef CAL_DEBUG
	printf("calendar_obs_number: ed0=%d, date=%s, t=%d, startday=%d, wkends=%d\n", 
	       (int) ed0, date, (int) t, startday, wkends);
#endif

	if (dset->pd == 5) {
	    t -= (2 * wkends);
	} else {
	    t -= wkends;
	}
    }

    return (int) t;
}

static int t_to_epoch_day (int t, long start, int wkdays)
{
    int startday = (((start - 1 + SATURDAY) - NUMBER_MISSING_DAYS) % 7);
    int wkends = (t + startday - 1) / wkdays;

    if (wkdays == 5) wkends *= 2;

    return start + t + wkends;
}

/**
 * calendar_date_string:
 * @str: string to be filled out.
 * @t: zero-based index of observation.
 * @dset: pointer to dataset information.
 * 
 * Writes to @str the calendar representation of the date of
 * observation @t, in the form YY[YY]/MM/DD.
 *
 * Returns: 0 on success, non-zero on error.
 */

int calendar_date_string (char *str, int t, const DATASET *dset)
{
    int rem, yr;
    int add, day, mo = 0, modays = 0;
    long yrstart, dfind;
    long d0 = (long) dset->sd0;

    if (dset->pd == 52) {
	dfind = d0 + 7 * t;
    } else if (dset->pd == 7) {
	dfind = d0 + t;
    } else {
	/* 5- or 6-day data */
	if (t == 0 && dset->pd == 5) {
	    int wday = weekday_from_epoch_day(d0);

	    if (wday == 0 || wday == 6) {
		gretl_errmsg_sprintf(_("Invalid starting date for %d-day data"), dset->pd);
		*str = '\0';
		return E_DATA;
	    }
	}
	dfind = t_to_epoch_day(t, d0, dset->pd);
    }

    yr = 1 + (double) dfind / 365.248; 
    
    yrstart = (long)(yr - 1) * 365 + leap_years_since_year_1(yr - 1);
    rem = dfind - yrstart;

    if (rem <= 0) {
	yr--;
	yrstart = (long)(yr - 1) * 365 + leap_years_since_year_1(yr - 1);
	rem = dfind - yrstart;
    }

    while (modays < rem) {
	mo++;
	add = days_in_month[leap_year(yr)][mo];
	if (modays + add < rem) modays += add;
	else break;
    }

    day = rem - modays;

    if (strlen(dset->stobs) == 8) {
	sprintf(str, YMD_WRITE_Y2_FMT, yr % 100, mo, day);
    } else {
	sprintf(str, YMD_WRITE_Y4_FMT, yr, mo, day);
    }

    return 0;
}

/**
 * MS_excel_date_string:
 * @date: date string to be filled out.
 * @mst: MS Excel-type date code: days since base.
 * @pd: periodicity of data (or 0 if unknown).
 * @d1904: set to 1 if the base is 1904/01/01; otherwise
 * the base is assumed to be 1899/12/31.
 * 
 * Writes to @date the calendar representation of the date of
 * observation @mst, in the form YYYY-MM-DD if @pd is 0, 5,
 * 6, 7 or 52 (unknown, daily, or weekly frequency), otherwise 
 * in the appropriate format for annual, quarterly or monthly data.
 * 
 * Returns: 0.
 */

int MS_excel_date_string (char *date, int mst, int pd, int d1904)
{
    int yr = (d1904)? 1904 : 1900;
    int day = (d1904)? 2 : 1;
    int mo = 1;
    int leap, drem;

    if (mst == 0) {
	if (d1904) {
	    day = 1;
	} else {
	    yr = 1899;
	    mo = 12;
	    day = 31;
	}
    } else if (mst > 0) {
	drem = mst + d1904;

	while (1) {
	    int yd = 365 + leap_year(yr);

	    /* MS tomfoolery */
	    if (yr == 1900) yd++;

	    if (drem > yd) {
		drem -= yd;
		yr++;
	    } else {
		break;
	    }
	}

	leap = leap_year(yr) + (yr == 1900);

	for (mo=1; mo<13; mo++) {
	    int md = days_in_month[leap][mo];

	    if (drem > md) {
		drem -= md;
	    } else {
		day = drem;
		break;
	    }
	}
    } else {
	/* mst < 0, date prior to base */
	drem = - (mst + d1904);

	yr = (d1904)? 1903 : 1899;

	while (1) {
	    int yd = 365 + leap_year(yr);

	    if (drem > yd) {
		drem -= yd;
		yr--;
	    } else {
		break;
	    }
	}

	leap = leap_year(yr);

	for (mo=12; mo>0; mo--) {
	    int md = days_in_month[leap][mo];

	    if (drem >= md) {
		drem -= md;
	    } else {
		day = md - drem;
		break;
	    }
	}
    }

    if (pd == 1) {
	sprintf(date, "%d", yr);
    } else if (pd == 12) {
	sprintf(date, "%d:%02d", yr, mo);
    } else if (pd == 4) {
	int qtr = 1 + mo / 3.25;

	sprintf(date, "%d:%d", yr, qtr);
    } else {
	sprintf(date, YMD_WRITE_Y4_FMT, yr, mo, day);
    }

    return 0;
}

/**
 * get_dec_date:
 * @date: calendar representation of date: YYYY-MM-DD.
 * 
 * Returns: representation of date as year plus fraction of year.
 */

double get_dec_date (const char *date)
{
    char tmp[OBSLEN];
    int yr, mo, day, n;
    long ed0, edn, edt;
    double dyr, frac;

    n = sscanf(date, YMD_READ_FMT, &yr, &mo, &day);

    if (n != 3 && strchr(date, '/') != NULL) {
	/* backward compatibility */
	n = sscanf(date, "%d/%d/%d", &yr, &mo, &day);
    }

    if (n != 3) {
	return NADBL;
    }

    edt = get_epoch_day(date);
    sprintf(tmp, "%04d-01-01", yr);
    ed0 = get_epoch_day(tmp);
    sprintf(tmp, "%04d-12-31", yr);
    edn = get_epoch_day(tmp);

    if (yr < 100) {
	yr = FOUR_DIGIT_YEAR(yr);
    }

    dyr = yr;
    frac = ((double) edt - ed0) / ((double) edn - ed0 + 1.0);

    return dyr + frac;
}

/**
 * day_of_week:
 * @yr: year, preferably 4-digit.
 * @mo: month, 1 to 12.
 * @day: day in month, 1 to 31.
 * @err: location to receive error code.
 *
 * Returns: the day of the week for the supplied date
 * (Sunday = 0, Monday = 1, ...) or %NADBL on failure 
 * (the date is invalid).
 */

double day_of_week (int yr, int mo, int day, int *err)
{
    if (yr < 0 || mo <= 0 || day <= 0) {
	*err = E_DATA;
	return NADBL;
    }

    if (mo > 12 || day > 31) {
	*err = E_DATA;
	return NADBL;
    } 

    if (day > days_in_month[leap_year(yr)][mo]) {
	*err = E_DATA;
	return NADBL;
    }	

    if (yr < 100) {
	yr = FOUR_DIGIT_YEAR(yr);
    }

    return day_of_week_from_ymd(yr, mo, day);
}

#define day_in_calendar(w, d) ((w == 6 && d != 0) || \
                               (w == 5 && d != 0 && d != 6))

/**
 * weekday_from_date:
 * @date: calendar representation of date, [YY]YY/MM/DD
 * 
 * Returns: day of week as integer, Sunday = 0.
 */

int weekday_from_date (const char *date)
{
    int yr, mo, day;

    if (sscanf(date, YMD_READ_FMT, &yr, &mo, &day) != 3) {
	return -1;
    }

    if (yr < 100) {
	yr = FOUR_DIGIT_YEAR(yr);
    }

    return day_of_week_from_ymd(yr, mo, day);
}

/**
 * day_starts_month:
 * @d: day of month, 1-based
 * @m: month number, 1-based
 * @y: 4-digit year
 * @wkdays: number of days in week (7, 6 or 5)
 * @pad: location to receive 1 if the first day of the month
 * can reasonably be padded by a missing value (Jan 1), or NULL.
 * 
 * Returns: 1 if the day is the "first day of the month", 
 * allowance made for the possibility of a 5- or 6-day week, 
 * else 0.
 */

int day_starts_month (int d, int m, int y, int wkdays, int *pad)
{
    int ret = 0;

    if (wkdays == 7) {
	if (d == 1) {
	    ret = 1;
	} else if (pad != NULL && m == 1 && d == 2) {
	    /* second of January */
	    *pad = 1;
	    ret = 1;
	}
    } else {
	/* 5- or 6-day week: check for first weekday or non-Sunday */
	int i, wd;

	for (i=1; i<6; i++) {
	   wd = day_of_week_from_ymd(y, m, i); 
	   if (day_in_calendar(wkdays, wd)) break;
	}
	if (d == i) {
	    ret = 1;
	} else if (pad != NULL && m == 1 && d == i + 1) {
	    /* January again */
	    *pad = 1;
	    ret = 1;
	}
    }
    
    return ret;
}

/**
 * day_ends_month:
 * @d: day of month, 1-based
 * @m: month number, 1-based
 * @y: 4-digit year
 * @wkdays: number of days in week (7, 6 or 5)
 * 
 * Returns: 1 if the day is the "last day of the month", 
 * allowance made for the possibility of a 5- or 6-day week, else 0.
 */

int day_ends_month (int d, int m, int y, int wkdays)
{
    int ret = 0;
    int leap = (m == 2)? leap_year(y) : 0;
    int dm = days_in_month[leap][m];

    if (wkdays == 7) {
	ret = (d == dm);
    } else {
	/* 5- or 6-day week: check for last weekday or non-Sunday */
	int i, wd;

	for (i=dm; i>0; i--) {
	    wd = day_of_week_from_ymd(y, m, i);
	    if (day_in_calendar(wkdays, wd)) break;
	}
	ret = (d == i);	
    } 

    return ret;
}

/**
 * get_days_in_month:
 * @mon: month number, 1-based
 * @yr: 4-digit year
 * @wkdays: number of days in week (7, 6 or 5)
 * 
 * Returns: the number of days in the month, allowance made 
 * for the possibility of a 5- or 6-day week.
 */

int get_days_in_month (int mon, int yr, int wkdays)
{
    int ret = 0;
    int leap = (mon == 2)? leap_year(yr) : 0;
    int dm = days_in_month[leap][mon];

    if (wkdays == 7) {
	ret = dm;
    } else {
	int i, wd;

	for (i=0; i<dm; i++) {
	    wd = day_of_week_from_ymd(yr, mon, i + 1);
	    if (day_in_calendar(wkdays, wd)) ret++;
	}	
    } 

    return ret;
}

/**
 * days_in_month_before:
 * @yr: 4-digit year
 * @mon: month number, 1-based
 * @day: day in month.
 * @wkdays: number of days in week (7, 6 or 5)
 * 
 * Returns: the number of relevant days in the month prior to
 * the supplied date, allowing for the possibility of a 5- or 
 * 6-day week.
 */

int days_in_month_before (int yr, int mon, int day, int wkdays)
{
    int ret = 0;

    if (wkdays == 7) {
	ret = day - 1;
    } else {
	int i, wd;

	for (i=1; i<day; i++) {
	    wd = day_of_week_from_ymd(yr, mon, i);
	    if (day_in_calendar(wkdays, wd)) ret++;
	}	
    } 

    return ret;    
}

/**
 * days_in_month_after:
 * @yr: 4-digit year
 * @mon: month number, 1-based
 * @day: day in month.
 * @wkdays: number of days in week (7, 6 or 5)
 * 
 * Returns: the number of relevant days in the month after
 * the supplied date, allowing for the possibility of a 5- or 
 * 6-day week.
 */

int days_in_month_after (int yr, int mon, int day, int wkdays)
{
    int ret = 0;
    int leap = (mon == 2)? leap_year(yr) : 0;
    int dm = days_in_month[leap][mon];

    if (wkdays == 7) {
	ret = dm - day;
    } else {
	int i, wd;

	for (i=dm; i>day; i--) {
	    wd = day_of_week_from_ymd(yr, mon, i);
	    if (day_in_calendar(wkdays, wd)) ret++;
	}
    }

    return ret;    
}

/**
 * n_hidden_missing_obs:
 * @dset: dataset information.
 *
 * For daily data with user-supplied data strings, 
 * determine the number of "hidden" missing observations,
 * i.e. the difference between the actual number of
 * observations and the number that should be there,
 * according to the calendar. Allowance is made for 
 * 5- or 6-day data, via the data frequency given
 * in @dset.
 *
 * Returns: number of hidden observations.
 */

int n_hidden_missing_obs (const DATASET *dset)
{
    int t1, t2;
    int cal_n;

    if (!dated_daily_data(dset) || dset->S == NULL) {
	return 0;
    }
    
    t1 = calendar_obs_number(dset->S[0], dset);
    t2 = calendar_obs_number(dset->S[dset->n - 1], dset);

    cal_n = t2 - t1 + 1;

    return cal_n - dset->n;
}

/**
 * guess_daily_pd:
 * @dset: dataset information.
 *
 * Based on user-supplied daily date strings recorded in
 * @dset, try to guess whether the number of observations
 * per week is 5, 6 or 7 (given that some observations 
 * may be missing).
 *
 * Returns: best quess at data frequency.
 */

int guess_daily_pd (const DATASET *dset)
{
    int t, wd, pd = 5;
    int wdbak = -1;
    int havesat = 0;
    int gotsat = 0, gotsun = 0;
    int contig = 0;

    wd = weekday_from_date(dset->S[0]);
    if (6 - wd < dset->n) {
	havesat = 1;
    }

    for (t=0; t<dset->n && t<28; t++) {
	wd = weekday_from_date(dset->S[t]);
	if (wd == 0) {
	    gotsun = 1;
	} else if (wd == 6) {
	    gotsat = 1;
	}
	if ((wdbak + 1) % 7 == wd) {
	    contig++;
	}
	wdbak = wd;
    }

    if (gotsat && gotsun) {
	pd = 7;
    } else if (contig > 10) {
	if (gotsun) pd = 7;
	else if (gotsat) pd = 6;
    } else if (dset->n > 7) {
	if (!gotsun && !gotsat) {
	    pd = 5;
	} else if (!gotsun) {
	    pd = 6;
	}
    } else if (havesat && !gotsat) {
	pd = 5;
    } else {
	pd = 7;
    }

    return pd;
}

/**
 * iso_basic_to_extended:
 * @b: source array of YYYYMMDD values.
 * @y: array to hold year values.
 * @m: array to hold month values.
 * @d: array to hold day-of-week values.
 * @n: length of all the above arrays.
 *
 * Given the array @b of ISO 8601 "basic" daily dates (YYYYMMDD as
 * doubles), fill out the arrays @y, @m and @d with year, month
 * and day.
 *
 * Returns: 0.
 */

int iso_basic_to_extended (const double *b, double *y, double *m, double *d,
			   int n)
{
    int bi, yi, mi;
    int i;

    for (i=0; i<n; i++) {
	if (na(b[i])) {
	    y[i] = m[i] = NADBL;
	    if (d != NULL) {
		d[i] = NADBL;
	    }
	} else {
	    bi = (int) b[i];
	    yi = bi / 10000;
	    mi = (bi - 10000*yi) / 100;
	    y[i] = yi; m[i] = mi;
	    if (d != NULL) {
		d[i] = bi - 10000*yi - 100*mi;
	    }
	}
    }

    return 0;
}

/**
 * easterdate:
 * @year: year for which we want Easter date (Gregorian).
 *
 * Algorithm taken from Wikipedia page 
 * "https://en.wikipedia.org/wiki/Computus"
 * under the heading "Anonymous Gregorian algorithm".
 *
 * Returns the date of Easter in the Gregorian calendar as 
 * (month + day/100). Note that April the 10th, is, under
 * this convention, 4.1; hence, 4.2 is April the 20th, not
 * April the 2nd (which would be 4.02).
 */

double easterdate (int year)
{
    int a = year % 19;
    int b = year / 100;
    int c = year % 100;
    int d = b / 4;
    int e = b % 4;
    int f = (b + 8) / 25;
    int g = (b - f + 1) / 3;
    int h = (19 * a + b - d - g + 15) % 30;
    int i = c / 4;
    int k = c % 4;
    int L = (32 + 2 * e + 2 * i - h - k) % 7;
    int m = (a + 11 * h + 22 * L) / 451 ;

    int month = (h + L - 7 * m + 114) / 31;
    int day = ((h + L - 7 * m + 114) % 31) + 1;

    return month + day * 0.01;
}

