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

#include "libgretl.h"

#ifdef WIN32
# include "gretl_win32.h" /* for strptime() */
#endif

#define CAL_DEBUG 0

/**
 * SECTION:calendar
 * @short_description: functions for working with dates
 * @title: Calendar
 * @include: libgretl.h
 *
 * Here we have various functions dealing with calendar dates;
 * for the most part these are designed for handling daily
 * time-series data.
 */

static int days_in_month[2][13] = {
    {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
};

/* Note: @w is the number of days in a week, and @d is the
   day-of-week number for a given day. This macro assumes
   that day-of-week runs from Monday = 1 to Sunday = 7.
*/
#define day_in_calendar(w, d) ((d) <= (w))

/* Jan 1, 0001, was a Monday on the proleptic Gregorian calendar */
#define DAY1 1

/* The number that must be added to an "epoch day" (see below)
   to obtain the astronomers' Julian Day Number
*/
#define JDN_ADJ 1721425

/* Note on the GLib API used below: where "julian" occurs in the names
   of GLib calendrical functions it refers to the "Julian day"; that is,
   the number of days since some fixed starting point, as used by
   astronomers. This is quite distinct from the Julian calendar.
   However, GLib takes Julian day 1 to be the first of January in
   AD 1, as opposed to the astronomical starting point in 4713 BC,
   so these are not strictly Julian days, and in our own functions
   which use the same concept we call them "epoch days".
*/

static int leap_year (int yr)
{
    return (!(yr % 4) && (yr % 100)) || !(yr % 400);
}

static int valid_ymd (int y, int m, int d, int julian)
{
    int ok = 1;

    if (!g_date_valid_dmy(d, m, y)) {
	ok = 0;
	if (julian && y > 0 && m == 2 && d == 29 && y%4 == 0) {
	    ok = 1;
	}
    }

    return ok;
}

/* Returns day of week on the GLib numbering, Monday=1 to Sunday=7,
   or 0 in case an invalid date was given.
*/

static int day_of_week_from_ymd (int y, int m, int d, int julian)
{
    GDate date;

    if (!valid_ymd(y, m, d, julian)) {
	return 0;
    }

    g_date_clear(&date, 1);

    if (julian) {
	guint32 ed = epoch_day_from_julian_ymd(y, m, d);

	g_date_set_julian(&date, ed);
    } else {
	g_date_set_dmy(&date, d, m, y);
    }

    return g_date_get_weekday(&date);
}

/**
 * epoch_day_from_ymd:
 * @y: year (1 <= y <= 9999).
 * @m: month (1 <= m <= 12).
 * @d: day of month (1 <= d <= 31).
 *
 * Returns: the epoch day number, which equals 1 for the first of
 * January in the year AD 1 on the proleptic Gregorian calendar,
 * or 0 on error.
 */

guint32 epoch_day_from_ymd (int y, int m, int d)
{
    GDate date;

    if (!g_date_valid_dmy(d, m, y)) {
	return 0;
    }

    g_date_clear(&date, 1);
    g_date_set_dmy(&date, d, m, y);

    return g_date_get_julian(&date);
}

/**
 * nearby_epoch_day:
 * @y: year (1 <= y <= 9999).
 * @m: month (1 <= m <= 12).
 * @d: day of month (1 <= d <= 31).
 * @wkdays: number of days per week.
 *
 * Returns: 0 if the supplied date is invalid, otherwise the
 * epoch day number for @y, @m, @d, or the first subsequent
 * epoch day if the date in question is not included in the
 * calendar, based on @wkdays (for example, the day is a
 * Sunday and @wkdays is 5 or 6).
 */

guint32 nearby_epoch_day (int y, int m, int d, int wkdays)
{
    GDate date;
    int dow;
    guint32 j;

    if (!g_date_valid_dmy(d, m, y)) {
	return 0;
    }

    g_date_clear(&date, 1);
    g_date_set_dmy(&date, d, m, y);
    j = g_date_get_julian(&date);
    dow = g_date_get_weekday(&date);

    if (wkdays != 7 && dow == G_DATE_SUNDAY) {
	j++;
    } else if (wkdays == 5 && dow == G_DATE_SATURDAY) {
	j += 2;
    }

    return j;
}

/**
 * ymd_bits_from_epoch_day:
 * @ed: epoch day (ed >= 1).
 * @y: location to receive year.
 * @m: location to receive month.
 * @m: location to receive day.
 *
 * Fills @y, @m and @d with year, month and day on the (possibly
 * proleptic) Gregorian calendar.
 *
 * Returns: 0 on success, non-zero on error.
 */

int ymd_bits_from_epoch_day (guint32 ed, int *y, int *m, int *d)
{
    GDate date;

    if (!g_date_valid_julian(ed)) {
	return E_INVARG;
    }

    g_date_clear(&date, 1);
    g_date_set_julian(&date, ed);

    *y = g_date_get_year(&date);
    *m = g_date_get_month(&date);
    *d = g_date_get_day(&date);

    return 0;
}

/**
 * julian_ymd_bits_from_epoch_day:
 * @ed: epoch day (ed >= 1).
 * @y: location to receive year.
 * @m: location to receive month.
 * @m: location to receive day.
 *
 * Given an epoch day, fills @y, @m and @d with year, month and day on
 * the (possibly proleptic) Julian calendar.
 *
 * Follows the algorithm of E.G. Richards (2013), "Calendars," In S.E.
 * Urban & P.K. Seidelmann, eds. Explanatory Supplement to the Astronomical
 * Almanac, 3rd ed. (pp. 585-624), Mill Valley, CA: University Science Books
 * (as set out on https://en.wikipedia.org/wiki/Julian_day).
 *
 * There are other algorithms for this purpose on the internet but they are
 * mostly wrong (at least, not right for all dates); many of them fail
 * the round-trip test (date -> epoch day -> date) for some dates.
 *
 * Returns: 0 on success, non-zero on error.
 */

int julian_ymd_bits_from_epoch_day (guint32 ed, int *y,
				    int *m, int *d)
{
    int x = 4716;
    int p = 1461;
    int f = ed + JDN_ADJ + 1401;
    int e = 4 * f + 3;
    int g = (e % p)/4;
    int h = 5 * g + 2;

    /* The addition of JDN_ADJ above translates from our
       "epoch day" to Julian Day Number; the addition of 1401
       to the JDN is specified by Richards' algorithm.
    */

    *d = (h % 153)/5 + 1;
    *m = (h/153 + 2) % 12 + 1;
    *y = e/p - x + (14 - *m)/12;

    return 0;
}

/**
 * epoch_day_from_julian_ymd:
 * @y: year (y >= 1).
 * @m: month (1 <= m <= 12).
 * @d: day of month (1 <= d <= 31).
 *
 * The @y, @m and @d arguments are assumed to refer to a date on
 * the Julian calendar. The conversion algorithm is taken from
 * https://en.wikipedia.org/wiki/Julian_day, where it appears to
 * be credited to the Department of Computer Science at UT, Austin.
 *
 * Returns: the epoch day number, which equals 1 for the first of
 * January in the year AD 1 on the proleptic Gregorian calendar,
 * or 0 on error.
 */

guint32 epoch_day_from_julian_ymd (int y, int m, int d)
{
    int a = (14 - m)/12;
    int jd;

    y = y + 4800 - a;
    m = m + 12*a - 3;

    jd = d + (153*m + 2)/5 + 365*y + y/4 - 32083;

    if (jd <= JDN_ADJ) {
	/* prior to AD 1 */
	return 0;
    } else {
	return (guint32) jd - JDN_ADJ;
    }
}

/**
 * ymd_extended_from_epoch_day:
 * @ed: epoch day (ed >= 1).
 * @julian: non-zero to use Julian calendar, otherwise Gregorian.
 * @err: location to receive error code.
 *
 * Returns: a string on the pattern YYYY-MM-DD (ISO 8601 extended
 * date format) given the epoch day number, which equals 1 for the
 * first of January in the year 1 AD, or NULL on error.
 */

char *ymd_extended_from_epoch_day (guint32 ed, int julian, int *err)
{
    char *ret = NULL;
    int y, m, d;
    int myerr;

    if (julian) {
	myerr = julian_ymd_bits_from_epoch_day(ed, &y, &m, &d);
    } else {
	myerr = ymd_bits_from_epoch_day(ed, &y, &m, &d);
    }

    if (!myerr) {
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
 * @julian: non-zero to use Julian calendar, otherwise Gregorian.
 * @err: location to receive error code.
 *
 * Returns: an 8-digit number on the pattern YYYYMMDD (ISO 8601 basic
 * date format) given the epoch day number, which equals 1 for the
 * first of January in the year 1 AD, or #NADBL on error.
 */

double ymd_basic_from_epoch_day (guint32 ed, int julian, int *err)
{
    int y = 0, m = 0, d = 0;

    if (julian) {
	*err = julian_ymd_bits_from_epoch_day(ed, &y, &m, &d);
    } else {
	*err = ymd_bits_from_epoch_day(ed, &y, &m, &d);
    }

    if (*err) {
	return NADBL;
    } else {
	return 10000*y + 100*m + d;
    }
}

/**
 * epoch_day_from_ymd_basic:
 * @ymd: number that is supposed to represent YYYYMMDD.
 *
 * Returns: the epoch day number corresponding to @ymd, interpreted
 * as YYYYMMDD (ISO 8601 "basic") if possible, or 0 on error.
 */

guint32 epoch_day_from_ymd_basic (double ymd)
{
    char bit[4] = {0};
    char tmp[9];
    int y, m, d, n;

    if (ymd != floor(ymd) || ymd < 10101 || ymd > 99991231) {
	return 0;
    }

    sprintf(tmp, "%.0f", ymd);
    n = strlen(tmp);

    /* get day from the last two digits */
    bit[0] = tmp[n-2];
    bit[1] = tmp[n-1];
    d = atoi(bit);

    /* get month from the preceding two digits */
    bit[0] = tmp[n-4];
    bit[1] = tmp[n-3];
    m = atoi(bit);

    /* get year from the leading digit(s) */
    tmp[n-4] = '\0';
    y = atoi(tmp);

    return epoch_day_from_ymd(y, m, d);
}

/**
 * weekday_from_epoch_day:
 * @ed: epoch day (ed >= 1).
 *
 * Returns: the weekday corrsponding to @ed, on the basis
 * of Monday = 1 to Sunday = 7, or 0 on error.
 */

int weekday_from_epoch_day (guint32 ed)
{
    GDate date;

    if (!g_date_valid_julian(ed)) {
	return 0;
    }

    g_date_clear(&date, 1);
    g_date_set_julian(&date, ed);

    return g_date_get_weekday(&date);
}

/**
 * get_epoch_day:
 * @datestr: string representation of calendar date, in form
 * YY[YY]-MM-DD.
 *
 * Returns: the epoch day number, or 0 on failure.
 */

guint32 get_epoch_day (const char *datestr)
{
    GDate date;
    int y, m, d, nf = 0;
    int ydigits = 0;

    if (strchr(datestr, '-')) {
	ydigits = strcspn(datestr, "-");
	nf = sscanf(datestr, YMD_READ_FMT, &y, &m, &d);
    } else if (strchr(datestr, '/')) {
	ydigits = strcspn(datestr, "/");
	nf = sscanf(datestr, "%d/%d/%d", &y, &m, &d);
    } else if (strlen(datestr) == 8) {
	ydigits = 4;
	nf = sscanf(datestr, "%4d%2d%2d", &y, &m, &d);
    }

    if (nf != 3 || (ydigits != 4 && ydigits != 2)) {
	return 0;
    }

    if (ydigits == 2) {
	y = FOUR_DIGIT_YEAR(y);
    }

    if (!g_date_valid_dmy(d, m, y)) {
	return 0;
    }

    g_date_clear(&date, 1);
    g_date_set_dmy(&date, d, m, y);

    return g_date_get_julian(&date);
}

/* Note that ed0 cannot be a Sunday, since we're working with a 5- or
   6-day calendar here. So the mod-7 operation will give us the
   correct day-of-week.
*/

static int subtract_irrelevant_days (guint32 ed0, int wkdays, int t)
{
    int dow0 = (ed0 - 1 + DAY1) % 7;
    int wkends = (t + dow0 - 1) / 7;

    if (wkdays == 5) {
	t -= 2 * wkends;
    } else {
	t -= wkends;
    }

    return t;
}

/**
 * calendar_obs_number:
 * @datestr: string representation of calendar date, in form
 * YY[YY]-MM-DD.
 * @dset: pointer to dataset.
 * @nolimit: allow @datestr to be outside of the range of
 * @dset.
 *
 * Returns: The zero-based observation number for the given
 * date within the current dataset.
 */

int calendar_obs_number (const char *datestr, const DATASET *dset,
			 int nolimit)
{
    guint32 ed0 = (guint32) dset->sd0;
    guint32 ut = get_epoch_day(datestr);
    int t = (int) ut;

#if CAL_DEBUG
    fprintf(stderr, "calendar_obs_number: '%s' -> epoch day %u, dow %d\n",
	    datestr, ut, weekday_from_epoch_day(ut));
#endif

    if (t <= 0 || (!nolimit && t < ed0)) {
	return -1;
    } else if (t == ed0) {
	return 0;
    }

    /* subtract starting day for dataset */
    t -= ed0;

    if (!nolimit && t <= 0) {
	return -1;
    }

    if (dset->pd == 52) {
	/* weekly data */
	t /= 7;
    } else if (dset->pd == 5 || dset->pd == 6) {
	/* check for out-of-calendar input day */
	if (weekday_from_epoch_day(ut) > dset->pd) {
	    fprintf(stderr, "Invalid date %s for %d-day data\n",
		    datestr, dset->pd);
	    return -1;
	} else {
	    t = subtract_irrelevant_days(ed0, dset->pd, t);
	}
    }

    return t;
}

/* Convert from 0-based @t in dataset to epoch day, assuming a
   complete 5- or 6-day calendar (specified by @wkdays).
*/

static guint32 t_to_epoch_day (int t, guint32 ed0, int wkdays)
{
    int startdow = (ed0 - 1 + DAY1) % 7;
    int wkends = (t + startdow - 1) / wkdays;

    if (wkdays == 5) {
	wkends *= 2;
    }

    return ed0 + t + wkends;
}

/**
 * epoch_day_from_t:
 * @t: 0-based observation index.
 * @dset: pointer to dataset.
 *
 * Returns: the epoch day based on calendrical information
 * in @dset. In case of error 0 is returned.
 */

guint32 epoch_day_from_t (int t, const DATASET *dset)
{
    guint32 ed0 = (guint32) dset->sd0;
    guint32 dt = 0;

    if (t == 0) {
	return ed0;
    } else if (dset->pd == 52) {
	dt = ed0 + 7 * t;
    } else if (dset->pd == 7) {
	dt = ed0 + t;
    } else {
	/* 5- or 6-day daily data */
	dt = t_to_epoch_day(t, ed0, dset->pd);
    }

    return dt;
}

/**
 * calendar_date_string:
 * @targ: string to be filled out.
 * @t: zero-based index of observation.
 * @dset: pointer to dataset.
 *
 * Writes to @targ the calendar representation of the date of
 * observation @t, in the form YY[YY]-MM-DD, on the assumption
 * that @dset contains calendrical information.
 *
 * Returns: 0 on success, non-zero on error.
 */

int calendar_date_string (char *targ, int t, const DATASET *dset)
{
    guint32 d0 = (guint32) dset->sd0;
    guint32 dt = 0;
    int y, m, d;
    int err = 0;

    if (d0 == 1) {
        strcpy(targ, "BAD DATE");
        gretl_errmsg_set("The dataset lacks calendrical information");
        return E_DATA;
    } else if (dataset_has_markers(dset)) {
        /* the calendar is presumably incomplete */
        strcpy(targ, dset->S[t]);
        return 0;
    }

    *targ = '\0';

    if (dset->pd == 52) {
	dt = d0 + 7 * t;
    } else if (dset->pd == 7) {
	dt = d0 + t;
    } else {
	if (t == 0 && (dset->pd == 5 || dset->pd == 6)) {
	    int wd = weekday_from_epoch_day(d0);

	    if (wd == 0 || wd > dset->pd) {
		strcpy(targ, "BAD DATE");
		gretl_errmsg_sprintf(_("Invalid starting date for %d-day data"), dset->pd);
		return E_DATA;
	    }
	}
	dt = t_to_epoch_day(t, d0, dset->pd);
    }

    err = ymd_bits_from_epoch_day(dt, &y, &m, &d);

    if (!err) {
	if (strlen(dset->stobs) == 8) {
	    sprintf(targ, YMD_WRITE_Y2_FMT, y % 100, m, d);
	} else {
	    sprintf(targ, YMD_WRITE_Y4_FMT, y, m, d);
	}
    }

    return err;
}

/**
 * MS_excel_date_string:
 * @targ: date string to be filled out.
 * @mst: MS Excel-type date code: days since base.
 * @pd: periodicity of data (or 0 if unknown).
 * @d1904: set to 1 if the base is 1904/01/01; otherwise
 * the base is assumed to be 1899/12/31.
 *
 * Writes to @targ the calendar representation of the date of
 * observation @mst, in the form YYYY-MM-DD if @pd is 0, 5,
 * 6, 7 or 52 (unknown, daily, or weekly frequency), otherwise
 * in the appropriate format for annual, quarterly or monthly
 * according to @pd.
 *
 * Returns: 0.
 */

int MS_excel_date_string (char *targ, int mst, int pd, int d1904)
{
    int y = (d1904)? 1904 : 1900;
    int d = (d1904)? 2 : 1;
    int m = 1;
    int leap, drem;

    *targ = '\0';

    if (mst == 0) {
	/* date coincident with base */
	if (d1904) {
	    d = 1;
	} else {
	    y = 1899;
	    m = 12;
	    d = 31;
	}
    } else if (mst > 0) {
	/* date subsequent to base */
	drem = mst + d1904;

	while (1) {
	    int yd = 365 + leap_year(y);

	    /* MS nincompoopery */
	    if (y == 1900) yd++;

	    if (drem > yd) {
		drem -= yd;
		y++;
	    } else {
		break;
	    }
	}

	leap = leap_year(y) + (y == 1900);

	for (m=1; m<=12; m++) {
	    int md = days_in_month[leap][m];

	    if (drem > md) {
		drem -= md;
	    } else {
		d = drem;
		break;
	    }
	}
    } else {
	/* mst < 0, date prior to base */
	drem = -(mst + d1904);

	y = (d1904)? 1903 : 1899;

	while (1) {
	    int yd = 365 + leap_year(y);

	    if (drem > yd) {
		drem -= yd;
		y--;
	    } else {
		break;
	    }
	}

	leap = leap_year(y);

	for (m=12; m>0; m--) {
	    int md = days_in_month[leap][m];

	    if (drem >= md) {
		drem -= md;
	    } else {
		d = md - drem;
		break;
	    }
	}
    }

    if (pd == 1) {
	sprintf(targ, "%d", y);
    } else if (pd == 12) {
	sprintf(targ, "%d:%02d", y, m);
    } else if (pd == 4) {
	int q = 1 + m / 3.25;

	sprintf(targ, "%d:%d", y, q);
    } else {
	sprintf(targ, YMD_WRITE_Y4_FMT, y, m, d);
    }

    return 0;
}

/**
 * get_dec_date:
 * @datestr: calendar representation of date: YYYY-MM-DD.
 *
 * Returns: representation of date as year plus fraction of year.
 */

double get_dec_date (const char *datestr)
{
    GDate date;
    int y, m, d, nf;
    int ydigits = 0;
    double yrday, frac;

    nf = sscanf(datestr, YMD_READ_FMT, &y, &m, &d);

    if (nf == 3) {
	ydigits = strcspn(datestr, "-");
    } else if (strchr(datestr, '/') != NULL) {
	/* backward compatibility */
	ydigits = strcspn(datestr, "/");
	nf = sscanf(datestr, "%d/%d/%d", &y, &m, &d);
    }

    if (nf != 3 || (ydigits != 4 && ydigits != 2)) {
	return NADBL;
    }

    if (ydigits == 2) {
	y = FOUR_DIGIT_YEAR(y);
    }

    if (!g_date_valid_dmy(d, m, y)) {
	return NADBL;
    }

    g_date_clear(&date, 1);
    g_date_set_dmy(&date, d, m, y);
    yrday = g_date_get_day_of_year(&date);
    frac = yrday / (365 + g_date_is_leap_year(y));

    return y + frac;
}

/**
 * legacy_day_of_week:
 * @y: year.
 * @m: month, 1 to 12.
 * @d: day in month, 1 to 31.
 * @julian: non-zero to use Julian calendar, otherwise Gregorian.
 * @err: location to receive error code.
 *
 * Returns: the day of the week for the supplied date on the
 * "lagacy" numbering of Sunday = 0 to Saturday = 6) or %NADBL
 * on failure (that is, the supplied date is invalid).
 */

double legacy_day_of_week (int y, int m, int d, int julian, int *err)
{
    int wd = day_of_week_from_ymd(y, m, d, julian);

    if (wd == 0) {
	return NADBL;
    } else {
	return wd == G_DATE_SUNDAY ? 0 : wd;
    }
}

/**
 * weekday_from_date:
 * @datestr: calendar representation of date, [YY]YY/MM/DD
 *
 * Returns: day of week as integer (Monday = 1 to Sunday = 7).
 */

int weekday_from_date (const char *datestr)
{
    int y, m, d;
    int ydigits;

    if (sscanf(datestr, YMD_READ_FMT, &y, &m, &d) != 3) {
	return -1;
    }

    ydigits = strcspn(datestr, "-");
    if (ydigits != 4 && ydigits != 2) {
	return -1;
    }

    if (ydigits == 2) {
	y = FOUR_DIGIT_YEAR(y);
    }

    return day_of_week_from_ymd(y, m, d, 0);
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
	int i, idx = day_of_week_from_ymd(y, m, 1, 0);

	for (i=1; i<6; i++) {
	   if (day_in_calendar(wkdays, idx)) {
	       break;
	   }
	   idx = idx == 7 ? 1 : idx+1;
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
	int i, idx = day_of_week_from_ymd(y, m, dm, 0);

	for (i=dm; i>0; i--) {
	    if (day_in_calendar(wkdays, idx)) {
		break;
	    }
	    idx = idx == 1 ? 7 : idx-1;
	}
	ret = (d == i);
    }

    return ret;
}

/**
 * get_days_in_month:
 * @m: month number, 1-based
 * @y: 4-digit year
 * @wkdays: number of days in week (7, 6 or 5)
 * @julian: non-zero for Julian calendar, otherwise Gregorian.
 *
 * Returns: the number of (relevant) days in the month, allowance
 * made for the possibility of a 5- or 6-day week.
 */

int get_days_in_month (int m, int y, int wkdays, int julian)
{
    int dm, leap = 0;
    int ret = 0;

    if (m == 2) {
	leap = julian ? (y%4 == 0) : leap_year(y);
    }
    dm = days_in_month[leap][m];

    if (wkdays == 7) {
	ret = dm;
    } else {
	int i, idx = day_of_week_from_ymd(y, m, 1, julian);

	for (i=0; i<dm; i++) {
	    if (day_in_calendar(wkdays, idx)) {
		ret++;
	    }
	    idx = idx == 7 ? 1 : idx+1;
	}
    }

    return ret;
}

/**
 * fill_monthlen_array:
 * @mlen: array to be filled.
 * @t1: starting index for fill.
 * @t2: stopping index for fill.
 * @wkdays: number of days in week (7, 6 or 5).
 * @mo: month (ignored if @movec is non-NULL).
 * @yr: year (ignored if @yrvec is non-NULL).
 * @movec: array of month values (or NULL).
 * @yrvec: array of year values (or NULL).
 * @julian: non-zero for Julian calendar, otherwise Gregorian.
 *
 * Fills @mlen from @t1 to @t2 with the number of days in
 * each month/year pair. It is assumed that at least one
 * of @movec and @yrvec is non-NULL. The various arrays
 * are doubles only because in context they will be
 * dataset series; in concept they are integers.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int fill_monthlen_array (double *mlen, int t1, int t2,
			 int wkdays, int mo, int yr,
			 const double *movec,
			 const double *yrvec,
			 int julian)
{
    int t, err = 0;

    if (movec == NULL && yrvec == NULL) {
	return E_INVARG;
    }

    for (t=t1; t<=t2; t++) {
	if (movec != NULL) {
	    mo = gretl_int_from_double(movec[t], &err);
	    if (err || mo < 1 || mo > 12) {
		err = E_INVARG;
		break;
	    }
	}
	if (yrvec != NULL) {
	    yr = gretl_int_from_double(yrvec[t], &err);
	    if (err) {
		break;
	    }
	    if (yr < 0) {
		yr = -yr;
		julian = 1;
	    } else {
		julian = 0;
	    }
	}
	mlen[t] = get_days_in_month(mo, yr, wkdays, julian);
    }

    return err;
}

/**
 * month_day_index:
 * @y: 4-digit year
 * @m: month number, 1-based
 * @d: day in month.
 * @wkdays: number of days in week (7, 6 or 5)
 *
 * Returns: the 1-based index of calendar day @d in month
 * @m of year @y, allowing for the possibility of a 5- or
 * 6-day week.
 */

int month_day_index (int y, int m, int d, int wkdays)
{
    int ret = 0;

    if (wkdays == 7) {
	ret = d;
    } else {
	int i, dow = day_of_week_from_ymd(y, m, 1, 0);

	for (i=1; i<=d; i++) {
	    if (day_in_calendar(wkdays, dow)) {
		ret++;
	    }
	    dow = dow == 7 ? 1 : dow+1;
	}
    }

    return ret;
}

/**
 * days_in_month_before:
 * @y: 4-digit year
 * @m: month number, 1-based
 * @d: day in month.
 * @wkdays: number of days in week (7, 6 or 5)
 *
 * Returns: the number of relevant days in the month prior to
 * the supplied date, allowing for the possibility of a 5- or
 * 6-day week.
 */

int days_in_month_before (int y, int m, int d, int wkdays)
{
    int ret = 0;

    if (wkdays == 7) {
	ret = d - 1;
    } else {
	int i, dow = day_of_week_from_ymd(y, m, 1, 0);

	for (i=1; i<d; i++) {
	    if (day_in_calendar(wkdays, dow)) {
		ret++;
	    }
	    dow = dow == 7 ? 1 : dow+1;
	}
    }

    return ret;
}

/**
 * days_in_month_after:
 * @y: 4-digit year
 * @m: month number, 1-based
 * @d: day in month.
 * @wkdays: number of days in week (7, 6 or 5)
 *
 * Returns: the number of relevant days in the month after
 * the supplied date, allowing for the possibility of a 5- or
 * 6-day week.
 */

int days_in_month_after (int y, int m, int d, int wkdays)
{
    int leap = (m == 2)? leap_year(y) : 0;
    int dm = days_in_month[leap][m];
    int ret = 0;

    if (wkdays == 7) {
	ret = dm - d;
    } else {
	int i, dow = day_of_week_from_ymd(y, m, dm, 0);

	for (i=dm; i>d; i--) {
	    if (day_in_calendar(wkdays, dow)) {
		ret++;
	    }
	    dow = dow == 7 ? 1 : dow-1;
	}
    }

    return ret;
}

/**
 * date_to_daily_index:
 * @datestr: date in format YYYY-MM-DD.
 * @wkdays: number of days in week (7, 6 or 5)
 *
 * Returns: the zero-based index of the specified day
 * within the specified month and year. In the case
 * of 5- or 6-day data index zero does not necessarily
 * correspond to the first day of the month but rather
 * to the first relevant day.
 */

int date_to_daily_index (const char *datestr, int wkdays)
{
    int y, m, d, seq = 0;

    if (sscanf(datestr, YMD_READ_FMT, &y, &m, &d) != 3) {
	return -1;
    }

    if (wkdays == 7) {
	seq = d - 1;
    } else {
	int leap = leap_year(y);
	int n = days_in_month[leap][m];
	int i, dow = day_of_week_from_ymd(y, m, 1, 0);

	for (i=1; i<=n; i++) {
	    if (d == i) {
		break;
	    }
	    if (day_in_calendar(wkdays, dow)) {
		seq++;
	    }
	    dow = dow == 7 ? 1 : dow+1;
	}
    }

    return seq;
}

/**
 * daily_index_to_date:
 * @targ: location to receive the date (YYYY-MM-DD).
 * @y: year.
 * @m: month.
 * @idx: zero-based index of day within month.
 * @wkdays: number of days in week (7, 6 or 5)
 *
 * Fills out @targ with the calendar data implied by
 * the specification of @y, @m, @seq and @wkdays,
 * provided this specification corresponds to an actual
 * calendar date.

 * Returns: 0 on successful completion, non-zero if
 * there is no such calendar date.
 */

int daily_index_to_date (char *targ, int y, int m, int idx,
			 int wkdays)
{
    int day = 0;

    *targ = '\0';

    if (m < 1 || m > 12 || idx < 0 || idx > 30) {
	fprintf(stderr, "daily_index_to_date: y=%d, m=%d, seq=%d\n",
		y, m, idx);
	return E_INVARG;
    }

    if (wkdays == 7) {
	day = idx + 1;
    } else {
	int leap = leap_year(y);
	int n = days_in_month[leap][m];
	int dow = day_of_week_from_ymd(y, m, 1, 0);
	int i, seq = 0;

	for (i=1; i<=n; i++) {
	    if (day_in_calendar(wkdays, dow)) {
		if (seq == idx) {
		    day = i;
		    break;
		}
		seq++;
	    }
	    dow = dow == 7 ? 1 : dow+1;
	}
    }

    if (day <= 0) {
	return E_DATA;
    } else {
	sprintf(targ, YMD_WRITE_FMT, y, m, day);
	return 0;
    }
}

/**
 * n_hidden_missing_obs:
 * @dset: dataset information.
 * @t1: first observation.
 * @t2: last observation.
 *
 * For daily data with user-supplied data strings,
 * determine the number of "hidden" missing observations
 * in the range @t1 to @t2 inclusive. This is
 * the difference between the actual number of
 * observations and the number that should be there,
 * according to the calendar. Allowance is made for
 * 5- or 6-day data, via the data frequency given
 * in @dset.
 *
 * Returns: number of hidden observations.
 */

int n_hidden_missing_obs (const DATASET *dset, int t1, int t2)
{
    int n_present = t2 - t1 + 1;
    int cal_n;

    if (!dated_daily_data(dset) || dset->S == NULL) {
	return 0;
    }

    t1 = calendar_obs_number(dset->S[t1], dset, 0);
    t2 = calendar_obs_number(dset->S[t2], dset, 0);

    cal_n = t2 - t1 + 1;

    return cal_n - n_present;
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
    int gotsat = 0;

    for (t=0; t<dset->n && t<40; t++) {
	wd = weekday_from_date(dset->S[t]);
	if (wd == 7) {
	    pd = 7;
	    break;
	} else if (wd == 6) {
	    gotsat = 1;
	}
    }

    if (pd == 5 && gotsat) {
	pd = 6;
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

int iso_basic_to_extended (const double *b, double *y, double *m,
			   double *d, int n)
{
    int bi, yi, mi, di;
    int i, julian;

    for (i=0; i<n; i++) {
	julian = 0;
	if (na(b[i])) {
	    y[i] = m[i] = NADBL;
	    if (d != NULL) {
		d[i] = NADBL;
	    }
	} else {
	    bi = (int) b[i];
	    if (bi < 0) {
		julian = 1;
		bi = -bi;
	    }
	    yi = bi / 10000;
	    mi = (bi - 10000*yi) / 100;
	    di = bi - 10000*yi - 100*mi;
	    /* now check for legit date */
	    if (!valid_ymd(yi, mi, di, julian)) {
		y[i] = m[i] = NADBL;
		if (d != NULL) {
		    d[i] = NADBL;
		}
	    } else {
		y[i] = yi;
		m[i] = mi;
		if (d != NULL) {
		    d[i] = di;
		}
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
 * https://en.wikipedia.org/wiki/Computus
 * under the heading "Anonymous Gregorian algorithm".
 *
 * Returns: the date of Easter in the Gregorian calendar as
 * (month + day/100). Note that April the 10th is, under
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

/**
 * dayspan:
 * @ed1: first epoch day.
 * @ed2: last epoch day.
 * @wkdays: relevant days per week (5, 6 or 7).
 * @err: location to receive error code.
 *
 * Returns: The number of days in the interval @ed1 to
 * @ed2, inclusive, taking account of the number of daily
 * observations per week, @wkdays. If @wkdays = 6 Sundays
 * are disregarded; if @wkdays = 5 both Saturdays and
 * Sundays are disregarded.
 */

int day_span (guint32 ed1, guint32 ed2, int wkdays, int *err)
{
    int n = 0;

    if (ed2 < ed1 || !g_date_valid_julian(ed1) ||
	!g_date_valid_julian(ed2)) {
	if (err != NULL) {
	    *err = E_INVARG;
	}
    } else if (wkdays == 7) {
	/* simple! */
	n = ed2 - ed1 + 1;
    } else {
	GDate date;
	guint32 i;
	int wd;

	g_date_clear(&date, 1);
	g_date_set_julian(&date, ed1);
	wd = g_date_get_weekday(&date);

	for (i=ed1; i<=ed2; i++) {
	    if (day_in_calendar(wkdays, wd)) {
		n++;
	    }
	    wd = wd == 7 ? 1 : wd+1;
	}
    }

    return n;
}

/**
 * iso_week_number:
 * @y: year.
 * @m: month (1 to 12).
 * @d: day on month.
 * @err: location to receive error code.
 *
 * Returns: The ISO 8601 week number (1 to 53) corresponding
 * to the Gregorian date specified by the first 3 arguments,
 * or -1 on invalid input.
 */

int iso_week_number (int y, int m, int d, int *err)
{
    GDate date;
    int wnum = -1;

    if (!g_date_valid_dmy(d, m, y)) {
	*err = E_INVARG;
    } else {
	g_date_clear(&date, 1);
	g_date_set_dmy(&date, d, m, y);
	wnum = g_date_get_iso8601_week_of_year(&date);
    }

    return wnum;
}

/**
 * iso_week_from_date:
 * @datestr: date in ISO 8601 format, YYYY-MM-DD.
 *
 * Returns: The ISO 8601 week number (1 to 53) corresponding
 * to the Gregorian date specified by @datestr, or -1 on
 * invalid input.
 */

int iso_week_from_date (const char *datestr)
{
    int y, m, d;
    int ydigits;
    int err = 0;

    if (sscanf(datestr, YMD_READ_FMT, &y, &m, &d) != 3) {
	return -1;
    }

    ydigits = strcspn(datestr, "-");

    if (ydigits != 4 && ydigits != 2) {
	return -1;
    } else if (ydigits == 2) {
	y = FOUR_DIGIT_YEAR(y);
    }

    return iso_week_number(y, m, d, &err);
}

/**
 * gretl_strfdate:
 * @s: target string.
 * @slen: length of target string.
 * @format: as per strftime().
 * @ed: days since 1 CE.
 *
 * If @ed is found to be valid, writes a string representing
 * the date of @ed to @s, governed by @format.
 *
 * Returns: The number of characters written to @s, or 0 in case
 * of invalid input.
 */

int gretl_strfdate (char *s, int slen, const char *format,
		    guint32 ed)
{
    int ret = 0;

    if (g_date_valid_julian(ed)) {
	GDate *date = g_date_new_julian(ed);

	ret = (int) g_date_strftime(s, (gsize) slen, format, date);
	g_date_free(date);
    }

    return ret;
}

/**
 * gretl_alt_strfdate:
 * @s: target string.
 * @slen: length of target string.
 * @ed: days since 1 CE.
 * @julian: 1 to use Julian calendar, 0 for Gregorian.
 *
 * If @ed is found to be valid, writes a string representing
 * the date of @ed to @s, as ISO 8601 extended.
 *
 * Returns: The number of characters written to @s, or 0 in case
 * of invalid input.
 */

int gretl_alt_strfdate (char *s, int slen, int julian,
			guint32 ed)
{
    int ret = 0;
    int y, m, d;
    int err;

    if (julian) {
	err = julian_ymd_bits_from_epoch_day(ed, &y, &m, &d);
    } else {
	err = ymd_bits_from_epoch_day(ed, &y, &m, &d);
    }

    if (!err) {
	sprintf(s, "%04d-%02d-%02d", y, m, d);
	ret = 10;
    }

    return ret;
}

static GDateTime *date_time_from_unix_offset (gint64 t, int off_secs)
{
    GDateTime *gdt0;

    gdt0 = g_date_time_new_from_unix_utc(t);

    if (off_secs == 0) {
	return gdt0;
    } else {
	GTimeZone *gtzo;
	GDateTime *gdt1;

#if GLIB_MINOR_VERSION < 58
	/* g_time_zone_new_offset() not available */
	int pm = off_secs < 0 ? -1 : 1;
	int a = abs(off_secs);
	int h = a / 3600;
	int m = (a - 3600*h) / 60;
	int s = a - 3600*h - 60*m;
	gchar *tmp;

	h *= pm;
	tmp = g_strdup_printf("%+02d:%02d:%02d", h, m, s);
	gtzo = g_time_zone_new(tmp);
	g_free(tmp);
#else
	gtzo = g_time_zone_new_offset(off_secs);
#endif
	gdt1 = g_date_time_to_timezone(gdt0, gtzo);
	g_date_time_unref(gdt0);
	g_time_zone_unref(gtzo);

	return gdt1;
    }
}

/**
 * gretl_strftime:
 * @s: target buffer.
 * @slen: length of target buffer.
 * @format: as per strftime(), or "8601" for ISO, or NULL for "%c".
 * @t: Unix time as 64-bit integer.
 * @off_secs: offset in seconds relative to UTC, or #NADBL.
 *
 * If @t and @format are found to be valid, writes a string representing
 * the date and time of @t to @s, governed by @format. If @off_secs is
 * #NADBL it is ignored, and the date and time on output are relative to
 * local time. Otherwise date and time are relative to the time zone
 * defined by @off_secs (and so a zero value gives UTC).
 *
 * Returns: The number of characters written to @s, or 0 in case
 * of invalid input.
 */

int gretl_strftime (char *s, int slen, const char *format,
		    gint64 t, double off_secs)
{
    GDateTime *gdt = NULL;
    int iso = 0;
    int ret = 0;

    if (format == NULL) {
	format = "%c";
    } else if (!strcmp(format, "8601")) {
	iso = 1;
    }

    if (!na(off_secs)) {
	int ioff = (int) floor(off_secs);

	if (abs(ioff) <= 12 * 3600) {
	    gdt = date_time_from_unix_offset(t, ioff);
	}
    } else {
	gdt = g_date_time_new_from_unix_local(t);
    }

    if (gdt != NULL) {
	gchar *tmp;

	if (iso) {
#if GLIB_MINOR_VERSION < 62
	    tmp = g_date_time_format(gdt, "%Y-%m-%dT%H%M%d%z");
#else
	    tmp = g_date_time_format_iso8601(gdt);
#endif
	} else {
	    tmp = g_date_time_format(gdt, format);
	}
	if (tmp != NULL) {
	    *s = '\0';
	    strncat(s, tmp, slen - 1);
	    ret = strlen(s);
	    g_free(tmp);
	}
	g_date_time_unref(gdt);
    }

    return ret;
}

/* Returns Unix time based on the information in @tm, inflected
   by time-zone information in @zs.
*/

static gint64 get_unix_time_GLib (struct tm *tm, const char *zs)
{
    GDateTime *gdt;
    GTimeZone *gtz;
    gdouble sec;
    gint64 t = 0;

    sec = tm->tm_sec >= 60 ? 59 : tm->tm_sec;
    gtz = g_time_zone_new(zs);

    gdt = g_date_time_new(gtz,
			  tm->tm_year + 1900,
			  tm->tm_mon + 1,
			  tm->tm_mday,
			  tm->tm_hour,
			  tm->tm_min,
			  sec);

    if (gdt != NULL) {
	t = g_date_time_to_unix(gdt);
	g_date_time_unref(gdt);
    }

    g_time_zone_unref(gtz);

    return t;
}

/* Returns a copy of @fmt with the time-zone format "%z" lopped off:
   by removing it we can get strptime() to return the time-zone
   portion of a date/time string as unprocessed remainder, which
   we can then pass to GLib to get a GTimeZone.
*/

static gchar *get_adjusted_format (const char *fmt, int *got_tz)
{
    gchar *ret = g_strdup(fmt);
    gchar *p = strstr(ret, "%z");

    if (p != NULL) {
	*p = '\0';
	*got_tz = 1;
    }

    return ret;
}

/* Wrapper for strptime() which handles time-zone information in
   the date/time string @s, if present. "Returns" Unix time as a
   double via the last argument.
*/

char *gretl_strptime (const char *s, const char *format, double *dt)
{
    struct tm tm = {0,0,0,1,0,0,0,0,-1};
    int got_tz = 0;
    char *rem;

    if (strstr(format, "%z") != NULL) {
	gchar *adjfmt = get_adjusted_format(format, &got_tz);

	rem = strptime(s, adjfmt, &tm);
	g_free(adjfmt);
    } else {
	rem = strptime(s, format, &tm);
    }

    if (got_tz && rem != NULL && *rem != '\0') {
	*dt = (double) get_unix_time_GLib(&tm, rem);
    } else {
#ifdef WIN32
	*dt = win32_mktime(&tm);
#else
	*dt = (double) mktime(&tm);
#endif
    }

    return rem;
}
