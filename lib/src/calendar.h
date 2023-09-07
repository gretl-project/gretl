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

#ifndef CALENDAR_H
#define CALENDAR_H

/**
 * FOUR_DIGIT_YEAR:
 * @y: year given to only two digits.
 *
 * Returns: a guess at the 4-digit year intended when
 * two digits were provided. But really, who knows?
 */
#define FOUR_DIGIT_YEAR(y) ((y < 50)? y + 2000 : y + 1900)

guint32 epoch_day_from_ymd (int y, int m, int d);

guint32 epoch_day_from_ymd_basic (double ymd);

guint32 epoch_day_from_julian_ymd (int y, int m, int d);

guint32 nearby_epoch_day (int y, int m, int d, int wkdays);

char *ymd_extended_from_epoch_day (guint32 ed, int julian, int *err);

double ymd_basic_from_epoch_day (guint32 ed, int julian, int *err);

int ymd_bits_from_epoch_day (guint32 ed, int *y, int *m, int *d);

int julian_ymd_bits_from_epoch_day (guint32 ed, int *y, int *m, int *d);

int iso_basic_to_extended (const double *b, double *y, double *m, 
			   double *d, int n);

guint32 get_epoch_day (const char *datestr);

guint32 epoch_day_from_t (int t, const DATASET *dset);

int weekday_from_date (const char *datestr);

int weekday_from_epoch_day (guint32 ed);

int day_starts_month (int d, int m, int y, int wkdays, int *pad);

int day_ends_month (int d, int m, int y, int wkdays);

int get_days_in_month (int m, int y, int wkdays, int julian);

int month_day_index (int y, int m, int d, int wkdays);

int days_in_month_before (int y, int m, int d, int wkdays);

int days_in_month_after (int y, int m, int d, int wkdays);

int date_to_daily_index (const char *datestr, int wkdays);

int daily_index_to_date (char *targ, int y, int m, int idx,
			 int wkdays);

int calendar_obs_number (const char *datestr, const DATASET *dset,
			 int nolimit);

int calendar_date_string (char *targ, int t, const DATASET *dset);

int MS_excel_date_string (char *targ, int mst, int pd, int d1904);

double get_dec_date (const char *datestr);

double legacy_day_of_week (int y, int m, int d, int julian, int *err);

int n_hidden_missing_obs (const DATASET *dset, int t1, int t2);

int guess_daily_pd (const DATASET *dset);

double easterdate (int year);

int day_span (guint32 ed1, guint32 ed2, int wkdays, int *err);

int iso_week_number (int y, int m, int d, int *err);

int iso_week_from_date (const char *datestr);

int fill_monthlen_array (double *mlen, int t1, int t2,
			 int wkdays, int mo, int yr,
			 const double *movec,
			 const double *yrvec,
			 int julian);

int gretl_strfdate (char *s, int slen, const char *format,
		    guint32 ed);

int gretl_alt_strfdate (char *s, int slen, int julian,
			guint32 ed);

int gretl_strftime (char *s, int slen, const char *format,
		    gint64 t, double off_secs);

char *gretl_strptime (const char *s, const char *format,
		      double *dt);

#endif /* CALENDAR_H */ 
