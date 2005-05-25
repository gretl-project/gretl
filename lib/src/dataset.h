/*
 *  Copyright (c) 2005 by Allin Cottrell
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

#ifndef DATASET_H
#define DATASET_H

/**
 * dataset_is_time_series:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains time series
 * data (1) or not (0).
 */
#define dataset_is_time_series(p) ((p)->structure == TIME_SERIES || \
				   (p)->structure == SPECIAL_TIME_SERIES)

/**
 * custom_time_series:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains time series
 * data with custom (non-standard) frequency (1) or not (0).
 */
#define custom_time_series(p) ((p)->structure == SPECIAL_TIME_SERIES)

/**
 * dataset_is_daily:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains daily time series
 * data (1) or not (0).
 */
#define dataset_is_daily(p) (p->structure == TIME_SERIES \
                             && (p->pd == 5 || p->pd == 6 || p->pd == 7))

/**
 * dataset_is_weekly:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains weekly time series
 * data (1) or not (0).
 */
#define dataset_is_weekly(p) (p->structure == TIME_SERIES \
                              && p->pd == 52)

/**
 * dataset_is_hourly:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains hourly time series
 * data (1) or not (0).
 */
#define dataset_is_hourly(p) (p->structure == TIME_SERIES \
                              && p->pd == 24)

/**
 * dataset_is_decennial:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains decennial time series
 * data (1) or not (0).
 */
#define dataset_is_decennial(p) (p->structure == TIME_SERIES \
                                 && p->pd == 10)

/**
 * dated_daily_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains dated daily time series
 * data (1) or not (0).
 */
#define dated_daily_data(p) (p->structure == TIME_SERIES \
                             && (p->pd == 5 || p->pd == 6 || p->pd == 7) \
                             && p->sd0 > 10000.0)

/**
 * dated_seven_day_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains dated daily 
 * (seven-day) time series data (1) or not (0).
 */
#define dated_seven_day_data(p) (p->structure == TIME_SERIES \
                                 && p->pd == 7 && \
                                 p->sd0 > 10000.0)

/**
 * dated_weekly_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains dated weekly 
 * time series data (1) or not (0).
 */
#define dated_weekly_data(p) (p->structure == TIME_SERIES \
                              && p->pd == 52 && \
                              p->sd0 > 10000.0)

/**
 * calendar_data:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set uses calendar
 * dates for observation strings (1) or not (0).
 */
#define calendar_data(p) (p->structure == TIME_SERIES && \
                          (p->pd == 5 || p->pd == 6 || p->pd == 7 \
                           || p->pd == 52) && p->sd0 > 10000.0) 
                          
/**
 * dataset_is_panel:
 * @p: pointer to data information struct.
 *
 * Attempt to determine whether a data set contains panel
 * data (1) or not (0).
 */
#define dataset_is_panel(p) ((p)->structure == STACKED_TIME_SERIES || \
                             (p)->structure == STACKED_CROSS_SECTION)


int dataset_add_observations (int newobs, double ***pZ, DATAINFO *pdinfo);

int dataset_add_series (int newvars, double ***pZ, DATAINFO *pdinfo);

int dataset_add_allocated_series (double *x, double ***pZ, 
				  DATAINFO *pdinfo);

int dataset_add_scalar (double ***pZ, DATAINFO *pdinfo);

int dataset_scalar_to_vector (int v, double ***pZ, DATAINFO *pdinfo);

int dataset_drop_listed_variables (const int *list, double ***pZ, 
				   DATAINFO *pdinfo, int *renumber);

int dataset_destroy_hidden_variables (double ***pZ, DATAINFO *pdinfo);

int dataset_drop_last_variables (int delvars, double ***pZ, DATAINFO *pdinfo);

int dataset_stack_variables (double ***pZ, DATAINFO *pdinfo, 
			     char *newvar, char *s);

int is_hidden_variable (int i, const DATAINFO *pdinfo);


#endif /* DATASET_H */
