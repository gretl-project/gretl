#ifndef LAGSELECT_H
#define LAGSELECT_H

enum {
    LAG_X = 1,    /* lags set for regular variable context */
    LAG_Y_X,      /* lags for dependent variable */
    LAG_W,        /* lags set for variable as instrument */
    LAG_Y_W,      /* lags for dependent var as instrument */
    LAG_Y_V       /* lags of endoenous vars in VAR */
} LagContext;

#define VDEFLT -1

/* FIXME globals */
extern selector *open_selector;
extern int y_x_lags_enabled;
extern int y_w_lags_enabled;

void destroy_lag_preferences (void);

int set_lag_prefs_from_list (int v, int *llist, char context,
			     int *changed);

int set_lag_prefs_from_minmax (int v, int lmin, int lmax,
			       char context, int *changed);

void set_null_lagpref (int v, char context, int *changed);

int set_lag_prefs_from_model (int dv, int *xlist, int *zlist);

int set_lag_prefs_from_VAR (const int *lags, int *xlist);

int *get_lag_pref_as_list (int v, char context);

int remove_specific_lag (int v, int lag, char context);

int is_lag_dummy (int v, int lag, char context);

const int *get_VAR_lags_list (void);

void set_VAR_max_lag (int lmax);

void get_lag_preference (int v, int *lmin, int *lmax, const int **laglist,
			 char context, selector *sr);

#endif
