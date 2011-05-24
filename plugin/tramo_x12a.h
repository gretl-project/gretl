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

enum tx_objects {
    TX_SA,    /* save seasonally adjusted series */
    TX_TR,    /* save trend/cycle */
    TX_IR,    /* save irregular component */
    TRIGRAPH, /* graph showing all of the above */
    TEXTOUT,  /* for full text output */
    TX_MAXOPT
};

typedef struct _common_opt_info common_opt_info;
typedef struct _x12a_opts x12a_opts;
typedef struct _tx_request tx_request;

struct _common_opt_info {
    GtkWidget *check;
    char save;
    unsigned short v;
    char savename[VNAMELEN];
};

struct _x12a_opts {
    int logtrans;
    int outliers;
    int trdays;
};    

struct _tx_request {
    int prog;          /* tramo vs x12arima */
    GtkWidget *dialog;
    void (*helpfunc);
    common_opt_info opts[TX_MAXOPT];
    char yname[VNAMELEN];
    void *gui;
    gretlopt *popt;
    int savevars;
    int pd;
    int seasonal_ok;
    x12a_opts xopt;
};

int add_tramo_options (tx_request *request, GtkWidget *vbox);

int print_tramo_options (tx_request *request, FILE *fp);

const char *get_tramo_save_string (int i);

void sensitize_tx_entry (GtkToggleButton *b, GtkWidget *w);

void update_tx_savename (GtkEntry *entry, char *name);
