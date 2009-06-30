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

#define N_COMMON_OPTS 5

enum tx_objects {
    D11,      /* seasonally adjusted series */
    D12,      /* trend/cycle */
    D13,      /* irregular component */
    TRIGRAPH, /* graph showing all of the above */
    TEXTOUT,  /* for full text output */
    XAXIS     /* x-axis (time) variable for graphing */
};

typedef struct _common_opt_info common_opt_info;
typedef struct _tx_request tx_request;

struct _common_opt_info {
    GtkWidget *check;
    char save;
    unsigned short v;
};

struct _tx_request {
    int code;          /* tramo vs x12arima */
    GtkWidget *dialog;
    common_opt_info opt[N_COMMON_OPTS];
    void *opts;
    int savevars;
    int pd;
};

int show_tramo_options (tx_request *request, GtkWidget *vbox);

int print_tramo_options (tx_request *request, FILE *fp);
