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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* lib.c for gretl -- main interface to libgretl functions */

#include "gretl.h"
#ifdef G_OS_WIN32 
# include "../lib/src/cmdlist.h"
# include <io.h>
#else
# include <unistd.h>
#endif

#include "selector.h"

extern DATAINFO *subinfo;
extern DATAINFO *fullinfo;
extern double **subZ;
extern double **fullZ;

/* ../cli/common.c */
static int data_option (int flag);
extern int loop_exec_line (LOOPSET *plp, const int round, 
			   const int cmdnum, PRN *prn);
/* boxplots.c */
extern int boxplots (int *list, char **bools, 
		     double ***pZ, const DATAINFO *pdinfo, 
		     int notches);
extern int boolean_boxplots (const char *str, double ***pZ, 
			     DATAINFO *pdinfo, int notches);

/* private functions */
static int gui_exec_line (char *line, 
			  LOOPSET *plp, int *plstack, int *plrun, 
			  SESSION *psession, SESSIONBUILD *rebuild,
			  PRN *prn, int exec_code, 
			  const char *myname); 
static void console_exec (void);
static void finish_genr (MODEL *pmod);
static gint stack_model (int gui);

int replay;                 /* shared, to indicate whether we're just
			       replaying old session commands or not */
int echo_off;               /* don't echo commands */

GtkItemFactoryEntry log_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>" },    
    { N_("/File/_Save As..."), NULL, file_save, SAVE_CMDS, NULL },
    { N_("/File/_Run"), NULL, do_run_script, SESSION_EXEC, NULL },
#if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("/File/_Print..."), NULL, window_print, 0, NULL },
#endif
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, NULL },
    { N_("/Edit/Copy _all"), NULL, text_copy, COPY_TEXT, NULL },
    { NULL, NULL, NULL, 0, NULL }
};

GtkItemFactoryEntry console_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>" }, 
    { N_("/File/Save _As..."), NULL, file_save, SAVE_CONSOLE, NULL },
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, NULL },
    { N_("/Edit/Copy _all"), NULL, text_copy, COPY_TEXT, NULL },
    { NULL, NULL, NULL, 0, NULL }
};

GtkItemFactoryEntry script_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>" }, 
    { N_("/File/Save _As..."), NULL, file_save, SAVE_SCRIPT, NULL },
    { N_("/File/_Run"), NULL, do_run_script, SCRIPT_EXEC, NULL },
#if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("/File/_Print..."), NULL, window_print, 0, NULL },
#endif
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, NULL },
    { N_("/Edit/Copy _all"), NULL, text_copy, COPY_TEXT, NULL },
    { N_("/Edit/_Paste"), NULL, text_paste, 0, NULL },
    { N_("/Edit/_Replace..."), NULL, text_replace, 0, NULL },
    { N_("/Edit/_Undo"), NULL, text_undo, 0, NULL },
    { NULL, NULL, NULL, 0, NULL }
};

GtkItemFactoryEntry sample_script_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>" },    
    { N_("/File/_Save As..."), NULL, file_save, SAVE_SCRIPT, NULL },
    { N_("/File/_Run"), NULL, do_run_script, SCRIPT_EXEC, NULL },
#if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("/File/_Print..."), NULL, window_print, 0, NULL },
#endif
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, NULL },
    { N_("/Edit/Copy _all"), NULL, text_copy, COPY_TEXT, NULL },
    { NULL, NULL, NULL, 0, NULL }
};

GtkItemFactoryEntry script_out_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>" },    
    { N_("/File/Save _As..."), NULL, file_save, SAVE_OUTPUT, NULL },
#if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("/File/_Print..."), NULL, window_print, 0, NULL },
#endif
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, NULL },
    { N_("/Edit/Copy _all"), NULL, text_copy, COPY_TEXT, NULL },
    { NULL, NULL, NULL, 0, NULL }
};

GtkItemFactoryEntry view_items[] = {
#if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("/_File"), NULL, NULL, 0, "<Branch>" },     
    { N_("/File/_Print..."), NULL, window_print, 0, NULL },
#endif
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, NULL },
    { N_("/Edit/Copy _all"), NULL, text_copy, COPY_TEXT, NULL },
    { NULL, NULL, NULL, 0, NULL }
};

const char *CANTDO = N_("Can't do this: no model has been estimated yet\n");

typedef struct {
    int ID, cmdnum;
    int n;
    char **cmds;
} model_stack;

/* file scope state variables */
static GtkWidget *console_view;
static int ignore;
static int oflag;
static char loopstorefile[MAXLEN];
static char errline[MAXLEN] = ""; /* for use with the console */
static model_stack *mstack;
static int n_mstacks;
static int model_count;
static char **cmd_stack;
static int n_cmds;
static MODELSPEC *modelspec;
static char *model_origin;
static char last_model = 's';

/* ........................................................... */

void register_graph (void)
{
#ifdef GNUPLOT_PNG
    extern int gnuplot_show_png(const char *fname);

    gnuplot_show_png(paths.plotfile);
#else
    graphmenu_state(TRUE);
#endif    
}

/* ........................................................... */

int quiet_sample_check (MODEL *pmod)
{
    double **checkZ;
    DATAINFO *pdinfo;

    if (fullZ == NULL) {
	checkZ = Z;
	pdinfo = datainfo;
    } else {
	checkZ = fullZ;
	pdinfo = fullinfo;
    }

    if (checkZ == NULL || pdinfo == NULL) return 1;

    if (model_sample_issue(pmod, NULL, checkZ, pdinfo)) return 1;
    else return 0;
}

/* ......................................................... */

static void set_sample_label_special (void)
{
    char labeltxt[80];

    sprintf(labeltxt, _("Undated: Full range n = %d; current sample"
	    " n = %d"), fullinfo->n, datainfo->n);
    gtk_label_set_text(GTK_LABEL(mdata->status), labeltxt);

}

/* ........................................................... */

void free_modelspec (void)
{
    int i = 0;

    if (modelspec != NULL) {
	while (modelspec[i].cmd != NULL) {
	    free(modelspec[i].cmd);
	    if (modelspec[i].subdum != NULL)
		free(modelspec[i].subdum);
	    i++;
	}
	free(modelspec);
	modelspec = NULL;
    }
}

/* ........................................................... */

void free_command_stack (void)
{
    int i, j;

    if (cmd_stack != NULL) {
	for (i=0; i<n_cmds; i++)
	    if (cmd_stack[i]) free(cmd_stack[i]);
	free(cmd_stack);
	cmd_stack = NULL;
    }
    n_cmds = 0;

    if (n_mstacks > 0 && mstack != NULL) {  
	for (i=0; i<n_mstacks; i++) {
	    for (j=0; j<mstack[i].n; j++)
		free(mstack[i].cmds[j]); 
	}
	free(mstack);
	mstack = NULL;
    }
    n_mstacks = 0;
}

/* ........................................................... */

void clear_data (int full)
{
    extern void clear_clist (GtkWidget *widget);

    if (full) 
	clear(paths.datfile, MAXLEN);
    restore_sample(NULL, 0, NULL);
    if (Z != NULL) free_Z(Z, datainfo); 
    clear_datainfo(datainfo, CLEAR_FULL);
    Z = NULL;
    fullZ = NULL;
    clear_clist(mdata->listbox);
    clear_sample_label();
    data_status = 0;
    orig_vars = 0;
    menubar_state(FALSE);

    /* clear everything out */
    clear_model(models[0], NULL);
    clear_model(models[1], NULL);
    clear_model(models[2], NULL);

    free_command_stack(); 
    free_modelspec();
    modelspec = NULL;

    stack_model(-1);
    model_count = 0;
}

/* ........................................................... */

char *user_fopen (const char *fname, char *fullname, PRN **pprn)
{
    strcpy(fullname, paths.userdir);
    strcat(fullname, fname);

    *pprn = gretl_print_new(GRETL_PRINT_FILE, fullname);
    if (*pprn == NULL) {
	errbox(_("Couldn't open file for writing"));
	return NULL;
    }
    return fullname;
}

/* ........................................................... */

gint bufopen (PRN **pprn)
{
    *pprn = gretl_print_new (GRETL_PRINT_BUFFER, NULL);
    if (*pprn == NULL) {
	errbox(_("Out of memory allocating output buffer"));
	return 1;
    }
    return 0;
}

/* ........................................................... */

PRN *bufopen_with_size (size_t sz)
{
    PRN *prn;

    prn = malloc(sizeof *prn);
    if (prn == NULL) {
	errbox(_("Out of memory allocating output buffer"));
	return NULL;
    }
    prn->fp = NULL;
    prn->buf = malloc(sz);
    if (prn->buf == NULL) {
	errbox(_("Out of memory allocating output buffer"));
	free(prn);
	return NULL;
    }
    return prn;
}

/* ........................................................... */

static int freq_error (FREQDIST *freq, PRN *prn)
{
    if (freq == NULL) {
	if (prn == NULL)
	    errbox(_("Out of memory in frequency distribution"));
	else
	    pprintf(prn, _("Out of memory in frequency distribution\n"));
	return 1;
    }
    if (get_gretl_errno()) {
	if (prn == NULL)
	    gui_errmsg(get_gretl_errno());
	else
	    errmsg(get_gretl_errno(), prn);
	free_freq(freq);
	return 1;
    }
    return 0;
}

/* ........................................................... */

gint check_cmd (char *line)
{
    strcpy(command.param, "");
    catchflag(line, &oflag);
    getcmd(line, datainfo, &command, &ignore, &Z, NULL); 
    if (command.errcode) {
	gui_errmsg(command.errcode);
	return 1;
    } 
    replay = 0; /* we're not just replaying saved session commands */
    return 0;
}

/* ........................................................... */

static void maybe_quote_filename (char *line, char *cmd)
{
    size_t len = strlen(cmd);

    if (strlen(line) > len + 1) {
	char *p = line + len + 1;

	if (*p == '"' || *p == '\'') return;
	
	if (strchr(p, ' ')) {
	    char tmp[MAXLEN];

	    *tmp = 0;
	    strcpy(tmp, p);
	    sprintf(line, "%s \"%s\"", cmd, tmp);
	}
    }
}

/* ........................................................... */

gint cmd_init (char *line)
{
    size_t len;
    PRN *echo;

#ifdef CMD_DEBUG
    fprintf(stderr, "cmd_init: got line: '%s'\n", line);
    fprintf(stderr, "command.cmd: '%s'\n", command.cmd);
    fprintf(stderr, "command.param: '%s'\n", command.param);
#endif

    if (n_cmds == 0) 
	cmd_stack = mymalloc(sizeof *cmd_stack);
    else 
	cmd_stack = myrealloc(cmd_stack, (n_cmds + 1) * sizeof *cmd_stack);
    if (cmd_stack == NULL) return 1;

    if (command.ci == OPEN || command.ci == RUN) 
	maybe_quote_filename(line, command.cmd);

    if (bufopen(&echo)) return 1;

    echo_cmd(&command, datainfo, line, 0, 1, oflag, echo);

    len = strlen(echo->buf);
    if ((cmd_stack[n_cmds] = mymalloc(len + 1)) == NULL)
	return 1;
    strcpy(cmd_stack[n_cmds], echo->buf);
#ifdef CMD_DEBUG
    fprintf(stderr, "cmd_init: cmd_stack[%d] now contains: '%s'\n", 
	    n_cmds, cmd_stack[n_cmds]);
#endif
    gretl_print_destroy(echo);
    n_cmds++;

    return 0;
}

/* ........................................................... */

static gint record_model_genr (char *line)
{
    size_t len = strlen(line);

    if (n_cmds == 0) 
	cmd_stack = mymalloc(sizeof(char *));
    else 
	cmd_stack = myrealloc(cmd_stack, (n_cmds + 1) * sizeof(char *));
    if (cmd_stack == NULL) return 1;

    if ((cmd_stack[n_cmds] = mymalloc(len + 1)) == NULL)
	return 1;
    strncpy(cmd_stack[n_cmds], line, len);
    n_cmds++;

    return 0;
}

/* ........................................................... */

static int grow_mstack (int i, int model_id)
{
    if (n_mstacks == 0) { 
#ifdef CMD_DEBUG
	fprintf(stderr, "grow_mstack: starting from scratch\n");
#endif
	mstack = mymalloc(sizeof *mstack);
    } else { 
#ifdef CMD_DEBUG
	fprintf(stderr, "grow_mstack: reallocating to %d stacks\n",
		n_mstacks+1);
#endif
	mstack = myrealloc(mstack, (n_mstacks+1) * sizeof *mstack);
    }
    if (mstack == NULL) {
	n_mstacks = 0;
	return 1;
    }
    mstack[i].ID = model_id;    
    mstack[i].cmdnum = n_cmds-1;
    mstack[i].n = 0;
#ifdef CMD_DEBUG
    fprintf(stderr, "mstack[%d]: ID=%d, cmdnum=%d\n", i, model_id, n_cmds-1);
#endif
    mstack[i].cmds = mymalloc(sizeof(char **));
    if (mstack[i].cmds == NULL) return 1;
    n_mstacks++;
#ifdef CMD_DEBUG
    fprintf(stderr, "grow_mstack: n_mstacks now = %d\n", n_mstacks);
#endif
    return 0;
}

/* ........................................................... */

static gint model_cmd_init (char *line, int ID)
     /* this makes a record of commands associated with
	a given model, so that they may be reconstructed later as
	part of the session mechanism */
{
    int i, sn;
    PRN *echo;
    size_t len;

    /* have we started this stuff at all, yet? */
    if (n_mstacks == 0) { /* no */
	if (grow_mstack(0, ID)) {
	    free(mstack);
	    return 1;
	}
    }

    /* have we already started a stack for this model? */
    sn = -1;
    for (i=0; i<n_mstacks; i++) {
	if (mstack[i].ID == ID) { /* yes */
	    sn = i;
	    break;
	}
    }
    if (sn == -1) { /* no, not yet */
	sn = n_mstacks;
	if (grow_mstack(sn, ID)) { 
	    free(mstack);
	    return 1;
	}
    } 

    if (mstack[sn].n > 0) { /* stack already underway for this model; 
			       make space for another command string */
#ifdef CMD_DEBUG
	fprintf(stderr, "model_cmd_init: realloc mstack[%d] for %d cmds\n",
		sn, mstack[sn].n+1);
#endif
	mstack[sn].cmds = myrealloc(mstack[sn].cmds,
				    (mstack[sn].n+1) * sizeof(char **));
	if (mstack[sn].cmds == NULL) {
	    /* do more stuff! */
	    return 1;
	}
    }

    if (bufopen(&echo)) return 1;
    echo_cmd(&command, datainfo, line, 0, 1, oflag, echo);

    len = strlen(echo->buf);

    mstack[sn].cmds[mstack[sn].n] = mymalloc(len + 1);
    if (mstack[sn].cmds[mstack[sn].n] == NULL) {
	gretl_print_destroy(echo);
	return 1;
    }
    strcpy(mstack[sn].cmds[mstack[sn].n], echo->buf);
    gretl_print_destroy(echo);

    mstack[sn].n += 1;

    return 0;
}

/* ........................................................... */

static gint stack_model (int gui)
{
    static int m;

    if (gui == -1) { /* code for reset, when changing datasets */
	m = 0;
	return 0;
    }

    /* record the way this model was estimated (GUI or not) */
    if (model_origin == NULL) 
	model_origin = malloc(sizeof *model_origin);
    else
	model_origin = myrealloc(model_origin, 
				 model_count * sizeof *model_origin);
    if (model_origin == NULL) return 1;
    last_model = (gui == 1)? 'g' : 's';
    model_origin[model_count - 1] = last_model;

    if (!gui) { /* Model estimated via console or script: unlike a gui
		   model, which is kept in memory so long as its window
		   is open, these models are immediately discarded.  So
		   if we want to be able to refer back to them later we
		   need to record their specification */
	if (modelspec == NULL) 
	    modelspec = mymalloc(2 * sizeof *modelspec);
	else 
	    modelspec = myrealloc(modelspec, (m+2) * sizeof *modelspec);
	if (modelspec == NULL) return 1;
	else {
	    modelspec[m].cmd = mymalloc(MAXLEN);
	    modelspec[m].subdum = NULL;
	    modelspec[m+1].cmd = NULL;
	    modelspec[m+1].subdum = NULL;
	    if (fullZ != NULL) {
		fullinfo->varname = datainfo->varname;
		fullinfo->label = datainfo->label;
		fullinfo->vector = datainfo->vector;
		attach_subsample_to_model(models[0], &fullZ, fullinfo);
	    }
	    save_model_spec(models[0], &modelspec[m], fullinfo);
	    m++;
	}
    }
    return 0;
}

/* ........................................................... */

static void dump_model_cmds (FILE *fp, int m)
{
    int i;

    fprintf(fp, "(* commands pertaining to model %d *)\n", mstack[m].ID);
    for (i=0; i<mstack[m].n; i++) 
	fprintf(fp, "%s", mstack[m].cmds[i]);
}

/* ........................................................... */

gint dump_cmd_stack (const char *fname)
     /* ship out the stack of commands entered in the current
	session */
{
    FILE *fp;
    int i, j;

    if (fname == NULL) return 0;

    if (!strcmp(fname, "stderr")) {
	fp = stderr;
	fprintf(fp, "dumping command stack:\n");
    } else {
	fp = fopen(fname, "w"); 
	if (fp == NULL) {
	    errbox(_("Couldn't open command file for writing"));
	    return 1;
	}
    }

    for (i=0; i<n_cmds; i++) {
	fprintf(fp, "%s", cmd_stack[i]);
	if (is_model_cmd(cmd_stack[i]) && mstack != NULL) {
#ifdef CMD_DEBUG
	    fprintf(stderr, "cmd_stack[%d]: looking for model commands\n", i);
#endif
	    for (j=0; j<n_mstacks; j++) { 
		if (mstack[j].cmdnum == i) {
		   dump_model_cmds(fp, j);
		   break;
		} 
	    }
	}
    }

    if (strcmp(fname, "stderr")) 
	fclose(fp);

    return 0;
}

/* ........................................................... */

void do_menu_op (gpointer data, guint action, GtkWidget *widget)
{
    PRN *prn;
    char title[48];
    int err = 0;
    windata_t *vwin;
    gpointer obj = NULL;
    gint hsize = 78, vsize = 380;

    clear(line, MAXLEN);
    strcpy(title, "gretl: ");

    switch (action) {
    case CORR:
	strcpy(line, "corr");
	strcat(title, _("correlation matrix"));
	break;
    case FREQ:
	sprintf(line, "freq %s", datainfo->varname[mdata->active_var]);
	strcat(title, _("frequency distribution"));
	vsize = 340;
	break;
    case RUNS:
	sprintf(line, "runs %s", datainfo->varname[mdata->active_var]);
	strcat(title, _("runs test"));
	vsize = 200;
	break;
    case SUMMARY:
	strcpy(line, "summary");
	strcat(title, _("summary statistics"));
	break;
    case VAR_SUMMARY:
	sprintf(line, "summary %s", datainfo->varname[mdata->active_var]);
	strcat(title, _("summary stats: "));
	strcat(title, datainfo->varname[mdata->active_var]);
	vsize = 300;
	break;
    default:
	break;
    }

    /* check the command and initialize output buffer */
    if (check_cmd(line) || cmd_init(line) || bufopen(&prn)) return;

    /* execute the command */
    switch (action) {
    case CORR:
	obj = corrlist(command.list, &Z, datainfo);
	if (obj == NULL) {
	    errbox(_("Failed to generate correlation matrix"));
	    gretl_print_destroy(prn);
	    return;
	} 
	/* printcorr(corr, datainfo, &prn); */
	matrix_print_corr(obj, datainfo, 0, prn);
	break;
    case FREQ:
	obj = freqdist(&Z, datainfo, mdata->active_var, 1);
	if (freq_error(obj, NULL)) {
	    gretl_print_destroy(prn);
	    return;
	} 
	printfreq(obj, prn);
	free_freq(obj);
	break;
    case RUNS:
	err = runs_test(command.list[1], Z, datainfo, prn);
	break;
    case SUMMARY:
    case VAR_SUMMARY:	
	obj = summary(command.list, &Z, datainfo, prn);
	if (obj == NULL) {
	    errbox(_("Failed to generate summary statistics"));
	    gretl_print_destroy(prn);
	    return;
	}	    
	print_summary(obj, datainfo, 0, prn);
	break;
    }
    if (err) gui_errmsg(err);

    vwin = view_buffer(prn, hsize, vsize, title, action, view_items);

    if (vwin && 
	(action == SUMMARY || action == VAR_SUMMARY || action == CORR)) 
	vwin->data = obj;
}

/* ........................................................... */

void do_coint (GtkWidget *widget, gpointer p)
{
    selector *sr = (selector *) p;
    char *edttext;
    PRN *prn;
    int err = 0, order = 0;

    edttext = sr->cmdlist;
    if (*edttext == 0) return;

    clear(line, MAXLEN);
    sprintf(line, "coint %s", edttext);

    /* check the command and initialize output buffer */
    if (check_cmd(line) || cmd_init(line) || bufopen(&prn)) return;

    order = atoi(command.param);
    if (!order) {
	errbox(_("Couldn't read cointegration order"));
	gretl_print_destroy(prn);
	return;
    }

    err = coint(order, command.list, &Z, datainfo, prn);
    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
	return;
    } 

    view_buffer(prn, 78, 400, _("gretl: cointegration test"), 
		COINT, view_items);
}

/* ........................................................... */

void do_dialog_cmd (GtkWidget *widget, dialog_t *ddata)
{
    char *edttext;
    PRN *prn;
    char title[48];
    int err = 0, order = 0, mvar = mdata->active_var;
    gint hsize = 78, vsize = 300;

    edttext = gtk_entry_get_text (GTK_ENTRY(ddata->edit));
    if (*edttext == '\0' && ddata->code != CORRGM) return;

    clear(line, MAXLEN);
    strcpy(title, "gretl: ");

    /* set up the command */
    switch (ddata->code) {
    case ADF:
	sprintf(line, "adf %s %s", edttext, datainfo->varname[mvar]);
	strcat(title, _("adf test"));
	vsize = 350;
	break;
    case COINT:
	sprintf(line, "coint %s", edttext);
	strcat(title, _("cointegration test"));
	vsize = 400;
	break;
    case SPEARMAN:
	sprintf(line, "spearman -o %s", edttext);
	strcat(title, _("rank correlation"));
	vsize = 400;
	break;
    case MEANTEST:
	sprintf(line, "meantest -o %s", edttext);
	strcat(title, _("means test"));
	break;
    case MEANTEST2:
	sprintf(line, "meantest %s", edttext);
	strcat(title, _("means test"));
	break;
    case VARTEST:
	sprintf(line, "vartest %s", edttext);
	strcat(title, _("variances test"));
	break;
    case CORRGM:
	if (*edttext != '\0') order = atoi(edttext);
	if (order) 
	    sprintf(line, "corrgm %s %d", 
		    datainfo->varname[mvar], order);
	else
	    sprintf(line, "corrgm %s", 
		    datainfo->varname[mvar]);
	strcat(title, _("correlogram"));
	break;
    default:
	dummy_call();
	return;
    }

    /* check the command and initialize output buffer */
    if (check_cmd(line) || cmd_init(line) || bufopen(&prn)) return;

    /* execute the command */
    switch (ddata->code) {
    case ADF:
    case COINT:
	order = atoi(command.param);
	if (!order) {
	    errbox((ddata->code == ADF)? 
		   _("Couldn't read ADF order") :
		   _("Couldn't read cointegration order"));
	    gretl_print_destroy(prn);
	    return;
	}
	if (ddata->code == ADF)
	    err = adf_test(order, command.list[1], &Z, datainfo, prn);
	else
	    err = coint(order, command.list, &Z, datainfo, prn);
	break;
    case SPEARMAN:
	err = spearman(command.list, Z, datainfo, 1, prn);
	break;
    case MEANTEST:
	err = means_test(command.list, Z, datainfo, 1, prn);;
	break;
    case MEANTEST2:
	err = means_test(command.list, Z, datainfo, 0, prn);;
	break;
    case VARTEST:
	err = vars_test(command.list, Z, datainfo, prn);
	break;	
    case CORRGM:
	err = corrgram(command.list[1], order, &Z, datainfo, &paths, 0, prn);
	break;
    default:
	dummy_call();
	return;
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	if (ddata->code == CORRGM) 
	    register_graph();
	view_buffer(prn, hsize, vsize, title, ddata->code, view_items);
    }
}

/* ........................................................... */

void open_info (gpointer data, guint edit, GtkWidget *widget)
{
    if (datainfo->descrip == NULL) {
	if (!yes_no_dialog(_("gretl: add info"), 
			   _("The data file contains no informative comments.\n"
			   "Would you like to add some now?"), 0)) {
	    edit_header(NULL, 0, NULL);
	}
    } else {
	PRN *prn;
	size_t sz = strlen(datainfo->descrip);

	prn = bufopen_with_size(sz + 1);
	if (prn != NULL) { 
	    strcpy(prn->buf, datainfo->descrip);
	    view_buffer(prn, 80, 400, _("gretl: data info"), INFO, view_items);
	}
    }
}

/* ........................................................... */

void view_log (void)
{
    char fname[MAXLEN];

    strcpy(fname, paths.userdir);
    strcat(fname, "session.inp");

    if (dump_cmd_stack(fname)) return;

    view_file(fname, 0, 0, 78, 370, VIEW_LOG, log_items);
}

/* ........................................................... */

void console (void)
{
    PRN *prn;
    char fname[MAXLEN];
    windata_t *vwin;

    if (console_view != NULL) {
	gdk_window_show(console_view->parent->window);
	gdk_window_raise(console_view->parent->window);
	return;
    }

    if (!user_fopen("console_tmp", fname, &prn)) return;

    pprintf(prn, _("gretl console: type 'help' for a list of commands\n? "));
    gretl_print_destroy(prn);
    vwin = view_file(fname, 1, 0, 78, 400, CONSOLE, console_items);
    console_view = vwin->w;
    gtk_signal_connect(GTK_OBJECT(console_view), "destroy",
		       GTK_SIGNAL_FUNC(gtk_widget_destroyed),
		       &console_view);
    
    gtk_text_set_point(GTK_TEXT(console_view), 52);

    GTK_TEXT(console_view)->cursor_pos_x = 
	2 * gdk_char_width(fixed_font, 'X');
    gtk_editable_set_position(GTK_EDITABLE(console_view), 52);
    gtk_editable_changed(GTK_EDITABLE(console_view));
    gtk_widget_grab_focus (console_view);

}

/* ........................................................... */

static int backkey (GdkEventKey *key)
{
    if (key->keyval == GDK_BackSpace || 
	key->keyval == GDK_Left) return 1;
    return 0;
}

/* ........................................................... */

static int last_console_line_len (void)
{
    int i, c, len;

    len = gtk_text_get_length(GTK_TEXT(console_view));
    for (i=len; i>0; i--) {
	c = GTK_TEXT_INDEX(GTK_TEXT(console_view), i - 1);
	if (c == '\n') break; 
    }
    return len - i - 2;
}

/* ........................................................... */

gboolean console_handler (GtkWidget *w, GdkEventKey *key, gpointer d)
{
    static int lastkey;
    int cw = gdk_char_width(fixed_font, 'X'); 
    int len = gtk_text_get_length (GTK_TEXT(console_view));
    int adjust, currpos, linelen, xpos, savekey = key->keyval;
    GdkModifierType mods;

    /* null action if not at prompt */
    currpos = GTK_EDITABLE(console_view)->current_pos;
    xpos = GTK_TEXT(console_view)->cursor_pos_x / cw; 
    linelen = last_console_line_len();

    if ((currpos < len - linelen) || 
	(backkey(key) && xpos < 3)) {
	gtk_signal_emit_stop_by_name(GTK_OBJECT(w), "key-press-event");
	return TRUE;
    }
    /* make return key execute the command */
    if (key->keyval == GDK_Return) {
	console_exec();
	key->keyval = GDK_End;
	return FALSE;
    }
    /* make up-arrow recall last command entered */
    if (key->keyval == GDK_Up) {
#ifdef CMD_DEBUG
	fprintf(stderr, "errline: '%s', n_cmds = %d\n", errline, n_cmds);
#endif
	key->keyval = GDK_VoidSymbol;
	if (lastkey != GDK_Up) {
	    if (errline[0] != '\0') {
		gtk_editable_insert_text(GTK_EDITABLE(console_view), 
					 errline, 
					 strlen(errline),
					 &len);
	    } else if (n_cmds > 0 && 
		       strncmp(cmd_stack[n_cmds-1], "open", 4)) {
		gtk_editable_insert_text(GTK_EDITABLE(console_view), 
					 cmd_stack[n_cmds-1], 
					 strlen(cmd_stack[n_cmds-1]) - 1,
					 &len);
	    }
	}
    }
    /* down-arrow clears line */
    if (lastkey == GDK_Up && key->keyval == GDK_Down) {
	key->keyval = GDK_VoidSymbol;
	adjust = GTK_TEXT(console_view)->cursor_pos_x / cw - 2;
	gtk_editable_delete_text(GTK_EDITABLE(console_view), 
				 len - adjust, len);
    }  
    /* response to Ctrl-A: go to start of typing area */
    gdk_window_get_pointer(console_view->window, NULL, NULL, &mods);
    if (mods & GDK_CONTROL_MASK && 
	gdk_keyval_to_upper(key->keyval) == GDK_A) {
	gtk_editable_set_position(GTK_EDITABLE(console_view), 
				  len - last_console_line_len());	
	gtk_signal_emit_stop_by_name(GTK_OBJECT(w), "key-press-event");
	return TRUE;
    }
    lastkey = savekey;
    return FALSE;
}

/* ........................................................... */

void console_exec (void)
{
    PRN *prn;
    int len, loopstack = 0, looprun = 0;
    gchar *c_line; 
    char execline[MAXLEN];
    extern GdkColor red;

    len = gtk_text_get_length (GTK_TEXT(console_view));
    c_line = gtk_editable_get_chars(GTK_EDITABLE(console_view), 
				    len - last_console_line_len(), len);
    top_n_tail(c_line);

    if (strcmp(c_line, "quit") == 0 || strcmp(c_line, "q") == 0) {
	gtk_widget_destroy(console_view->parent->parent->parent);
	g_free(c_line);
	return;
    }

    if (bufopen(&prn)) {
	g_free(c_line);
	return;
    }

    strncpy(execline, c_line, MAXLEN - 1);
    g_free(c_line);
    gui_exec_line(execline, NULL, &loopstack, &looprun, NULL, NULL, 
		  prn, CONSOLE_EXEC, NULL);

    /* put results into console window */
    gtk_text_freeze(GTK_TEXT(console_view));
    gtk_editable_insert_text(GTK_EDITABLE(console_view), 
			     "\n", 1, &len);

    gtk_text_insert(GTK_TEXT(console_view), fixed_font, 
		    NULL, NULL, prn->buf, strlen(prn->buf));
    gretl_print_destroy(prn);

    gtk_text_insert(GTK_TEXT(console_view), fixed_font,
		    &red, NULL, "\n? ", 3);
    gtk_text_thaw(GTK_TEXT(console_view));
    len = gtk_text_get_length (GTK_TEXT(console_view));
    gtk_editable_set_position(GTK_EDITABLE(console_view), len);
}

/* ........................................................... */

void gui_errmsg (const int errcode)
{
    char *msg = get_gretl_errmsg();

    if (msg[0] != '\0') 
	errbox(msg);
    else {
	msg = get_errmsg(errcode, errtext, NULL);
	if (msg != NULL)
	    errbox(msg);
	else
	    errbox(_("Unspecified error"));
    }
}

/* ........................................................... */

void change_sample (GtkWidget *widget, dialog_t *ddata) 
{
    char *edttext;
    int err;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;

    clear(line, MAXLEN);
    sprintf(line, "smpl %s", edttext);
    if (check_cmd(line) || cmd_init(line)) return;

    err = set_sample(line, datainfo);
    if (err) gui_errmsg(err);
    else {
	set_sample_label(datainfo);
	restore_sample_state(TRUE);
    }
}
/* ........................................................... */

void bool_subsample (gpointer data, guint opt, GtkWidget *w)
     /* opt = 0     -- drop all obs with missing data values 
	opt = OPT_O -- sample using dummy variable
	opt = OPT_R -- sample using boolean expression
     */
{
    int err = 0;

    restore_sample(NULL, 0, NULL);
    if ((subinfo = mymalloc(sizeof *subinfo)) == NULL) 
	return;

    if (opt == 0)
	err = set_sample_dummy(NULL, &Z, &subZ, datainfo, subinfo, OPT_O);
    else
	err = set_sample_dummy(line, &Z, &subZ, datainfo, subinfo, opt);
    if (err) {
	gui_errmsg(err);
	return;
    }

    /* save the full data set for later use */
    fullZ = Z;
    fullinfo = datainfo;
    datainfo = subinfo;
    Z = subZ;

    set_sample_label_special();
    restore_sample_state(TRUE);
    if (opt == 0)
	infobox(_("Sample now includes only complete observations"));
    else
	infobox(_("Sub-sampling done"));
}

/* ........................................................... */

void do_samplebool (GtkWidget *widget, dialog_t *ddata)
{
    char *edttext = NULL;

    edttext = gtk_entry_get_text(GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;

    clear(line, MAXLEN);
    sprintf(line, "smpl %s -r", edttext); 
    if (check_cmd(line) || cmd_init(line)) return;

    bool_subsample(NULL, OPT_R, NULL);
}

/* ........................................................... */

void do_sampledum (GtkWidget *widget, dialog_t *ddata)
{
    char *edttext = NULL, dumv[9];

    edttext = gtk_entry_get_text(GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;
    sscanf(edttext, "%8s", dumv);
    dumv[8] = '\0';
	
    clear(line, MAXLEN);
    sprintf(line, "smpl %s -o", dumv);
    if (check_cmd(line) || cmd_init(line)) return;

    bool_subsample(NULL, OPT_O, NULL);
}

/* ........................................................... */

void do_setobs (GtkWidget *widget, dialog_t *ddata)
{
    char *edttext, pdstr[8], stobs[9];
    int err, opt;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;

    sscanf(edttext, "%7s %8s", pdstr, stobs);
	
    clear(line, MAXLEN);
    sprintf(line, "setobs %s %s ", pdstr, stobs);
    catchflag(line, &opt);
    if (check_cmd(line) || cmd_init(line)) return;

    err = set_obs(line, datainfo, opt);
    if (err) {
	errbox(get_gretl_errmsg());
	return;
    } else {
	char msg[80];

	sprintf(msg, _("Set data frequency to %d, starting obs to %s"),
		datainfo->pd, datainfo->stobs);
	infobox(msg);
	set_sample_label(datainfo);
	data_status |= MODIFIED_DATA;
    }
}

/* ........................................................... */

void count_missing (void)
{
    PRN *prn;

    if (bufopen(&prn)) return;
    if (count_missing_values(&Z, datainfo, prn)) {
	view_buffer(prn, 77, 300, _("gretl: missing values info"), 
		    SMPL, view_items);
    } else {
	infobox(_("No missing data values"));
	gretl_print_destroy(prn);
    }
}

/* ........................................................... */

void do_add_markers (GtkWidget *widget, dialog_t *ddata) 
{
    char *edttext;
    char fname[MAXLEN];

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;
    strcpy(fname, edttext);

    if (add_case_markers(datainfo, fname)) 
	errbox(_("Failed to add case markers"));
    else {
	infobox(_("Case markers added"));
	data_status |= MODIFIED_DATA; 
    }
}

/* ........................................................... */

void do_forecast (GtkWidget *widget, dialog_t *ddata) 
{
    windata_t *mydata = ddata->data;
    MODEL *pmod = mydata->data;
    char *edttext;
    PRN *prn;
    int err;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;
    
    clear(line, MAXLEN);
    sprintf(line, "fcasterr %s", edttext);
    if (check_cmd(line) || cmd_init(line) || bufopen(&prn)) return;

    err = fcast_with_errs(line, pmod, &Z, datainfo, prn,
			  &paths, 1); 
    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
	return;
    }

    register_graph();

    view_buffer(prn, 78, 350, _("gretl: forecasts"), FCAST, NULL);    
}

/* ........................................................... */

void do_add_omit (GtkWidget *widget, gpointer p)
{
    selector *sr = (selector *) p;
    windata_t *vwin = sr->data;
    char *edttext;
    PRN *prn;
    char title[26];
    MODEL *orig, *pmod;
    gint err;

    orig = vwin->data;
    edttext = sr->cmdlist;
    if (*edttext == '\0') return;
    
    clear(line, MAXLEN);
    if (sr->code == ADD) 
        sprintf(line, "addto %d %s", orig->ID, edttext);
    else 
        sprintf(line, "omitfrom %d %s", orig->ID, edttext);

    if (check_cmd(line) || bufopen(&prn)) return;

    pmod = gretl_model_new(datainfo);
    if (pmod == NULL) {
	errbox(_("Out of memory"));
	gretl_print_destroy(prn);
	return;
    }

    if (sr->code == ADD) 
        err = auxreg(command.list, orig, pmod, &model_count, 
                     &Z, datainfo, AUX_ADD, prn, NULL);
    else 
        err = omit_test(command.list, orig, pmod, &model_count, 
			&Z, datainfo, prn);

    if (err) {
        gui_errmsg(err);
        gretl_print_destroy(prn);
        clear_model(pmod, NULL); 
        return;
    }

    if (cmd_init(line) || stack_model(1)) {
	errbox(_("Error saving model information"));
	return;
    }

    /* update copy of most recently estimated model */
    if (copy_model(models[2], pmod, datainfo))
	errbox(_("Out of memory copying model"));

    /* record sub-sample info (if any) with the model */
    if (fullZ != NULL) {
	fullinfo->varname = datainfo->varname;
	fullinfo->label = datainfo->label;	
	attach_subsample_to_model(pmod, &fullZ, fullinfo);
    }

    sprintf(title, _("gretl: model %d"), model_count);
    view_model(prn, pmod, 78, 400, title);
}

/* ........................................................... */

static gint add_test_to_model (GRETLTEST *test, MODEL *pmod)
{
    int i, nt = pmod->ntests;

    if (nt == 0) {
        pmod->tests = malloc(sizeof(GRETLTEST));
    } else {
        for (i=0; i<nt; i++) {
            if (strcmp(test->type, pmod->tests[i].type) == 0) {
                return -1;
	    }
	}
        pmod->tests = myrealloc(pmod->tests, (nt + 1) * sizeof(GRETLTEST));
    }
    if (pmod->tests == NULL) return 1;

    strcpy(pmod->tests[nt].type, test->type);
    strcpy(pmod->tests[nt].h_0, test->h_0);
    pmod->tests[nt].teststat = test->teststat;
    pmod->tests[nt].value = test->value;
    pmod->tests[nt].dfn = test->dfn;
    pmod->tests[nt].dfd = test->dfd;
    pmod->tests[nt].pvalue = test->pvalue;

    pmod->ntests += 1;

    return 0;
}

/* ........................................................... */

static void print_test_to_window (GRETLTEST *test, GtkWidget *w)
{
    if (w == NULL) {
        return;
    } else {
        char test_str[64], pval_str[64];
        gchar *tempstr;

        get_test_stat_string (test, test_str, GRETL_PRINT_FORMAT_PLAIN);
        get_test_pval_string (test, pval_str, GRETL_PRINT_FORMAT_PLAIN);

        tempstr = g_strdup_printf("%s -\n"
                                  "  %s: %s\n"
                                  "  %s: %s\n"
                                  "  %s = %s\n\n",
                                  _(test->type), 
                                  _("Null hypothesis"), _(test->h_0), 
                                  _("Test statistic"), test_str, 
                                  _("with p-value"), pval_str);

	gtk_text_freeze(GTK_TEXT (w));
	gtk_text_insert(GTK_TEXT (w), fixed_font, NULL, NULL, tempstr, 
			strlen(tempstr));
	gtk_text_thaw(GTK_TEXT (w));
	g_free(tempstr);
    }
}

/* ........................................................... */

void do_lmtest (gpointer data, guint aux_code, GtkWidget *widget)
{
    int err;
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    PRN *prn;
    char title[40];
    GRETLTEST test;

    if (bufopen(&prn)) return;
    strcpy(title, _("gretl: LM test "));
    clear(line, MAXLEN);

    if (aux_code == AUX_WHITE) {
	strcpy(line, "lmtest -c");
	err = whites_test(pmod, &Z, datainfo, prn, &test);
	if (err) {
	    gui_errmsg(err);
	    gretl_print_destroy(prn);
	    return;
	} else {
	    strcat(title, _("(heteroskedasticity)"));
	    if (add_test_to_model(&test, pmod) == 0)
		print_test_to_window(&test, mydata->w);
	}
    } 
    else {
	if (aux_code == AUX_SQ) 
	    strcpy(line, "lmtest -s");
	else
	    strcpy(line, "lmtest -l");
	clear_model(models[0], NULL);
	err = auxreg(NULL, pmod, models[0], &model_count, 
		     &Z, datainfo, aux_code, prn, &test);
	if (err) {
	    gui_errmsg(err);
	    clear_model(models[0], NULL);
	    gretl_print_destroy(prn);
	    return;
	} else {
	    clear_model(models[0], NULL); 
	    model_count--;
	    strcat(title, _("(non-linearity)"));
	    if (add_test_to_model(&test, pmod) == 0)
		print_test_to_window(&test, mydata->w);
	} 
    }

    if (check_cmd(line) || model_cmd_init(line, pmod->ID)) return;

    view_buffer(prn, 77, 400, title, LMTEST, view_items); 
}

/* ........................................................... */

void gui_set_panel_structure (gpointer data, guint u, GtkWidget *w)
{
    extern GtkWidget *open_dialog;
    extern void panel_structure_dialog (DATAINFO *, GtkWidget *,
					void (*)(), void (*)());    

    if (open_dialog != NULL) {
	gdk_window_raise(open_dialog->window);
	return;
    }

    panel_structure_dialog(datainfo, open_dialog, 
			   destroy_dialog_data, context_help);
}

/* ........................................................... */

void do_panel_diagnostics (gpointer data, guint u, GtkWidget *w)
{
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    void *handle;
    void (*panel_diagnostics)(MODEL *, double ***, DATAINFO *, PRN *);
    PRN *prn;

    if (!balanced_panel(datainfo)) {
	errbox(_("Sorry, can't do this test on an unbalanced panel.\n"
	       "You need to have the same number of observations\n"
	       "for each cross-sectional unit"));
	return;
    }

    if (gui_open_plugin("panel_data", &handle)) return;
    panel_diagnostics = get_plugin_function("panel_diagnostics", handle);

    if (panel_diagnostics == NULL) {
	errbox(_("Couldn't load plugin function"));
	close_plugin(handle);
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return;
    }	
	
    (*panel_diagnostics)(pmod, &Z, datainfo, prn);

    close_plugin(handle);

    view_buffer(prn, 77, 400, _("gretl: panel model diagnostics"), 
		PANEL, view_items);
}



/* ........................................................... */

static void do_chow_cusum (gpointer data, int code)
{
    windata_t *mydata;
    dialog_t *ddata = NULL;
    MODEL *pmod;
    char *edttext;
    PRN *prn;
    GRETLTEST test;
    gint err;

    if (code == CHOW) {
	ddata = (dialog_t *) data;
	mydata = ddata->data;
    } else
	mydata = (windata_t *) data;

    pmod = mydata->data;
    if (pmod->ci != OLS) {
	errbox(_("This test only implemented for OLS models"));
	return;
    }

    if (code == CHOW) {
	edttext = gtk_entry_get_text (GTK_ENTRY(ddata->edit));
	if (*edttext == '\0') return;
	clear(line, MAXLEN);
	sprintf(line, "chow %s", edttext);
    } else 
	strcpy(line, "cusum");

    if (bufopen(&prn)) return;

    if (code == CHOW)
	err = chow_test(line, pmod, &Z, datainfo, prn, &test);
    else
	err = cusum_test(pmod, &Z, datainfo, prn, &paths, &test);
    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
	return;
    } 

    if (add_test_to_model(&test, pmod) == 0)
	print_test_to_window(&test, mydata->w);

    if (check_cmd(line) || model_cmd_init(line, pmod->ID))
	return;

    view_buffer(prn, 77, 400, (code == CHOW)?
		_("gretl: Chow test output"): _("gretl: CUSUM test output"),
		code, view_items);
}

/* ........................................................... */

void do_chow (GtkWidget *widget, dialog_t *ddata)
{
    do_chow_cusum((gpointer) ddata, CHOW);
}    


/* ........................................................... */

void do_cusum (gpointer data, guint u, GtkWidget *widget)
{
    do_chow_cusum(data, CUSUM);
}

/* ........................................................... */

void do_autocorr (GtkWidget *widget, dialog_t *ddata)
{
    windata_t *mydata = ddata->data;
    MODEL *pmod = mydata->data;
    GRETLTEST test;
    char *edttext;
    PRN *prn;
    char title[40];
    int order, err;

    edttext = gtk_entry_get_text (GTK_ENTRY(ddata->edit));
    if (*edttext == '\0') return;

    order = atoi(edttext);

    if (bufopen(&prn)) return;
    strcpy(title, _("gretl: LM test (autocorrelation)"));

    clear(line, MAXLEN);
    sprintf(line, "lmtest -m %d", order);

    err = autocorr_test(pmod, order, &Z, datainfo, prn, &test);
    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
	return;
    } else {
	if (add_test_to_model(&test, pmod) == 0)
	    print_test_to_window(&test, mydata->w);
    }

    if (check_cmd(line) || model_cmd_init(line, pmod->ID)) return;

    view_buffer(prn, 77, 400, title, LMTEST, view_items); 
}

/* ........................................................... */

void do_arch (GtkWidget *widget, dialog_t *ddata)
{
    windata_t *mydata = ddata->data;
    MODEL *pmod = mydata->data;
    GRETLTEST test;
    char *edttext;
    PRN *prn;
    char tmpstr[26];
    int order, err, i;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;

    clear(line, MAXLEN);
    sprintf(line, "arch %s ", edttext);
    for (i=1; i<=pmod->list[0]; i++) {
	sprintf(tmpstr, "%d ", pmod->list[i]);
	strcat(line, tmpstr);
    }
    if (check_cmd(line) || cmd_init(line)) return;

    order = atoi(command.param);
    if (!order) {
	errbox(_("Couldn't read ARCH order"));
	return;
    }

    if (bufopen(&prn)) return;

    clear_model(models[1], NULL);
    exchange_smpl(pmod, datainfo);
    *models[1] = arch(order, pmod->list, &Z, datainfo, 
		     NULL, prn, &test);
    if ((err = (models[1])->errcode)) 
	errmsg(err, prn);
    else {
	if (add_test_to_model(&test, pmod) == 0)
	    print_test_to_window(&test, mydata->w);
	if (oflag) outcovmx(models[1], datainfo, 0, prn);
    }
    clear_model(models[1], NULL);
    exchange_smpl(pmod, datainfo);

    view_buffer(prn, 78, 400, _("gretl: ARCH test"), ARCH, view_items);
}

/* ........................................................... */

void set_storelist (GtkWidget *widget, dialog_t *ddata)
{
    char *edttext;
    size_t len;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (edttext == NULL || *edttext == '\0') return;

    len = strlen(edttext);
    if (len >= MAXLEN || (storelist && !strcmp(storelist, edttext))) {
	free(storelist);
	storelist = NULL; /* default list -- don't need to be explicit */
    } else
	strcpy(storelist, edttext);
}

/* ........................................................... */

static int model_error (const MODEL *pmod)
{
    if (pmod->errcode) {
	gui_errmsg(pmod->errcode);
	return 1;
    }
    return 0;
}

/* ........................................................... */

static int model_output (MODEL *pmod, PRN *prn)
{
    if (model_error(pmod)) return 1;

    ++model_count;
    pmod->ID = model_count;
    if (printmodel(pmod, datainfo, prn))
	pmod->errcode = E_NAN; /* some statistics were NAN */

    return 0;
}

/* ........................................................... */

static gint check_model_cmd (char *line, char *modelgenr)
{
    PRN *getgenr;

    if (bufopen(&getgenr)) return 1;
    strcpy(command.param, "");
    catchflag(line, &oflag);
    getcmd(line, datainfo, &command, &ignore, &Z, getgenr); 
    if (command.errcode) {
	gui_errmsg(command.errcode);
	return 1;
    }
    if (strlen(getgenr->buf)) strcpy(modelgenr, getgenr->buf);
    gretl_print_destroy(getgenr);
    return 0;
}

/* ........................................................... */

#ifdef ENABLE_GMP

void do_mp_ols (GtkWidget *widget, gpointer p)
{
    char *edttext, estimator[9];
    void *handle;
    int (*mplsq)(const int *, const int *,
		 double ***, DATAINFO *, PRN *, char *, mp_results *);
    int err, action;
    selector *sr = (selector *) p;
    PRN *prn;

    action = sr->code;
    strcpy(estimator, commands[action]);

    edttext = sr->cmdlist;    
    if (*edttext == '\0') return;

    clear(line, MAXLEN);
    sprintf(line, "%s %s", estimator, edttext);

    if (check_cmd(line) || cmd_init(line) || bufopen(&prn)) return;

    if (gui_open_plugin("mp_ols", &handle)) return;
    mplsq = get_plugin_function("mplsq", handle);

    if (mplsq == NULL) {
	errbox(_("Couldn't load plugin function"));
	close_plugin(handle);
	return;
    }

    err = (*mplsq)(command.list, NULL, &Z, datainfo, prn, errtext, NULL);

    close_plugin(handle);

    if (err) {
	if (errtext[0] != 0) errbox(errtext);
	else errbox(get_errmsg(err, errtext, NULL));
	gretl_print_destroy(prn);
	return;
    }

    view_buffer(prn, 78, 400, _("gretl: high precision estimates"), 
                MPOLS, view_items);
}

#endif /* ENABLE_GMP */

/* ........................................................... */

void do_model (GtkWidget *widget, gpointer p) 
{
    char *edttext;
    PRN *prn;
    char title[26], estimator[9], modelgenr[80];
    int order, err = 0, action;
    double rho;
    MODEL *pmod = NULL;
    selector *sr = (selector *) p;  

    action = sr->code;
    strcpy(estimator, commands[action]);

    edttext = sr->cmdlist;    
    if (*edttext == '\0') return;

    clear(line, MAXLEN);
    sprintf(line, "%s %s", estimator, edttext);
    modelgenr[0] = '\0';
    if (check_model_cmd(line, modelgenr)) return;
    echo_cmd(&command, datainfo, line, 0, 1, oflag, NULL);
    if (command.ci == 999) {
	errbox(_("A variable was duplicated in the list of regressors"));
	return;
    }

    if (bufopen(&prn)) return;

    if (action != VAR) {
	pmod = gretl_model_new(datainfo);
	if (pmod == NULL) {
	    errbox(_("Out of memory"));
	    return;
	}
    }

    switch (action) {

    case CORC:
    case HILU:
	err = hilu_corc(&rho, command.list, &Z, datainfo, action, prn);
	if (err) {
	    errmsg(err, prn);
	    break;
	}
	*pmod = lsq(command.list, &Z, datainfo, action, 1, rho);
	err = model_output(pmod, prn);
	break;

    case OLS:
    case WLS:
    case POOLED:
	*pmod = lsq(command.list, &Z, datainfo, action, 1, 0.0);
	if ((err = model_output(pmod, prn))) break;
	if (oflag) outcovmx(pmod, datainfo, 0, prn);
	break;

    case HSK:
	*pmod = hsk_func(command.list, &Z, datainfo);
	if ((err = model_output(pmod, prn))) break;
	if (oflag) outcovmx(pmod, datainfo, 0, prn);
	break;

    case HCCM:
	*pmod = hccm_func(command.list, &Z, datainfo);
	if ((err = model_output(pmod, prn))) break;
	if (oflag) print_white_vcv(pmod, prn);
	break;

    case TSLS:
	*pmod = tsls_func(command.list, atoi(command.param), 
				&Z, datainfo);
	if ((err = model_output(pmod, prn))) break;
	if (oflag) outcovmx(pmod, datainfo, 0, prn);
	break;

    case AR:
	*pmod = ar_func(command.list, atoi(command.param), 
			      &Z, datainfo, &model_count, prn);
	if ((err = model_error(pmod))) break;
	if (oflag) outcovmx(pmod, datainfo, 0, prn);
	break;

    case VAR:
	/* requires special treatment: doesn't return model */
	sscanf(edttext, "%d", &order);
	err = var(order, command.list, &Z, datainfo, 0, prn);
	if (err) errmsg(err, prn);
	view_buffer(prn, 78, 450, _("gretl: vector autoregression"), 
		    VAR, view_items);
	return;

    case LOGIT:
    case PROBIT:
	*pmod = logit_probit(command.list, &Z, datainfo, action);
	err = model_output(pmod, prn);
	break;	

    default:
	errbox(_("Sorry, not implemented yet!"));
	break;
    }

    if (err) {
	gretl_print_destroy(prn);
	return;
    }

    if (modelgenr[0] && record_model_genr(modelgenr)) {
	errbox(_("Error saving model information"));
	return;
    }
    if (cmd_init(line) || stack_model(1)) {
	errbox(_("Error saving model information"));
	return;
    }
    /* make copy of most recent model */
    if (copy_model(models[2], pmod, datainfo))
	errbox(_("Out of memory copying model"));

    /* record sub-sample info (if any) with the model */
    if (fullZ != NULL) {
	fullinfo->varname = datainfo->varname;
	fullinfo->label = datainfo->label;	
	attach_subsample_to_model(pmod, &fullZ, fullinfo);
    }
    
    /* record the fact that the last model was estimated via GUI */
    sprintf(title, _("gretl: model %d"), pmod->ID);

    /* fprintf(stderr, "do_model: calling view_model\n"); */
    view_model(prn, pmod, 78, 400, title); 
}

/* ........................................................... */

void do_sim (GtkWidget *widget, dialog_t *ddata)
{
    char *edttext, varname[9], info[24];
    int err;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;

    clear(line, MAXLEN);
    sprintf(line, "sim %s", edttext);
    if (check_cmd(line) || cmd_init(line)) return;

    sscanf(line, "%*s %*s %*s %s", varname);
    sprintf(info, _("%s redefined OK"), varname);

    err = simulate(line, &Z, datainfo);
    if (err) gui_errmsg(err);
    else infobox(info);
} 

/* ........................................................... */

void do_simdata (GtkWidget *widget, dialog_t *ddata) 
{
    char *edttext;
    int err, nulldata_n;
    PRN *prn;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;

    clear(line, MAXLEN);
    sprintf(line, "nulldata %s", edttext);
    if (check_cmd(line) || cmd_init(line)) return;

    nulldata_n = atoi(command.param);
    if (nulldata_n < 2) {
	errbox(_("Data series length missing or invalid"));
	return;
    }
    if (nulldata_n > 1000000) {
	errbox(_("Data series too long"));
	return;
    }
    
    prn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);
    if (prn == NULL) return;
    err = open_nulldata(&Z, datainfo, data_status, nulldata_n, prn);
    if (err) { 
	errbox(_("Failed to create empty data set"));
	return;
    }

    err = plotvar(&Z, datainfo, "index");
    if (err) 
	errbox(_("Error generating index variable"));

    infobox(prn->buf);
    gretl_print_destroy(prn);
    paths.datfile[0] = '\0';
    populate_clist(mdata->listbox, datainfo);
    data_status = HAVE_DATA | GUI_DATA | MODIFIED_DATA;
    set_sample_label(datainfo);
    orig_vars = datainfo->v;
    menubar_state(TRUE);
}

/* ........................................................... */

void do_genr (GtkWidget *widget, dialog_t *ddata) 
{
    char *edttext;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;

    clear(line, MAXLEN);
    sprintf(line, "genr %s", edttext);
    if (check_cmd(line) || cmd_init(line)) return;

    finish_genr(NULL);
}

/* ........................................................... */

void do_model_genr (GtkWidget *widget, dialog_t *ddata) 
{
    char *edttext;
    windata_t *mydata = (windata_t *) ddata->data;
    MODEL *pmod = mydata->data;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;

    clear(line, MAXLEN);
    sprintf(line, "genr %s", edttext);
    if (check_cmd(line) || model_cmd_init(line, pmod->ID)) return;

    finish_genr(pmod);
}
/* ........................................................... */

void do_random (GtkWidget *widget, dialog_t *ddata) 
{
    char *edttext;
    char tmp[32], vname[9];
    float f1, f2;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;

    if (sscanf(edttext, "%s %f %f", tmp, &f1, &f2) != 3) {
	if (ddata->code == GENR_NORMAL) 
	    errbox(_("Specification is malformed\n"
		   "Should be like \"foo 1 2.5\""));
	else
	    errbox(_("Specification is malformed\n"
		   "Should be like \"foo 0 10\""));
	return;
    }
    if (ddata->code == GENR_NORMAL && f2 < 0) {
	errbox(_("Can't have a negative standard deviation!"));
	return;
    } else if (ddata->code == GENR_UNIFORM && f1 >= f2) {
	errbox(_("Range is non-positive!"));
	return;
    }
	
    strncpy(vname, tmp, 8);
    vname[8] = '\0';
    if (validate_varname(vname)) return;

    clear(line, MAXLEN);

    if (ddata->code == GENR_NORMAL) {
	if (f1 != 0. || f2 != 1.)
	    sprintf(line, "genr %s = %.3f * normal() + %.3f", 
		    vname, f2, f1);
	else sprintf(line, "genr %s = normal()", vname); 
    } else {
	if (f1 != 0. || f2 != 100.)
	    sprintf(line, "genr %s = %.3f + (uniform() * %.3f)", 
		    vname, f1, (f2 - f1)/100.);
	else sprintf(line, "genr %s = uniform()", vname); 
    }

    if (check_cmd(line) || cmd_init(line)) return;

    finish_genr(NULL);
}

/* ........................................................... */

void do_seed (GtkWidget *widget, dialog_t *ddata)
{
    char *edttext;
    char tmp[32];

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;

    sscanf(edttext, "%31s", tmp);
	
    clear(line, MAXLEN);
    sprintf(line, "seed %s", tmp); 
    if (check_cmd(line) || cmd_init(line)) return;

    srand((unsigned) atoi(tmp));
}

/* ........................................................... */

static void finish_genr (MODEL *pmod)
{
    int err;

    if (pmod != NULL)
	err = generate(&Z, datainfo, line, model_count, 
		       pmod, 0); 
    else
	err = generate(&Z, datainfo, line, model_count, 
		       (last_model == 's')? models[0] : models[2], 0); 
    if (err) {
	gui_errmsg(err);
	free(cmd_stack[n_cmds-1]);
	n_cmds--;
    } else {
	populate_clist(mdata->listbox, datainfo);
	data_status |= MODIFIED_DATA;
    }
}

/* ........................................................... */

void do_rename_var (GtkWidget *widget, dialog_t *ddata) 
{
    char *edttext;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;
    
    if (validate_varname(edttext)) return;
    strcpy(datainfo->varname[mdata->active_var], edttext);
    populate_clist(mdata->listbox, datainfo);
    data_status |= MODIFIED_DATA; 
}

/* ........................................................... */

static int real_do_setmiss (double missval, int varno) 
{
    int i, t, count = 0;
    int start = 1, end = datainfo->v;

    if (varno) {
	start = varno;
	end = varno + 1;
    }

    for (i=start; i<end; i++) {
	if (!datainfo->vector[i]) continue;
	for (t=0; t<datainfo->n; t++) {
	    if (Z[i][t] == missval) {
		Z[i][t] = NADBL;
		count++;
	    }
	}	
    }
    return count;
}

void do_global_setmiss (GtkWidget *widget, dialog_t *ddata)
{
    char *edttext;
    double missval;
    int count;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;

    missval = atof(edttext);
    count = real_do_setmiss(missval, 0);

    if (count) {
	sprintf(errtext, _("Set %d values to \"missing\""), count);
	infobox(errtext);
	data_status |= MODIFIED_DATA;
    } else 
	errbox(_("Didn't find any matching observations"));	
}

void do_variable_setmiss (GtkWidget *widget, dialog_t *ddata)
{
    char *edttext;
    double missval;
    int count;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;

    if (!datainfo->vector[mdata->active_var]) {
	errbox(_("This variable is a scalar"));
	return;
    }

    missval = atof(edttext);
    count = real_do_setmiss(missval, mdata->active_var);

    if (count) {
	sprintf(errtext, _("Set %d observations to \"missing\""), count);
	infobox(errtext);
	data_status |= MODIFIED_DATA;
    } else 
	errbox(_("Didn't find any matching observations"));
}

/* ........................................................... */

void delete_var (void)
{
    if (datainfo->v <= 1) {
	errbox(_("Can't delete last variable"));
	return;
    }
    if (dataset_drop_vars(1, &Z, datainfo)) {
	errbox(_("Failed to shrink the data set"));
	return;
    }
    populate_clist(mdata->listbox, datainfo);
    data_status |= MODIFIED_DATA; 
}

/* ........................................................... */

void do_edit_label (GtkWidget *widget, dialog_t *ddata) 
{
    char *edttext;

    edttext = gtk_entry_get_text (GTK_ENTRY (ddata->edit));
    if (*edttext == '\0') return;
    
    strncpy(datainfo->label[mdata->active_var], edttext, MAXLABEL-1);
    datainfo->label[mdata->active_var][MAXLABEL-1] = '\0';
    populate_clist(mdata->listbox, datainfo);
    data_status |= MODIFIED_DATA; 
}

/* ........................................................... */

static void normal_test (GRETLTEST *test, FREQDIST *freq)
{
    strcpy(test->type, _("Test for normality of residual"));
    strcpy(test->h_0, _("error is normally distributed"));
    test->teststat = GRETL_TEST_CHISQ;
    test->value = freq->chisqu;
    test->dfn = 2;
    test->pvalue = chisq(freq->chisqu, 2);
}

/* ........................................................... */

void do_resid_freq (gpointer data, guint action, GtkWidget *widget)
{
    FREQDIST *freq;
    PRN *prn;
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    GRETLTEST test;

    if (bufopen(&prn)) return;

    if (genr_fit_resid(pmod, &Z, datainfo, GENR_RESID, 1)) {
	errbox(_("Out of memory attempting to add variable"));
	return;
    }

    freq = freqdist(&Z, datainfo, datainfo->v - 1, pmod->ncoeff);
    dataset_drop_vars(1, &Z, datainfo);
    if (freq_error(freq, NULL)) {
	gretl_print_destroy(prn);
	return;
    }
    
    normal_test(&test, freq);

    if (add_test_to_model(&test, pmod) == 0)
	print_test_to_window(&test, mydata->w);

    clear(line, MAXLEN);
    strcpy(line, "testuhat");
    if (check_cmd(line) || model_cmd_init(line, pmod->ID)) return;
 
    printfreq(freq, prn);
    free_freq(freq);

    view_buffer(prn, 77, 300, _("gretl: residual dist."), TESTUHAT,
		NULL);
}

/* ........................................................... */

void do_freqplot (gpointer data, guint dist, GtkWidget *widget)
{
    FREQDIST *freq;

    if (mdata->active_var < 0) return;
    if (mdata->active_var == 0) {
	errbox(_("This command is not applicable to the constant"));
	return;
    }

    clear(line, MAXLEN);
    sprintf(line, "freq %s", datainfo->varname[mdata->active_var]);
    if (check_cmd(line) || cmd_init(line)) return;

    freq = freqdist(&Z, datainfo, mdata->active_var, 1);

    if (!freq_error(freq, NULL)) { 
	if (dist == GAMMA && freq->midpt[0] < 0.0 && freq->f[0] > 0) {
	    errbox(_("Data contain negative values: gamma distribution not "
		   "appropriate"));
	} else {
	    if (plot_freq(freq, &paths, dist))
		errbox(_("gnuplot command failed"));
	    else
		register_graph();
	}
	free_freq(freq);
    }
}

/* ........................................................... */

#ifdef USE_TRAMO

extern char tramodir[];

void do_tramo (gpointer data, guint opt, GtkWidget *widget)
{
    gint err;
    void *handle;
    int (*write_tramo_data) (char *, int, double **, const DATAINFO *, 
			     const char *);
    char fname[MAXLEN];

    if (!datainfo->vector[mdata->active_var]) {
	errbox(_("Can't do this analysis on a scalar"));
	return;
    }

    if (datainfo->pd == 1 || !dataset_is_time_series(datainfo)) {
	errbox(_("This analysis is applicable only to seasonal time series"));
	return;
    }

    if (gui_open_plugin("tramo-seats", &handle)) return;
    write_tramo_data = get_plugin_function("write_tramo_data", handle);

    if (write_tramo_data == NULL) {
	errbox(_("Couldn't load plugin function"));
	close_plugin(handle);
	return;
    }

    err = write_tramo_data (fname, mdata->active_var, Z, datainfo, tramodir);

    close_plugin(handle);

    if (err) {
	errbox(_("TRAMO command failed"));
	gretl_print_destroy(prn);
	return;
    }

    view_file(fname, 0, 0, 120, 500, _("gretl: TRAMO analysis"), TRAMO, NULL);
}
#endif

/* ........................................................... */

void do_range_mean (gpointer data, guint opt, GtkWidget *widget)
{
    gint err;
    void *handle;
    int (*range_mean_graph) (int, double **, const DATAINFO *, 
			     PRN *, PATHS *);
    PRN *prn;

    if (gui_open_plugin("range-mean", &handle)) return;
    range_mean_graph = get_plugin_function("range_mean_graph", handle);

    if (range_mean_graph == NULL) {
	errbox(_("Couldn't load plugin function"));
	close_plugin(handle);
	return;
    }

    if (bufopen(&prn)) {
	close_plugin(handle);
	return; 
    }

    err = range_mean_graph (mdata->active_var, Z, datainfo, 
			    prn, &paths);

    close_plugin(handle);

    if (!err) {
	gnuplot_display(&paths);
	register_graph();
    }

    view_buffer(prn, 60, 350, _("gretl: range-mean statistics"), RANGE_MEAN, NULL);
}

/* ........................................................... */

void do_pergm (gpointer data, guint opt, GtkWidget *widget)
{
    gint err;
    PRN *prn;

    if (bufopen(&prn)) return;

    clear(line, MAXLEN);
    if (opt)
	sprintf(line, "pergm %s -o", datainfo->varname[mdata->active_var]);
    else
	sprintf(line, "pergm %s", datainfo->varname[mdata->active_var]);

    if (check_cmd(line) || cmd_init(line)) {
	gretl_print_destroy(prn);
	return;
    }

    err = periodogram(command.list[1], &Z, datainfo, &paths, 0, opt, prn);
    if (err) {
	errbox(_("Periodogram command failed"));
	gretl_print_destroy(prn);
	return;
    }
    register_graph();

    view_buffer(prn, 60, 400, _("gretl: periodogram"), PERGM, NULL);
}

/* ........................................................... */

void do_coeff_intervals (gpointer data, guint i, GtkWidget *w)
{
    PRN *prn;
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;

    if (bufopen(&prn)) return;

    print_model_confints(pmod, datainfo, prn);

    view_buffer(prn, 77, 300, _("gretl: coefficient confidence intervals"), 
		CONFINT, view_items);
}

/* ........................................................... */

void do_outcovmx (gpointer data, guint action, GtkWidget *widget)
{
    PRN *prn;
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;

    if (Z == NULL || datainfo == NULL) {
	errbox(_("Data set is gone"));
	return;
    }

    if (bufopen(&prn)) return;

    if (pmod->ci == HCCM) print_white_vcv(pmod, prn);
    else outcovmx(pmod, datainfo, 0, prn); 

    view_buffer(prn, 77, 300, _("gretl: coefficient covariances"), 
		COVAR, view_items);
}

/* ......................................................... */

void add_dummies (gpointer data, guint panel, GtkWidget *widget)
{
    gint err;

    clear(line, MAXLEN);

    if (panel) {
	if (datainfo->time_series == STACKED_TIME_SERIES)
	    sprintf(line, "genr paneldum");
	else if (datainfo->time_series == STACKED_CROSS_SECTION)
	    sprintf(line, "genr paneldum -o");
	else {
	    errbox(_("Data set is not recognized as a panel.\n"
		   "Please use \"Sample/Set frequency, startobs\"."));
	    return;
	}
    } else 
	sprintf(line, "genr dummy");

    if (check_cmd(line) || cmd_init(line)) return;

    if (panel) 
	err = paneldum(&Z, datainfo, 
		       (datainfo->time_series == STACKED_TIME_SERIES)? 0 : 1);
    else 
	err = dummy(&Z, datainfo);

    if (err) gui_errmsg(err);
    else populate_clist(mdata->listbox, datainfo);
}

/* ......................................................... */

void add_time (gpointer data, guint index, GtkWidget *widget)
{
    gint err;

    clear(line, MAXLEN);
    if (index) sprintf(line, "genr index");
    else sprintf(line, "genr time");
    if (check_cmd(line) || cmd_init(line)) return;

    err = plotvar(&Z, datainfo, (index)? "index" : "time");
    if (err) 
	errbox((index)? _("Error generating index variable") : 
	       _("Error generating time trend"));
    else populate_clist(mdata->listbox, datainfo);
}

/* ......................................................... */

void add_logs_etc (GtkWidget *widget, gpointer p)
{
    gint err = 0;
    char *edttext, msg[80];
    selector *ddata = (selector *) p;

    edttext = ddata->cmdlist;
    if (*edttext == '\0') return;

    line[0] = '\0';
    msg[0] = '\0';
    sprintf(line, "%s %s", commands[ddata->code], edttext);

    if (check_cmd(line) || cmd_init(line)) return;

    if (ddata->code == LAGS)
	err = lags(command.list, &Z, datainfo);
    else if (ddata->code == LOGS) {
	/* returns number of terms created */
	err = logs(command.list, &Z, datainfo);
	if (err < command.list[0]) err = 1;
	else err = 0;
    }
    else if (ddata->code == SQUARE) {
	/* returns number of terms created */
	err = xpxgenr(command.list, &Z, datainfo, 0, 1);
	if (err <= 0) err = 1;
	else err = 0;
    } 
    else if (ddata->code == DIFF)
	err = list_diffgenr(command.list, &Z, datainfo);
    else if (ddata->code == LDIFF)
	err = list_ldiffgenr(command.list, &Z, datainfo);

    if (err) {
	if (msg[0]) errbox(msg);
	else errbox(_("Error adding variables"));
    }
    else populate_clist(mdata->listbox, datainfo);
}

/* ......................................................... */

int add_fit_resid (MODEL *pmod, const int code, const int undo)
   /* If undo = 1, don't bother with the label, don't update
   the var display in the main window, and don't add to
   command log. */
{
    if (genr_fit_resid(pmod, &Z, datainfo, code, undo)) {
	errbox(_("Out of memory attempting to add variable"));
	return 1;
    }

    if (!undo) {
	int v;
	char line[32];

	v = datainfo->v - 1;
	populate_clist(mdata->listbox, datainfo);
	if (code == 0)
	    sprintf(line, "genr %s = uhat", datainfo->varname[v]);
	else if (code == 1)
	    sprintf(line, "fcast %s", datainfo->varname[v]);
	else if (code == 2)
	    sprintf(line, "genr %s = uhat*uhat", datainfo->varname[v]);
	check_cmd(line);
	model_cmd_init(line, pmod->ID);
    }
    return 0;
}

/* ......................................................... */

void add_model_stat (MODEL *pmod, const int which)
{
    char vname[9], vlabel[MAXLABEL], cmdstr[MAXLEN];
    int i, n;

    if (dataset_add_scalar(&Z, datainfo)) {
	errbox(_("Out of memory attempting to add variable"));
	return;
    }

    i = datainfo->v - 1;
    n = datainfo->n;

    switch (which) {
    case ESS:
	sprintf(vname, "ess_%d", pmod->ID);
	sprintf(vlabel, _("error sum of squares from model %d"), 
		pmod->ID);
	Z[i][0] = pmod->ess;
	sprintf(cmdstr, "genr ess_%d = $ess", pmod->ID);
	break;
    case R2:
	sprintf(vname, "r2_%d", pmod->ID);
	sprintf(vlabel, _("R-squared from model %d"), pmod->ID);
	Z[i][0] = pmod->rsq;
	sprintf(cmdstr, "genr r2_%d = $rsq", pmod->ID);
	break;
    case TR2:
	sprintf(vname, "trsq%d", pmod->ID);
	sprintf(vlabel, _("T*R-squared from model %d"), pmod->ID);
	Z[i][0] = pmod->nobs * pmod->rsq;
	sprintf(cmdstr, "genr trsq%d = $trsq", pmod->ID);
	break;
    case DF:
	sprintf(vname, "df_%d", pmod->ID);
	sprintf(vlabel, _("degrees of freedom from model %d"), 
		pmod->ID);
	Z[i][0] = (double) pmod->dfd;
	sprintf(cmdstr, "genr df_%d = $df", pmod->ID);
	break;
    case SIGMA:
	sprintf(vname, "sgma_%d", pmod->ID);
	sprintf(vlabel, _("std err of residuals from model %d"), 
		pmod->ID);
	Z[i][0] = pmod->sigma;
	sprintf(cmdstr, "genr sgma_%d = $sigma", pmod->ID);
	break;
    case LNL:
	sprintf(vname, "lnl_%d", pmod->ID);
	sprintf(vlabel, _("log likelihood from model %d"), 
		pmod->ID);
	Z[i][0] = pmod->lnL;
	sprintf(cmdstr, "genr lnl_%d = $lnl", pmod->ID);
	break;	
    }

    strcpy(datainfo->varname[i], vname);
    strcpy(datainfo->label[i], vlabel);
    populate_clist(mdata->listbox, datainfo);
    check_cmd(cmdstr);
    model_cmd_init(cmdstr, pmod->ID);
    infobox(_("variable added"));
}

/* ........................................................... */

void resid_plot (gpointer data, guint xvar, GtkWidget *widget)
{
    int err, origv = datainfo->v, plot_list[4], lines[1];
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    int ts = dataset_is_time_series(datainfo);
    int pdum = mydata->active_var; 

    /* add residuals to data set temporarily */
    if (add_fit_resid(pmod, 0, 1)) return;

    plot_list[0] = 2;
    plot_list[1] = datainfo->v - 1; /* last var added */
    strcpy(datainfo->varname[plot_list[1]], 
	   datainfo->varname[pmod->list[1]]);

    if (xvar) { /* plot against specified xvar */
	plot_list[2] = xvar;
	lines[0] = 0;
    } else {    /* plot against obs index or time */
	err = plotvar(&Z, datainfo, (ts)? "time" : "index");
	if (err) {
	    errbox(_("Failed to add plotting index variable"));
	    dataset_drop_vars(1, &Z, datainfo);
	    return;
	}
	plot_list[2] = varindex(datainfo, (ts)? "time" : "index");
	lines[0] = (ts)? 1 : 0;
    } 

    /* plot separated by dummy variable? */
    if (pdum) {
	plot_list[0] = 3;
	plot_list[3] = pdum;
    }

    /* generate graph */
    err = gnuplot(plot_list, lines, &Z, datainfo,
		  &paths, &plot_count, 0, 1, 
		  (pdum)? OPT_RESIDZ : OPT_RESID);
    if (err < 0) errbox(_("gnuplot command failed"));
    else register_graph();
    
    dataset_drop_vars(datainfo->v - origv, &Z, datainfo);
}

/* ........................................................... */

void fit_actual_plot (gpointer data, guint xvar, GtkWidget *widget)
{
    int err, origv = datainfo->v, plot_list[4], lines[2];
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;
    int ts = dataset_is_time_series(datainfo);

    /* add fitted values to data set temporarily */
    if (add_fit_resid(pmod, 1, 1)) return;

    /* common part of list setup */
    plot_list[0] = 3;
    plot_list[1] = datainfo->v - 1; /* last var added */
    plot_list[2] = pmod->list[1];   /* depvar from regression */

    if (xvar) {  /* plot against specified xvar */
	plot_list[3] = xvar;
	lines[0] = (pmod->list[0] == 3)? 1 : 0;
	lines[1] = 0;
    } else { /* plot against obs */
	err = plotvar(&Z, datainfo, (ts)? "time" : "index");
	if (err) {
	    errbox(_("Failed to add plotting index variable"));
	    dataset_drop_vars(1, &Z, datainfo);
	    return;
	}
	plot_list[3] = varindex(datainfo, (ts)? "time" : "index");
	lines[0] = (ts)? 1 : 0; 
	lines[1] = (ts)? 1 : 0;
    } 

    err = gnuplot(plot_list, lines, &Z, datainfo,
		  &paths, &plot_count, 0, 1, OPT_FA);
    if (err < 0) errbox(_("gnuplot command failed"));
    else register_graph();

    dataset_drop_vars(datainfo->v - origv, &Z, datainfo);
}

/* ........................................................... */

#define MAXDISPLAY 4096
/* max number of observations for which we expect to be able to 
   use the buffer approach for displaying data, as opposed to
   disk file */

void display_data (gpointer data, guint u, GtkWidget *widget)
{
    int err;
    PRN *prn;

    if (datainfo->v * datainfo->n > MAXDISPLAY) { /* use file */
	char fname[MAXLEN];

	if (!user_fopen("data_display_tmp", fname, &prn)) return;

	err = printdata(NULL, &Z, datainfo, 0, 1, prn);
	gretl_print_destroy(prn);
	view_file(fname, 0, 1, 77, 350, VIEW_DATA, NULL);
    } else { /* use buffer */
	if (bufopen(&prn)) return;

	err = printdata(NULL, &Z, datainfo, 0, 1, prn);
	if (err) {
	    errbox(_("Out of memory in display buffer"));
	    gretl_print_destroy(prn);
	    return;
	}
	view_buffer(prn, 77, 350, _("gretl: display data"), PRINT, NULL);
    }
}

/* ........................................................... */

void display_selected (GtkWidget *widget, gpointer p)
{
    selector *sr = (selector *) p;
    char *edttext; 
    PRN *prn;
    int ig = 0;
    CMD prcmd;

    prcmd.list = malloc(sizeof(int));
    prcmd.param = malloc(1);
    if (prcmd.list == NULL || prcmd.param == NULL) {
	errbox(_("Out of memory!"));
	return;
    }
    edttext = sr->cmdlist;
    if (*edttext == '\0') return;

    clear(line, MAXLEN);
    sprintf(line, "print %s", edttext);
    getcmd(line, datainfo, &prcmd, &ig, &Z, NULL);
    if (prcmd.errcode) {
	gui_errmsg(prcmd.errcode);
	return;
    }    

    if (prcmd.list[0] * datainfo->n > MAXDISPLAY) { /* use disk file */
	char fname[MAXLEN];

	if (!user_fopen("data_display_tmp", fname, &prn)) return;

	printdata(prcmd.list, &Z, datainfo, 0, 1, prn);
	gretl_print_destroy(prn);
	view_file(fname, 0, 1, 77, 350, VIEW_DATA, NULL);
    } else { /* use buffer */
	int err;

	if (bufopen(&prn)) return;
	err = printdata(prcmd.list, &Z, datainfo, 0, 1, prn);
	if (err) {
	    errbox(_("Out of memory in display buffer"));
	    gretl_print_destroy(prn);
	    return;
	}
	view_buffer(prn, 77, 350, _("gretl: display data"), PRINT, NULL);
    }
    free(prcmd.list);
    free(prcmd.param);
}

/* ........................................................... */

void display_fit_resid (gpointer data, guint code, GtkWidget *widget)
{
    PRN *prn;
    int err;
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;

    if (bufopen(&prn)) return;
    err = print_fit_resid(pmod, &Z, datainfo, prn);

    if (err) {
	errbox(_("Failed to generate fitted values"));
	gretl_print_destroy(prn);
    } else 
	view_buffer(prn, 77, 350, _("gretl: display data"), PRINT, NULL);    
}

/* ........................................................... */

void do_graph_var (void)
{
    int err, lines[1];

    if (mdata->active_var < 0) return;
    clear(line, MAXLEN);
    sprintf(line, "gnuplot %s time", datainfo->varname[mdata->active_var]);
    if (check_cmd(line) || cmd_init(line)) return;

    lines[0] = 1;
    err = gnuplot(command.list, lines, &Z, datainfo,
		  &paths, &plot_count, 0, 1, 0);
    if (err == -999)
	errbox(_("No data were available to graph"));
    else if (err < 0) 
	errbox(_("gnuplot command failed"));
    else register_graph();
}

/* ........................................................... */

void do_boxplot_var (void)
{
    if (mdata->active_var < 0) return;
    clear(line, MAXLEN);
    sprintf(line, "boxplot %s", datainfo->varname[mdata->active_var]);
    if (check_cmd(line) || cmd_init(line)) return;

    if (boxplots(command.list, NULL, &Z, datainfo, 0)) 
	errbox (_("boxplot command failed"));
}

/* ........................................................... */

void do_scatters (GtkWidget *widget, gpointer p)
{
    selector *sr = (selector *) p;
    char *edttext;
    gint err; 

    edttext = sr->cmdlist;
    if (*edttext == 0) return;

    clear(line, MAXLEN);
    sprintf(line, "scatters %s", edttext);
    if (check_cmd(line) || cmd_init(line)) return;
    err = multi_scatters(command.list, atoi(command.param), &Z, 
			 datainfo, &paths);
    if (err < 0) errbox(_("gnuplot command failed"));
    else register_graph();
}

/* ........................................................... */

void do_box_graph_trad (GtkWidget *widget, dialog_t *ddata)
{
    char *edttext;
    gint err, code = ddata->code; 

    edttext = gtk_entry_get_text (GTK_ENTRY(ddata->edit));
    if (edttext == NULL || *edttext == 0) return;

    if (strchr(edttext, '('))
	err = boolean_boxplots(edttext, &Z, datainfo, (code == GR_NBOX));
    else {
	clear(line, MAXLEN);
	sprintf(line, "boxplot %s%s", (code == GR_NBOX)? "-o " : "", edttext);

	if (check_cmd(line) || cmd_init(line)) return;
	err = boxplots(command.list, NULL, &Z, datainfo, (code == GR_NBOX));
    }
    if (err) errbox(_("boxplot command failed"));
}

/* ........................................................... */

void do_box_graph (GtkWidget *widget, gpointer p)
{
    selector *sr = (selector *) p;
    char *edttext;
    gint err, code = sr->code; 

    edttext = sr->cmdlist;
    if (*edttext == 0) return;

    if (strchr(edttext, '('))
	err = boolean_boxplots(edttext, &Z, datainfo, (code == GR_NBOX));
    else {
	clear(line, MAXLEN);
	sprintf(line, "boxplot %s%s", (code == GR_NBOX)? "-o " : "", edttext);

	if (check_cmd(line) || cmd_init(line)) return;
	err = boxplots(command.list, NULL, &Z, datainfo, (code == GR_NBOX));
    }
    if (err) errbox(_("boxplot command failed"));
}

/* ........................................................... */

void do_dummy_graph (GtkWidget *widget, gpointer p)
     /* X, Y scatter with separation by dummy (factor) */
{
    selector *sr = (selector *) p;
    char *edttext;
    gint err, lines[1] = {0}; 

    edttext = sr->cmdlist;
    if (*edttext == 0) return;

    clear(line, MAXLEN);
    sprintf(line, "gnuplot -z %s", edttext);

    if (check_cmd(line) || cmd_init(line)) return;

    if (command.list[0] != 3 || 
	!isdummy(command.list[3], datainfo->t1, datainfo->t2, Z)) {
	errbox(_("You must supply three variables, the last\nof which "
	       "is a dummy variable (values 1 or 0)"));
	return;
    }

    err = gnuplot(command.list, lines, &Z, datainfo,
		  &paths, &plot_count, 0, 1, OPT_Z);

    if (err < 0) errbox(_("gnuplot command failed"));
    else register_graph();
}

/* ........................................................... */

#ifdef GNUPLOT_PNG
extern int gnuplot_show_png (char *fname); /* gpt_control.c */
#endif

void do_graph_from_selector (GtkWidget *widget, gpointer p)
{
    selector *sr = (selector *) p;
    char *edttext;
    gint i, err, *lines = NULL;
    gint imp = (sr->code == GR_IMP);

    edttext = sr->cmdlist;
    if (*edttext == '\0') return;

    clear(line, MAXLEN);
    sprintf(line, "gnuplot %s%s", (imp)? "-m " : "", edttext);
    if (sr->code == GR_PLOT) { 
	strcat(line, " time");
    }

    if (check_cmd(line) || cmd_init(line)) return;
    lines = mymalloc(command.list[0] - 1);
    if (lines == NULL) return;
    for (i=0; i<command.list[0]-1 ; i++) {
	if (sr->code == GR_PLOT) lines[i] = 1;
	else lines[i] = 0;
    }

    if (imp) {
	err = gnuplot(command.list, NULL, &Z, datainfo,
		      &paths, &plot_count, 0, 1, OPT_M);
    } else {
	err = gnuplot(command.list, lines, &Z, datainfo,
		      &paths, &plot_count, 0, 1, 0);
    }
    if (err == -999)
	errbox(_("No data were available to graph"));
    else if (err < 0) errbox(_("gnuplot command failed"));
    else register_graph();

    free(lines);
}

/* ........................................................... */

void display_var (void)
{
    int list[2];
    PRN *prn;

    list[0] = 1;
    list[1] = mdata->active_var;
    if (bufopen(&prn)) return;
    printdata(list, &Z, datainfo, 0, 1, prn);
    view_buffer(prn, 24, 350, _("gretl: display data"), PRINT, NULL);    
}

/* ........................................................... */

void do_run_script (gpointer data, guint code, GtkWidget *w)
{
    PRN *prn;
    char *runfile = NULL, fname[MAXLEN];
    int err;

#ifdef notdef
    fprintf(stderr, "do_run_script(): data=%p, code=%d, w=%p\n",
	    (void *) data, code, (void *) w);
#endif

    if (!user_fopen("gretl_output_tmp", fname, &prn)) return;

    if (code == SCRIPT_EXEC) runfile = scriptfile;
    else if (code == SESSION_EXEC) runfile = cmdfile;

    if (data != NULL) { /* get commands from file view buffer */
	GdkCursor *cursor;
	windata_t *mydata = (windata_t *) data;
	gchar *buf = gtk_editable_get_chars(GTK_EDITABLE(mydata->w), 0, -1);

	if (buf == NULL || !strlen(buf)) {
	    errbox("No commands to execute");
	    gretl_print_destroy(prn);
	    if (buf) g_free(buf);
	    return;
	}

	cursor = gdk_cursor_new(GDK_WATCH);
	gdk_window_set_cursor(mydata->dialog->window, cursor);
	gdk_window_set_cursor(GTK_TEXT(mydata->w)->text_area, cursor);
	gdk_cursor_destroy(cursor);	

	err = execute_script(NULL, buf, NULL, NULL, prn, code);
	g_free(buf);

	gdk_window_set_cursor(mydata->dialog->window, NULL);
	gdk_window_set_cursor(GTK_TEXT(mydata->w)->text_area, NULL);
    } else
	err = execute_script(runfile, NULL, NULL, NULL, prn, code);

    gretl_print_destroy(prn);

    if (err == -1) return;

    refresh_data();

    view_file(fname, 1, 1, 77, 450, SCRIPT_OUT, script_out_items);
}

/* ........................................................... */

void do_open_script (GtkWidget *w, GtkFileSelection *fs)
{
    int ret, n = strlen(paths.scriptdir);

    if (fs) {
	if (isdir(gtk_file_selection_get_filename(GTK_FILE_SELECTION(fs))))
	    return;
	strncpy(scriptfile, 
		gtk_file_selection_get_filename(GTK_FILE_SELECTION(fs)), 
		MAXLEN-1);
	gtk_widget_destroy(GTK_WIDGET (fs)); 
    } else {
	strcpy(scriptfile, tryscript); /* might cause problems? */
    }

    /* is this a "session" file? */
    ret = saved_objects(scriptfile);
    if (ret == -1) {
	sprintf(errtext, _("Couldn't open %s"), tryscript);
	errbox(errtext);
	delete_from_filelist(2, tryscript);
	delete_from_filelist(3, tryscript);
	return;
    }
    else if (ret > 0) {
	verify_open_session(NULL);
	return;
    }

    /* or just an "ordinary" script */
    mkfilelist(3, scriptfile);

    if (strncmp(scriptfile, paths.scriptdir, n)) 
	view_file(scriptfile, 1, 0, 78, 370, EDIT_SCRIPT, script_items);
    else 
	view_file(scriptfile, 0, 0, 78, 370, VIEW_SCRIPT, sample_script_items);
}

/* ........................................................... */

void do_new_script (gpointer data, guint loop, GtkWidget *widget) 
{
    PRN *prn;
    char fname[MAXLEN];

    if (!user_fopen("script_tmp", fname, &prn)) return;
    if (loop) pprintf(prn, "loop 1000\n\nendloop\n");
    gretl_print_destroy(prn);
    strcpy(scriptfile, fname);
    
    view_file(scriptfile, 1, 0, 77, 350, EDIT_SCRIPT, script_items);
}

/* ........................................................... */

void do_open_csv_box (char *fname, int code, int append)
{
    int err;
    PRN *prn;
    char buf[30];

    if (bufopen(&prn)) return;

    if (code == OPEN_BOX)
	err = import_box(&Z, datainfo, fname, prn);
    else
	err = import_csv(&Z, datainfo, fname, prn); 

    sprintf(buf, _("gretl: import %s data"), 
	    (code == OPEN_BOX)? "BOX" : "CSV");
    view_buffer(prn, 77, 350, buf, IMPORT, NULL); 

    if (err) return;

    data_status |= IMPORT_DATA;
    strcpy(paths.datfile, fname);

    register_data(fname, !append);
}

/* ........................................................... */

int do_store (char *mydatfile, const int opt, int overwrite)
{
    char f = getflag(opt);
    gchar *msg;
    FILE *fp;
    int showlist = 1;

    line[0] = '\0';

    if (storelist == NULL) showlist = 0;

    if (f) { /* not a standard native save */
	sprintf(line, "store '%s' %s -%c", mydatfile, 
		(showlist)? storelist : "", f);
    } else {
	if (!overwrite) {
	    fp = fopen(mydatfile, "r");
	    if (fp != NULL) {
		fclose(fp);
		if (yes_no_dialog(_("gretl: save data"), 
				  _("There is already a data file of this name.\n"
				  "OK to overwrite it?"), 0)) {
		    return 1;
		}
	    }
	}
	sprintf(line, "store '%s' %s", mydatfile, 
	       (showlist)? storelist : "");   
	strcpy(paths.datfile, mydatfile);
    }

    if (check_cmd(line) || cmd_init(line)) return 1; 

    /* back up existing datafile if need be */
    if ((fp = fopen(mydatfile, "r")) && fgetc(fp) != EOF &&
	fclose(fp) == 0) {
	char backup[MAXLEN];

	sprintf(backup, "%s~", mydatfile);
	if (copyfile(mydatfile, backup)) {
	    errbox(_("Couldn't make backup of data file"));
	    return 1;
	}
    }

    if (write_data(mydatfile, command.list, Z, datainfo, 
		   data_option(opt), &paths)) {
	sprintf(errtext, _("Write of data file failed\n%s"),
		get_gretl_errmsg());
	errbox(errtext);
	return 1;
    }

    if (opt != OPT_M && opt != OPT_R && opt != OPT_R_ALT)
	mkfilelist(1, mydatfile);

    msg = g_strdup_printf(_("%s written OK"), mydatfile);
    infobox(msg);
    g_free(msg);

    /* record that data have been saved */
    if (!f) {
	data_status = (HAVE_DATA|USER_DATA);
	set_sample_label(datainfo);
    }

    return 0;
}

#ifdef G_OS_WIN32

static void win_show_error (void)
{
    LPVOID buf;

    FormatMessage( 
		  FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		  FORMAT_MESSAGE_FROM_SYSTEM | 
		  FORMAT_MESSAGE_IGNORE_INSERTS,
		  NULL,
		  GetLastError(),
		  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		  (LPTSTR) &buf,
		  0,
		  NULL 
		  );
    MessageBox(NULL, (LPCTSTR) buf, "Error", MB_OK | MB_ICONERROR);
    LocalFree(buf);
}

static int get_latex_path (char *latex_path)
{
    int ret;
    char *p;

    ret = SearchPath(NULL, "latex.exe", NULL, MAXLEN, latex_path, &p);

    return (ret == 0);
}

int winfork (char *cmdline, const char *dir, int wshow)
{
    int child;
    STARTUPINFO si;
    PROCESS_INFORMATION pi; 

    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    si.dwFlags = STARTF_USESHOWWINDOW;
    si.wShowWindow = wshow;

    ZeroMemory(&pi, sizeof(pi));    

    /* zero return means failure */
    child = CreateProcess(NULL, cmdline, 
			  NULL,NULL, FALSE,
			  CREATE_NEW_CONSOLE | HIGH_PRIORITY_CLASS,
			  NULL, dir,
			  &si, &pi);

    if (!child) {
	win_show_error();
	return 1;
    }

    WaitForSingleObject(pi.hProcess, INFINITE);    
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);

    return 0;
}

#endif

/* ........................................................... */

void view_latex (gpointer data, guint prn_code, GtkWidget *widget)
{
    char texfile[MAXLEN], texbase[MAXLEN], tmp[MAXLEN];
    int dot, err;
    windata_t *mydata = (windata_t *) data;
    MODEL *pmod = (MODEL *) mydata->data;

    if (pmod->errcode == E_NAN) {
	errbox("Sorry, can't format this model");
	return;
    }

    if (prn_code)
	err = eqnprint(pmod, datainfo, &paths, texfile, model_count, 1);
    else 
	err = tabprint(pmod, datainfo, &paths, texfile, model_count, 1);
	
    if (err) {
	errbox(_("Couldn't open tex file for writing"));
	return;
    }

    dot = dotpos(texfile);
    *texbase = 0;
    strncat(texbase, texfile, dot);     

#ifdef G_OS_WIN32
    {
	static char latex_path[MAXLEN];
	char *texshort = strrchr(texbase, SLASH) + 1;

	if (*latex_path == 0 && get_latex_path(latex_path)) {
	    win_show_error();
	    return;
	}

	sprintf(tmp, "\"%s\" %s", latex_path, texshort);
	if (winfork(tmp, paths.userdir, SW_SHOWMINIMIZED)) {
	    return;
	} else {
	    sprintf(tmp, "\"%s\" %s.dvi", viewdvi, texbase);
	    if (WinExec(tmp, SW_SHOWNORMAL) < 32)
		win_show_error();	
	}
    }
#else
    sprintf(tmp, "cd %s && latex %s", paths.userdir, texbase);
    err = system(tmp);
    if (err) 
	errbox("Failed to run latex");
    else 
	gretl_fork(viewdvi, texbase);
#endif

    remove(texfile);
#ifndef G_OS_WIN32
    sleep(2); /* let forked xdvi get the DVI file */
    sprintf(tmp, "%s.dvi", texbase);
    remove(tmp);
#endif
    sprintf(tmp, "%s.log", texbase);
    remove(tmp);
    sprintf(tmp, "%s.aux", texbase);
    remove(tmp);
}

/* ........................................................... */

void do_save_tex (char *fname, const int code, MODEL *pmod)
{
    PRN *texprn;

    texprn = gretl_print_new(GRETL_PRINT_FILE, fname);
    if (texprn == NULL) {
	errbox(_("Couldn't open tex file for writing"));
	return;
    }  

    if (code == SAVE_TEX_EQ)
	tex_print_equation(pmod, datainfo, 1, texprn);
    else if (code == SAVE_TEX_TAB)
	tex_print_model(pmod, datainfo, 1, texprn);
    else if (code == SAVE_TEX_EQ_FRAG)
	tex_print_equation(pmod, datainfo, 0, texprn);
    else if (code == SAVE_TEX_TAB_FRAG)
	tex_print_model(pmod, datainfo, 0, texprn);

    gretl_print_destroy(texprn);

    infobox(_("LaTeX file saved"));
}

/* ........................................................... */

static char *bufgets (char *s, int size, const char *buf)
{
    int i;
    static const char *p;

    /* mechanism for resetting p */
    if (s == NULL || size == 0) {
	p = NULL;
	return 0;
    }

    /* start at beginning of buffer */
    if (p == NULL) p = buf;

    /* signal that we've reached the end of the buffer */
    if (p && *p == 0) return NULL;

    *s = 0;
    /* advance to newline, end of buffer, or maximum size,
       whichever comes first */
    for (i=0; i<size; i++) {
	s[i] = p[i];
	if (p[i] == 0) break;
	if (p[i] == '\n') {
	    /* throw away newlines */
	    s[i] = 0;
	    break;
	}
    }
    /* advance the buffer pointer */
    p += i + (p[i] != 0);
    return s;
}

/* ........................................................... */

int execute_script (const char *runfile, const char *buf,
		    SESSION *psession, SESSIONBUILD *rebuild,
		    PRN *prn, int exec_code)
     /* run commands from runfile or buf, output to prn */
{
    FILE *fb = NULL;
    int cont, exec_err = 0;
    int i, j = 0, loopstack = 0, looprun = 0;
    char tmp[MAXLEN];
    LOOPSET loop;            /* struct for monte carlo loop */

    if (runfile != NULL) { /* we'll get commands from file */
	fb = fopen(runfile, "r");
	if (fb == NULL) {
	    errbox(_("Couldn't open script"));
	    return -1;
	}

	/* check that the file has something in it */
	cont = 0;
	while (fgets(tmp, MAXLEN-1, fb)) {
	    if (strlen(tmp)) {
		for (i=0; i<strlen(tmp); i++) {
		    if (!isspace(tmp[i])) {
			cont = 1;
			break;
		    }
		}
	    }
	    if (cont) break;
	}
	fclose(fb);

	if (!cont) {
	    errbox(_("No commands to execute"));
	    return -1;
	}
    } else { /* no runfile, commands from buffer */
	if (buf == NULL || !strlen(buf)) {
	    errbox(_("No commands to execute"));
	    return -1;	
	}
    }

    if (runfile != NULL) fb = fopen(runfile, "r");
    else bufgets(NULL, 0, buf);

    /* reset model count to 0 if starting/saving session */
    if (exec_code == SESSION_EXEC || exec_code == REBUILD_EXEC ||
	exec_code == SAVE_SESSION_EXEC) 
	model_count = 0;

    /* monte carlo struct */
    loop.lines = NULL;
    loop.models = NULL;
    loop.lmodels = NULL;
    loop.prns = NULL;
    loop.storename = NULL;
    loop.storelbl = NULL;
    loop.storeval = NULL;
    loop.nmod = 0;

    /* Put the action of running this script into the command log? */
    if (exec_code == SCRIPT_EXEC && runfile != NULL) {
	char runcmd[MAXLEN];

	sprintf(runcmd, "run %s", runfile);
	check_cmd(runcmd);
	cmd_init(runcmd);
    }

    command.cmd[0] = '\0';

    while (strcmp(command.cmd, "quit")) {
	if (looprun) { /* Are we doing a Monte Carlo simulation? */
	    if (!loop.ncmds) {
		pprintf(prn, _("No commands in loop\n"));
		looprun = 0;
		continue;
	    }
	    i = 0;
	    while (j != MAXLOOP && loop_condition(i, &loop, Z, datainfo)) {
		if (loop.type == FOR_LOOP && !echo_off)
		    pprintf(prn, "loop: i = %d\n\n", i + 1);
		for (j=0; j<loop.ncmds; j++) {
		    if (loop_exec_line(&loop, i, j, prn)) {
			pprintf(prn, _("Error in command loop: aborting\n"));
			j = MAXLOOP - 1;
			i = loop.ntimes;
		    }
		}
		i++;
	    }
	    if (j != MAXLOOP) 
		print_loop_results(&loop, datainfo, prn, &paths, 
				   &model_count, loopstorefile);
	    looprun = 0;
	    monte_carlo_free(&loop);
	    if (j == MAXLOOP) return 1;
	} else { /* end if Monte Carlo stuff */
	    line[0] = '\0';
	    if ((fb && fgets(line, MAXLEN, fb) == NULL) ||
		(fb == NULL && bufgets(line, MAXLEN, buf) == NULL)) 
		goto endwhile;
	    while ((cont = top_n_tail(line))) {
		if (cont == E_ALLOC) {
		    errbox(_("Out of memory loading command line"));
		    return 1;
		}
		*tmp = '\0';
		if (fb)
		    fgets(tmp, MAXLEN-1, fb);
		else
		    bufgets(tmp, MAXLEN-1, buf); 
		strcat(line, tmp);
		compress_spaces(line);
	    }
	    if (!strncmp(line, "noecho", 6)) echo_off = 1;
	    if (strncmp(line, "(* saved objects:", 17) == 0) 
		strcpy(line, "quit"); 
	    else if (!echo_off) {
		if ((line[0] == '(' && line[1] == '*') ||
		    (line[strlen(line)-1] == ')' && 
		     line[strlen(line)-2] == '*')) 
		    pprintf(prn, "\n%s\n", line);
		else 
		    pprintf(prn, "\n? %s\n", line);	
	    }
	    oflag = 0;
	    strcpy(tmp, line);
	    exec_err = gui_exec_line(line, &loop, &loopstack, 
				     &looprun, psession, rebuild, 
				     prn, exec_code, runfile);
	    if (exec_err) {
		pprintf(prn, _("\nError executing script: halting\n"));
		pprintf(prn, "> %s\n", tmp);
		return 1;
	    }
	} /* end alternative to Monte Carlo stuff */
    } /* end while() */
 endwhile:
    if (fb) fclose(fb);
    if (psession && rebuild) /* recreating a gretl session */
	save_model_copy(&models[0], psession, rebuild, datainfo);
    return 0;
}

/* ........................................................... */

static int ready_for_command (char *line)
{
    if (*line == 'q' || *line == 'x' ||
        *line == '\0' ||
        strncmp(line, "open", 4) == 0 ||
        strncmp(line, "run", 3) == 0 ||
        strncmp(line, "nulldata", 6) == 0 ||
        strncmp(line, "import", 4) == 0 ||
        strncmp(line, "pvalue", 6) == 0 ||
        strncmp(line, "!", 1) == 0 ||
        strncmp(line, "(*", 2) == 0 ||
        strncmp(line, "man ", 4) == 0 ||
        strncmp(line, "help", 4) == 0 ||
	strncmp(line, "noecho", 6) == 0 ||
        strncmp(line, "critical", 8) == 0)
        return 1;
    return 0;
}

/* ........................................................... */

static int script_model_test (const int id, PRN *prn, const int ols_only)
{
    /* need to work in terms of modelspec here, _not_ model_count */

    int m = (id)? id - 1 : 0;

    if (model_count == 0) { 
	pprintf(prn, _("Can't do this: no model has been estimated yet\n"));
	return 1;
    }
    if (id > model_count) { 
	pprintf(prn, _("Can't do this: there is no model %d\n"), id);
	return 1;
    }

    /* ID == 0 -> no model specified -> look for last script model */
    if (modelspec != NULL && id == 0) {
	m = model_count - 1;
	while (m) { 
	    if (model_origin[m] == 's') break;
	    m--;
	}
    }
    if (modelspec == NULL || model_origin[m] == 'g') {
	pprintf(prn, _("Sorry, can't do this.\nTo operate on a model estimated "
		"via the graphical interface, please use the\nmenu items in "
		"the model window.\n"));
	return 1;
    }    
    if (ols_only && strncmp(modelspec[m].cmd, "ols", 3)) {
	pprintf(prn, _("This command is only available for OLS models "
		"at present\n"));
	return 1;
    }
    if (model_sample_issue(NULL, &modelspec[m], (fullZ == NULL)? Z : fullZ, 
			   (fullZ == NULL)? datainfo : fullinfo)) {
	pprintf(prn, _("Can't do: the current data set is different from "
		"the one on which\nthe reference model was estimated\n"));
	return 1;
    }
    return 0;
}

/* ........................................................... */

static int gui_exec_line (char *line, 
			  LOOPSET *plp, int *plstack, int *plrun, 
			  SESSION *psession, SESSIONBUILD *rebuild,
			  PRN *prn, int exec_code, 
			  const char *myname) 
{
    int i, err = 0, check = 0, rebuilding = 0, order, nulldata_n, lines[1];
    double rho;
    char runfile[MAXLEN], datfile[MAXLEN];
    char linecopy[1024];
    char texfile[MAXLEN];
    MODEL tmpmod;
    FREQDIST *freq;             /* struct for freq distributions */
    GRETLTEST test;             /* struct for model tests */
    GRETLTEST *ptest;

    if (!data_status && !ready_for_command(line)) {
	pprintf(prn, _("You must open a data file first\n"));
	return 1;
    }

    if (psession != NULL && rebuild != NULL) rebuilding = 1;

#ifdef CMD_DEBUG
    fprintf(stderr, "gui_exec_line: '%s'\n", line);
#endif

    /* parse the command line */
    strncpy(linecopy, line, 1023);
    linecopy[1023] = '\0';
    catchflag(line, &oflag);
    /* but if we're stacking commands for a loop, parse "lightly" */
    if (*plstack) 
	get_cmd_ci(line, &command);
    else 
	getcmd(line, datainfo, &command, &ignore, &Z, cmds);
    if (command.ci == -2) { /* line was a comment, pass */
#ifdef notdef
 	cmds->fp = fopen(cmdfile, "a");
 	if (cmds->fp) {
 	    pprintf(cmds, "%s\n", linecopy);
 	    fclose(cmds->fp);
 	}
#endif
	return 0;
    }
    if (command.ci < 0) return 0; /* nothing there */
    if (command.errcode) {
        errmsg(command.errcode, prn);
	if (exec_code == CONSOLE_EXEC) {
	    strcpy(errline, linecopy);
	}
        return 1;
    }
    if (*plstack) {  /* accumulating loop commands */
	if (!ok_in_loop(command.ci)) {
            pprintf(prn, _("Sorry, this command is not available in loop mode\n"));
            return 1;
        } else {
            if (!echo_off) echo_cmd(&command, datainfo, line, 1, 1, oflag, cmds);
            if (command.ci != ENDLOOP) {
                if (add_to_loop(plp, line, command.ci, oflag)) {
                    pprintf(prn, _("Failed to add command to loop stack\n"));
		    return 1;
                }
                return 0;
            } 
        }
    } 

    /* if rebuilding a session, add tests back to models */
    if (rebuild) ptest = &test;
    else ptest = NULL;

    /* if rebuilding a session, put the commands onto the stack */
    if (rebuild) cmd_init(line);

#ifdef notdef /* ??? */
     if (is_model_ref_cmd(command.ci)) {
 	if (model_sample_issue(models[0], &Z, datainfo)) {
 	    pprintf(prn, _("Can't do: the current data set is different from "
 		   "the one on which\nthe reference model was estimated\n"));
 	    return 1;
 	}
     }
#endif

    switch (command.ci) {

    case ADF: case COINT:
    case CORR:
    case CRITERIA: case CRITICAL:
    case DIFF: case LDIFF: case LAGS: case LOGS:
    case MULTIPLY:
    case GRAPH: case PLOT:
    case INFO: case LABELS: case VARLIST:
    case PRINT:
    case SUMMARY:
    case MEANTEST: case VARTEST:
    case RUNS: case SPEARMAN:
	err = simple_commands(&command, line, &Z, datainfo, &paths,
			      0, oflag, prn);
	break;

    case ADD:
    case OMIT:
	if ((err = script_model_test(0, prn, 0))) break;
    plain_add_omit:
	clear_model(models[1], NULL);
	if (command.ci == ADD || command.ci == ADDTO)
	    err = auxreg(command.list, models[0], models[1], &model_count, 
			 &Z, datainfo, AUX_ADD, prn, NULL);
	else
	    err = omit_test(command.list, models[0], models[1],
			    &model_count, &Z, datainfo, prn);
	if (err) {
	    errmsg(err, prn);
	    clear_model(models[1], NULL);
	} else {
	    /* for command-line use, we keep a stack of 
	       two models, and recycle the places */
	    swap_models(&models[0], &models[1]);
	    clear_model(models[1], NULL);
	    if (oflag) outcovmx(models[0], datainfo, 0, prn);
	}
	break;	

    case ADDTO:
    case OMITFROM:
	i = atoi(command.param);
	if ((err = script_model_test(i, prn, 0))) break;
	if (i == (models[0])->ID) goto plain_add_omit;
	err = re_estimate(modelspec[i-1].cmd, &tmpmod, &Z, datainfo);
	if (err) {
	    pprintf(prn, _("Failed to reconstruct model %d\n"), i);
	    break;
	} 
	clear_model(models[1], NULL);
	tmpmod.ID = i;
	if (command.ci == ADDTO)
	    err = auxreg(command.list, &tmpmod, models[1], &model_count, 
			 &Z, datainfo, AUX_ADD, prn, NULL);
	else
	    err = omit_test(command.list, &tmpmod, models[1],
			    &model_count, &Z, datainfo, prn);
	if (err) {
	    errmsg(err, prn);
	    clear_model(models[1], NULL);
	    break;
	} else {
	    swap_models(&models[0], &models[1]);
	    clear_model(models[1], NULL);
	    if (oflag) outcovmx(models[0], datainfo, 0, prn);
	}
	clear_model(&tmpmod, NULL);
	break;

    case AR:
	if (rebuilding)
	    save_model_copy(&models[0], psession, rebuild, datainfo);
	else
	    clear_model(models[0], datainfo); 
	*models[0] = ar_func(command.list, atoi(command.param), &Z, 
			    datainfo, &model_count, prn);
	if ((err = (models[0])->errcode)) { 
	    errmsg(err, prn); 
	    break;
	}
	if (oflag) outcovmx(models[0], datainfo, 0, prn);
	break;

    case ARCH:
	order = atoi(command.param);
	clear_model(models[1], NULL);
	*models[1] = arch(order, command.list, &Z, datainfo, 
			  &model_count, prn, ptest);
	if ((err = (models[1])->errcode)) 
	    errmsg(err, prn);
	if ((models[1])->ci == ARCH) {
	    swap_models(&models[0], &models[1]);
	    if (oflag) outcovmx(models[0], datainfo, 0, prn);
	} else if (rebuild)
	    add_test_to_model(ptest, models[0]);
	clear_model(models[1], NULL);
	break;

    case BXPLOT:
	if (exec_code == REBUILD_EXEC || exec_code == SAVE_SESSION_EXEC) 
	    break;
	if (command.nolist) 
	    err = boolean_boxplots(line, &Z, datainfo, (oflag != 0));
	else
	    err = boxplots(command.list, NULL, &Z, datainfo, (oflag != 0));
	break;

    case CHOW:
	if ((err = script_model_test(0, prn, 1))) break;
	err = chow_test(line, models[0], &Z, datainfo, prn, ptest);
	if (err) errmsg(err, prn);
	else if (rebuild) 
	    add_test_to_model(ptest, models[0]);
	break;

    case CUSUM:
	if ((err = script_model_test(0, prn, 1))) break;
	err = cusum_test(models[0], &Z, datainfo, prn, 
			 &paths, ptest);
	if (err) errmsg(err, prn);
	else if (rebuild) 
	    add_test_to_model(ptest, models[0]);
	break;

    case CORC:
    case HILU:
	err = hilu_corc(&rho, command.list, &Z, datainfo, command.ci, prn);
	if (err) {
	    errmsg(err, prn);
	    break;
	}
	if (rebuilding)
	    save_model_copy(&models[0], psession, rebuild, datainfo);
	else
	    clear_model(models[0], datainfo); 
	*models[0] = lsq(command.list, &Z, datainfo, command.ci, 1, rho);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    break;
	}
	++model_count;
	(models[0])->ID = model_count;
	if (printmodel(models[0], datainfo, prn))
	    (models[0])->errcode = E_NAN;
	if (oflag) outcovmx(models[0], datainfo, 0, prn);
	break;

    case CORRGM:
	order = atoi(command.param);
	err = corrgram(command.list[1], order, &Z, datainfo, &paths,
		       1, prn);
	if (err) pprintf(prn, _("Failed to generate correlogram\n"));
	break;

    case DELEET:
	if (fullZ != NULL) {
	    pprintf(prn, _("Can't delete last variable when in sub-sample"
		    " mode\n"));
	    break;
	}
	if (datainfo->v <= 1 || dataset_drop_vars(1, &Z, datainfo)) 
	    pprintf(prn, _("Failed to shrink the data set"));
	else varlist(datainfo, prn);
	break;

    case ENDLOOP:
	if (*plstack != 1) {
	    pprintf(prn, _("You can't end a loop here, "
		    "you haven't started one\n"));
	    break;
	}
	*plstack = 0;
	*plrun = 1;
	break;

    case EQNPRINT:
    case TABPRINT:
	if ((models[0])->errcode == E_NAN) {
	    pprintf(prn, _("Couldn't format model\n"));
	    break;
	}
	if ((err = script_model_test(0, prn, 1))) break;
	if (command.ci == EQNPRINT)
	    err = eqnprint(models[0], datainfo, &paths, 
			   texfile, model_count, oflag);
	else
	    err = tabprint(models[0], datainfo, &paths, 
			   texfile, model_count, oflag);
	if (err) 
	    pprintf(prn, _("Couldn't open tex file for writing\n"));
	else 
	    pprintf(prn, _("Model printed to %s\n"), texfile);
	break;

    case FCAST:
	if ((err = script_model_test(0, prn, 0))) break;
	err = fcast(line, models[0], datainfo, &Z);
	if (err < 0) {
	    err *= -1;
	    printf(_("Error retrieving fitted values\n"));
	    errmsg(err, prn);
	    break;
	}
	err = 0;
	varlist(datainfo, prn);
	break;

    case FCASTERR:
	if ((err = script_model_test(0, prn, 0))) break;
	err = fcast_with_errs(line, models[0], &Z, datainfo, prn,
			      &paths, oflag); 
	if (err) errmsg(err, prn);
	break;

    case FIT:
	if ((err = script_model_test(0, prn, 0))) break;
	err = fcast("fcast autofit", models[0], datainfo, &Z);
	if (err < 0) {
	    err *= -1;
	    errmsg(err, prn);
	    break;
	}
	err = 0;
	pprintf(prn, _("Retrieved fitted values as \"autofit\"\n"));
	varlist(datainfo, prn); 
	if (dataset_is_time_series(datainfo)) {
	    plotvar(&Z, datainfo, "time");
	    command.list = myrealloc(command.list, 4 * sizeof(int));
	    command.list[0] = 3; 
	    command.list[1] = (models[0])->list[1];
	    command.list[2] = varindex(datainfo, "autofit");
	    command.list[3] = varindex(datainfo, "time");
	    lines[0] = oflag;
	    err = gnuplot(command.list, lines, &Z, datainfo,
			  &paths, &plot_count, 1, 0, 0);
	    if (err < 0) pprintf(prn, _("gnuplot command failed\n"));
	    else register_graph();
	}
	break;
		
    case FREQ:
	freq = freqdist(&Z, datainfo, command.list[1], 1);
	if ((err = freq_error(freq, prn))) {
	    break;
	}
	printfreq(freq, prn);
	if (exec_code == CONSOLE_EXEC) {
	    if (plot_freq(freq, &paths, NORMAL))
		pprintf(prn, _("gnuplot command failed\n"));
	}
	free_freq(freq);
	break;

    case GENR:
	err = generate(&Z, datainfo, line, model_count,
		       (last_model == 's')? models[0] : models[2], 
		       oflag);
	if (err) 
	    errmsg(err, prn);
	else {
	    pprintf(prn, "%s\n", get_gretl_msg()); 
	    if (exec_code == CONSOLE_EXEC)
		populate_clist(mdata->listbox, datainfo);
	}
	break;

    case GNUPLOT:
	if (exec_code == SAVE_SESSION_EXEC || exec_code == REBUILD_EXEC)
	    break;
	if (plp != NULL) {
	    pprintf(prn, _("script mode: gnuplot command ignored\n"));
	    break;
	}
	if (oflag == OPT_M) { /* plot with impulses */
	    err = gnuplot(command.list, NULL, &Z, datainfo,
			  &paths, &plot_count, 1, 0, OPT_M);
	} else {	
	    lines[0] = oflag;
	    err = gnuplot(command.list, lines, &Z, datainfo,
			  &paths, &plot_count, 0, 0, 0);
	}
	if (err < 0) pprintf(prn, _("gnuplot command failed\n"));
	else register_graph();
	break;

    case HAUSMAN:
	if ((err = script_model_test(0, prn, 0))) break;
	err = hausman_test(models[0], &Z, datainfo, prn);
	break;

    case HCCM:
    case HSK:
	if (rebuilding)
	    save_model_copy(&models[0], psession, rebuild, datainfo);
	else
	    clear_model(models[0], datainfo); 
	if (command.ci == HCCM)
	    *models[0] = hccm_func(command.list, &Z, datainfo);
	else
	    *models[0] = hsk_func(command.list, &Z, datainfo);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    break;
	}
	++model_count;
	(models[0])->ID = model_count;
	if (printmodel(models[0], datainfo, prn))
	    (models[0])->errcode = E_NAN;
	if (!oflag) break;
	if (command.ci == HCCM) 
	    print_white_vcv(models[0], prn);
	else
	    outcovmx(models[0], datainfo, 0, prn);
	break;

    case HELP:
	if (strlen(command.param)) 
	    help(command.param, paths.cmd_helpfile, prn);
	else help(NULL, paths.cmd_helpfile, prn);
	break;

    case IMPORT:
	if (exec_code == SAVE_SESSION_EXEC) break;
        err = getopenfile(line, datfile, &paths, 0, 0);
        if (err) {
            pprintf(prn, _("import command is malformed\n"));
            break;
        }
	if (data_status & HAVE_DATA)
	    close_session();
        if (oflag)
            err = import_box(&Z, datainfo, datfile, prn);
        else
            err = import_csv(&Z, datainfo, datfile, prn);
        if (!err) { 
	    data_status |= IMPORT_DATA;
	    register_data(datfile, (exec_code != REBUILD_EXEC));
            print_smpl(datainfo, 0, prn);
            varlist(datainfo, prn);
            pprintf(prn, _("You should now use the \"print\" command "
		    "to verify the data\n"));
            pprintf(prn, _("If they are OK, use the  \"store\" command "
		    "to save them in gretl format\n"));
        }
        break;

    case OPEN:
	if (exec_code == SAVE_SESSION_EXEC) break;
	err = getopenfile(line, datfile, &paths, 0, 0);
	if (err) {
	    errbox(_("'open' command is malformed"));
	    break;
	}
#ifdef CMD_DEBUG
	fprintf(stderr, "OPEN in gui_exec_line, datfile='%s'\n", datfile);
#endif
	if (data_status & HAVE_DATA)
	    close_session();
	check = detect_filetype(datfile, &paths, prn);
	if (check == GRETL_CSV_DATA)
	    err = import_csv(&Z, datainfo, datfile, prn);
	else if (check == GRETL_BOX_DATA)
	    err = import_box(&Z, datainfo, datfile, prn);
	else if (check == GRETL_XML_DATA)
	    err = get_xmldata(&Z, datainfo, datfile, &paths, data_status, prn, 0);
	else
	    err = get_data(&Z, datainfo, datfile, &paths, data_status, prn);
	if (err) {
	    gui_errmsg(err);
	    break;
	}
	strncpy(paths.datfile, datfile, MAXLEN-1);
	if (check == GRETL_CSV_DATA || check == GRETL_BOX_DATA)
	    data_status |= IMPORT_DATA;
	register_data(paths.datfile, (exec_code != REBUILD_EXEC));
	varlist(datainfo, prn);
	paths.currdir[0] = '\0'; 
	break;

    case LMTEST:
	if ((err = script_model_test(0, prn, 1))) break;
	/* non-linearity (squares) */
	if (oflag == OPT_S || oflag == OPT_O || !oflag) {
	    err = auxreg(NULL, models[0], models[1], &model_count, 
			 &Z, datainfo, AUX_SQ, prn, ptest);
	    clear_model(models[1], NULL);
	    model_count--;
	    if (err) errmsg(err, prn);
	}
	/* non-linearity (logs) */
	if (oflag == OPT_L || oflag == OPT_O || !oflag) {
	    err = auxreg(NULL, models[0], models[1], &model_count, 
			 &Z, datainfo, AUX_LOG, prn, ptest);
	    clear_model(models[1], NULL);
	    model_count--;
	    if (err) errmsg(err, prn);
	}
	/* autocorrelation or heteroskedasticity */
	if (oflag == OPT_M || oflag == OPT_O) {
	    int order = atoi(command.param);

	    err = autocorr_test(models[0], order, &Z, datainfo, prn, ptest);
	    if (err) errmsg(err, prn);
	    /* FIXME: need to respond? */
	} 
	if (oflag == OPT_C || !oflag) {
	    err = whites_test(models[0], &Z, datainfo, prn, ptest);
	    if (err) errmsg(err, prn);
	}
	if (rebuild)
	    add_test_to_model(ptest, models[0]);
	break;

    case LOGIT:
    case PROBIT:
	if (rebuilding)
	    save_model_copy(&models[0], psession, rebuild, datainfo);
	else
	    clear_model(models[0], datainfo); 
	*models[0] = logit_probit(command.list, &Z, datainfo, command.ci);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn);
	    break;
	}
	++model_count;
	(models[0])->ID = model_count;
	if (printmodel(models[0], datainfo, prn))
	    (models[0])->errcode = E_NAN;
	if (oflag) outcovmx(models[0], datainfo, 0, prn); 
	break;

    case LOOP:
	if (plp == NULL) {
	    pprintf(prn, _("Sorry, Monte Carlo loops not available "
		    "in this mode\n"));
	    break;
	}
	if ((err = parse_loopline(line, plp, datainfo))) {
	    pprintf(prn, "%s\n", get_gretl_errmsg());
	    break;
	}
	if (plp->lvar == 0 && plp->ntimes < 2) {
	    printf(_("Loop count missing or invalid\n"));
	    monte_carlo_free(plp);
	    break;
	}
	*plstack = 1; 
	break;

    case NOECHO:
	echo_off = 1;
	break;

    case NULLDATA:
	nulldata_n = atoi(command.param);
	if (nulldata_n < 2) {
	    pprintf(prn, _("Data series length count missing or invalid\n"));
	    err = 1;
	    break;
	}
	if (nulldata_n > 1000000) {
	    pprintf(prn, _("Data series too long\n"));
	    err = 1;
	    break;
	}
	err = open_nulldata(&Z, datainfo, data_status, nulldata_n, prn);
	if (err) { 
	    pprintf(prn, _("Failed to create empty data set\n"));
	    break;
	}
	paths.datfile[0] = '\0';
	populate_clist(mdata->listbox, datainfo);
	data_status = HAVE_DATA | GUI_DATA | MODIFIED_DATA;
	set_sample_label(datainfo);
	orig_vars = datainfo->v;
	menubar_state(TRUE);
	break;

    case OLS:
    case WLS:
    case POOLED:
	if (rebuilding)
	    save_model_copy(&models[0], psession, rebuild, datainfo);
	else
	    clear_model(models[0], datainfo); 
	*models[0] = lsq(command.list, &Z, datainfo, command.ci, 1, 0.0);
	if ((err = (models[0])->errcode)) {
	    errmsg(err, prn); 
	    break;
	}
	++model_count;
	(models[0])->ID = model_count;
	if (printmodel(models[0], datainfo, prn))
	    (models[0])->errcode = E_NAN;
	if (oflag) outcovmx(models[0], datainfo, 0, prn); 
	break;

#ifdef ENABLE_GMP
    case MPOLS:
	err = mp_ols(command.list, command.param, &Z, datainfo, prn);
	break;
#endif

    case PANEL:
	err = set_panel_structure(oflag, datainfo, prn);
	break;

    case PERGM:
	err = periodogram(command.list[1], &Z, datainfo, &paths,
			  1, oflag, prn);
	if (err) pprintf(prn, _("Failed to generate periodogram\n"));
	break;

    case PVALUE:
	if (strcmp(line, "pvalue") == 0)
	    help("pvalue", paths.cmd_helpfile, prn);	    
	else
	    err = (batch_pvalue(line, Z, datainfo, prn) == NADBL);
	break;

    case QUIT:
	if (plp) pprintf(prn, _("Script done\n"));
	else pprintf(prn, _("Please use the Close button to exit\n"));
	break;

    case RHODIFF:
	if (!command.list[0]) {
	    pprintf(prn, _("This command requires a list of variables\n"));
	    err = 1;
	    break;
	}
	err = rhodiff(command.param, command.list, &Z, datainfo);
	if (err) errmsg(err, prn);
	else varlist(datainfo, prn);
	break;

    case RUN:
	err = getopenfile(line, runfile, &paths, 1, 1);
	if (err) { 
	    pprintf(prn, _("Run command failed\n"));
	    break;
	}
	if (myname != NULL && strcmp(runfile, myname) == 0) { 
	    pprintf(prn, _("Infinite loop detected in script\n"));
	    return 1;
	}
	execute_script(runfile, NULL, NULL, NULL, prn, SESSION_EXEC);
	/* pprintf(prn, "run command not available in script mode\n"); */
	break;

    case SCATTERS:
        if (plp != NULL) /* fixme? */
            pprintf(prn, _("scatters command not available in batch mode\n"));
        else {
            err = multi_scatters(command.list, atoi(command.param), &Z, 
                                 datainfo, &paths);
            if (err) pprintf(prn, _("scatters command failed\n"));
        }               
        break;

    case SEED:
	srand((unsigned) atoi(command.param));
	pprintf(prn, _("Pseudo-random number generator seeded with %d\n"),
		atoi(command.param));
	break;

    case SETOBS:
	err = set_obs(line, datainfo, oflag);
	if (err) errmsg(err, prn);
	else {
	    set_sample_label(datainfo);
	    print_smpl(datainfo, 0, prn);
	}
	break;	

    case SETMISS:
        set_miss(command.list, command.param, Z, datainfo, prn);
        break;

    case SHELL:
	pprintf(prn, _("shell command not implemented in script mode\n"));
	break;

    case SIM:
	err = simulate(line, &Z, datainfo);
	if (err) errmsg(err, prn);
	break;

    case SMPL:
	if (oflag) {
	    restore_sample(NULL, 0, NULL);
	    if ((subinfo = malloc(sizeof *subinfo)) == NULL) 
		err = E_ALLOC;
	    else 
		err = set_sample_dummy(line, &Z, &subZ, datainfo, 
				       subinfo, oflag);
	    if (!err) {
		/* save the full data set for later use */
		fullZ = Z;
		fullinfo = datainfo;
		datainfo = subinfo;
		Z = subZ;		
	    }
	} 
	else if (strcmp(line, "smpl full") == 0) {
	    restore_sample(NULL, 0, NULL);
	    restore_sample_state(FALSE);
	    check = 1;
	} else 
	    err = set_sample(line, datainfo);
	if (err) errmsg(err, prn);
	else {
	    print_smpl(datainfo, (oflag)? fullinfo->n : 0, prn);
	    if (oflag) 
		set_sample_label_special();
	    else
		set_sample_label(datainfo);
	    if (!check) restore_sample_state(TRUE);
	}
	break;

    case SQUARE:
	if (oflag) check = xpxgenr(command.list, &Z, datainfo, 1, 1);
	else check = xpxgenr(command.list, &Z, datainfo, 0, 1);
	if (check < 0) {
	    pprintf(prn, _("Failed to generate squares\n"));
	    err = 1;
	} else {
	    pprintf(prn, _("Squares generated OK\n"));
	    varlist(datainfo, prn);
	}
	break;

    case STORE:
	if ((err = command.errcode)) {
	    errmsg(command.errcode, prn);
	    break;
	}
	if (strlen(command.param)) {
	    if (oflag == OPT_Z && !has_gz_suffix(command.param))
		pprintf(prn, _("store: using filename %s.gz\n"), command.param);
	    else
		pprintf(prn, _("store: using filename %s\n"), command.param);
	} else {
	    pprintf(prn, _("store: no filename given\n"));
	    break;
	}
	if (write_data(command.param, command.list, 
		       Z, datainfo, data_option(oflag), NULL)) {
	    pprintf(prn, _("write of data file failed\n"));
	    err = 1;
	    break;
	}
	pprintf(prn, _("Data written OK.\n"));
	if ((oflag == OPT_O || oflag == OPT_S) && datainfo->markers) 
	    pprintf(prn, _("Warning: case markers not saved in "
		    "binary datafile\n"));
	break;

    case TESTUHAT:
	if ((err = script_model_test(0, prn, 0))) break;
	if (genr_fit_resid(models[0], &Z, datainfo, GENR_RESID, 1)) {
	    pprintf(prn, _("Out of memory attempting to add variable\n"));
	    err = 1;
	    break;
	}
	freq = freqdist(&Z, datainfo, datainfo->v - 1, (models[0])->ncoeff);
	dataset_drop_vars(1, &Z, datainfo);
	if (!(err = freq_error(freq, prn))) {
	    if (rebuild) {
		normal_test(ptest, freq);
		add_test_to_model(ptest, models[0]);
	    }
	    printfreq(freq, prn); 
	    free_freq(freq);
	}
	break;

    case TSLS:
	if (rebuilding)
	    save_model_copy(&models[0], psession, rebuild, datainfo);
	else
	    clear_model(models[0], datainfo); 
	*models[0] = tsls_func(command.list, atoi(command.param), 
			       &Z, datainfo);
	if ((err = (models[0])->errcode)) {
	    errmsg((models[0])->errcode, prn);
	    break;
	}
	++model_count;
	(models[0])->ID = model_count;
	if (printmodel(models[0], datainfo, prn))
	    (models[0])->errcode = E_NAN;
	/* is this OK? */
	if (oflag) outcovmx(models[0], datainfo, 0, prn); 
	break;		

    case VAR:
	order = atoi(command.param);
	err = var(order, command.list, &Z, datainfo, 0, prn);
	break;

    case 999:
	err = 1;
	break;

    default:
	pprintf(prn, _("Sorry, the %s command is not yet implemented "
		"in libgretl\n"), command.cmd);
	break;
    } /* end of command switch */

    /* log the specific command? */
    if (exec_code == CONSOLE_EXEC) {
	if (err) 
	    strcpy(errline, line);
	else {
	    cmd_init(line);
	    errline[0] = '\0';
	}
    }

    if (is_model_cmd(command.cmd) && !err) 
	err = stack_model(0);
		
    if (err) return 1;
    else return 0;
}

/* ........................................................... */

void view_script_default (void)
     /* for "session" use */
{
    if (dump_cmd_stack(cmdfile)) return;

    view_file(cmdfile, 0, 0, 77, 350, EDIT_SCRIPT, NULL);
}

/* .................................................................. */

struct search_replace {
    GtkWidget *w;
    GtkWidget *f_entry;
    GtkWidget *r_entry;
    gchar *f_text;
    gchar *r_text;
};

/* .................................................................. */

static void replace_string_callback (GtkWidget *widget, 
				     struct search_replace *s)
{
    s->f_text = 
	gtk_editable_get_chars(GTK_EDITABLE(s->f_entry), 0, -1);
    s->r_text = 
	gtk_editable_get_chars(GTK_EDITABLE(s->r_entry), 0, -1);
    gtk_widget_destroy(s->w);
}

/* .................................................................. */

static void trash_replace (GtkWidget *widget, 
			   struct search_replace *s)
{
    s->f_text = NULL;
    s->r_text = NULL;
    gtk_widget_destroy(s->w);
}

/* .................................................................. */

static void replace_string_dialog (struct search_replace *s)
{
    GtkWidget *label, *button, *hbox;

    s->w = gtk_dialog_new();

    gtk_window_set_title (GTK_WINDOW (s->w), _("gretl: replace"));
    gtk_container_border_width (GTK_CONTAINER (s->w), 5);

    /* Find part */
    hbox = gtk_hbox_new(TRUE, TRUE);
    label = gtk_label_new(_("Find:"));
    gtk_widget_show (label);
    s->f_entry = gtk_entry_new();
    gtk_widget_show (s->f_entry);
    gtk_box_pack_start (GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_box_pack_start (GTK_BOX(hbox), s->f_entry, TRUE, TRUE, 0);
    gtk_widget_show (hbox);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (s->w)->vbox), 
                        hbox, TRUE, TRUE, 5);

    /* Replace part */
    hbox = gtk_hbox_new(TRUE, TRUE);
    label = gtk_label_new(_("Replace with:"));
    gtk_widget_show (label);
    s->r_entry = gtk_entry_new();
    gtk_signal_connect(GTK_OBJECT (s->r_entry), 
			"activate", 
		       GTK_SIGNAL_FUNC (replace_string_callback), s);
    gtk_widget_show (s->r_entry);
    gtk_box_pack_start (GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_box_pack_start (GTK_BOX(hbox), s->r_entry, TRUE, TRUE, 0);
    gtk_widget_show (hbox);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (s->w)->vbox), 
                        hbox, TRUE, TRUE, 5);

    gtk_box_set_spacing(GTK_BOX (GTK_DIALOG (s->w)->action_area), 15);
    gtk_box_set_homogeneous(GTK_BOX 
			     (GTK_DIALOG (s->w)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW (s->w), GTK_WIN_POS_MOUSE);

    gtk_signal_connect(GTK_OBJECT(s->w), "destroy",
		       gtk_main_quit, NULL);

    /* replace button -- make this the default */
    button = gtk_button_new_with_label (_("Replace all"));
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (s->w)->action_area), 
		       button, TRUE, TRUE, FALSE);
    gtk_signal_connect(GTK_OBJECT (button), "clicked",
		       GTK_SIGNAL_FUNC (replace_string_callback), s);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* cancel button */
    button = gtk_button_new_with_label (_("Cancel"));
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (s->w)->action_area), 
		       button, TRUE, TRUE, FALSE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(trash_replace), s);
    gtk_widget_show(button);

    gtk_widget_grab_focus(s->f_entry);
    gtk_widget_show (s->w);
    gtk_main();
}

/* ........................................................... */

void text_replace (windata_t *mydata, guint u, GtkWidget *widget)
{
    gchar *buf;
    int count = 0;
    gint pos = 0;
    guint sel_start, sel_end;
    size_t sz, fullsz, len, diff;
    char *replace = NULL, *find = NULL;
    char *modbuf, *p, *q;
    gchar *old;
    struct search_replace *s;
    GtkEditable *gedit = GTK_EDITABLE(mydata->w);

    s = mymalloc(sizeof *s);
    if (s == NULL) return;

    replace_string_dialog(s);

    if (s->f_text == NULL || s->r_text == NULL) {
	free(s);
	return;
    }

    find = s->f_text;
    replace = s->r_text;

    if (!strlen(find)) {
	free(find);
	free(replace);
	free(s);
	return;
    }

    if (gedit->has_selection) {
	sel_start = gedit->selection_start_pos;
	sel_end = gedit->selection_end_pos;
    } else {
	sel_start = 0;
	sel_end = 0;
    }

    buf = gtk_editable_get_chars(gedit, sel_start, (sel_end)? sel_end : -1);
    if (buf == NULL || !(sz = strlen(buf))) 
	return;

    fullsz = gtk_text_get_length(GTK_TEXT(mydata->w));
    len = strlen(find);
    diff = strlen(replace) - len;

    p = buf;
    while (*p) {
	if ((q = strstr(p, find))) {
	    count++;
	    p = q + 1;
	}
	else break;
    }
    if (count) {
	fullsz += count * diff;
    } else {
	errbox(_("String to replace was not found"));
	free(buf);
	return;
    }

    modbuf = mymalloc(fullsz + 1);
    if (modbuf == NULL) {
	free(find);
	free(replace);
	free(s);
	return;
    }

    *modbuf = '\0';

    if (sel_start) {
	gchar *tmp = gtk_editable_get_chars(gedit, 0, sel_start);

	strcat(modbuf, tmp);
	g_free(tmp);
    }

    p = buf;
    while (*p) {
	if ((q = strstr(p, find))) {
	    strncat(modbuf, p, q - p);
	    strcat(modbuf, replace);
	    p = q + len;
	} else {
	    strcat(modbuf, p);
	    break;
	}
    }

    if (sel_end) {
	gchar *tmp = gtk_editable_get_chars(gedit, sel_end, -1);

	strcat(modbuf, tmp);
	g_free(tmp);
    }    

    /* save original buffer for "undo" */
    old = gtk_object_get_data(GTK_OBJECT(mydata->w), "undo");
    if (old != NULL) {
	g_free(old);
	gtk_object_remove_data(GTK_OBJECT(mydata->w), "undo");
    }
    gtk_object_set_data(GTK_OBJECT(mydata->w), "undo", 
			gtk_editable_get_chars(gedit, 0, -1));

    /* now insert the modified buffer */
    gtk_text_freeze(GTK_TEXT(mydata->w));
    gtk_editable_delete_text(gedit, 0, -1);
    gtk_editable_insert_text(gedit, modbuf, strlen(modbuf), &pos);
    gtk_text_thaw(GTK_TEXT(mydata->w));

    /* and clean up */
    free(find);
    free(replace);
    free(s);
    free(modbuf);
    g_free(buf);
}

#include "../cli/common.c"
