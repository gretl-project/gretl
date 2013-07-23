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

/* gretl audio graph plugin */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "libgretl.h"
#include "version.h"

#include "miditypes.h"
#include "midi_utils.h"

#include <glib.h>
#include <gtk/gtk.h>

#define ADEBUG 0

#ifdef HAVE_FLITE
# include <flite/flite.h>
extern cst_voice *register_cmu_us_kal (void);
#endif

#ifdef WIN32_SAPI
# define COBJMACROS
# include <sapi.h>
#endif

#if defined(HAVE_FLITE) || defined(WIN32_SAPI)
# define DO_SPEECH
# include "libgretl.h"
#endif

const char *track_hdr = "MTrk";

enum dataset_comments {
    TS_COMMENT = 0,
    YLABEL_COMMENT,
    XLABEL_COMMENT,
    XRANGE_COMMENT,
    N_COMMENTS
};

#define DEFAULT_FORCE 96 /* volume of MIDI notes */

#ifndef G_OS_WIN32
# define DEFAULT_MIDI_PROG "timidity" /* MIDI player to be used */
# define DEFAULT_MIDI_OPT  "-A500"    /* option to pass to MIDI player */
#endif
    
typedef struct _datapoint datapoint;
typedef struct _dataset dataset;
typedef struct _midi_spec midi_spec;
typedef struct _midi_track midi_track;
typedef struct _note_event note_event;

struct _datapoint {
    double x;
    double y;
};

struct _dataset {
    int pd;
    int n;
    int series2;
    double intercept;
    double slope;
    datapoint *points;
    double *y2;
    gchar *comments[N_COMMENTS];
};

struct _midi_spec {
    int ntracks;
    int nticks;
    int nsecs;
    dataset *dset;
    FILE *fp;
};

struct _midi_track {
    unsigned char channel;
    unsigned char patch;
    int n_notes;
    note_event *notes; 
};

struct _note_event {
    double dtime;
    double duration;
    unsigned char pitch;
    unsigned char force;
    unsigned char channel;
};

#define delta_time_zero(f) (putc(0x00, f))

#define YMAX 60

static int write_note_event (note_event *event, midi_spec *spec)
{
    int len;

    len = delta_time(event->dtime, spec->nticks, spec->fp);

    putc(MIDI_NOTE_ON + event->channel, spec->fp);
    putc(event->pitch, spec->fp);
    putc(event->force, spec->fp);
    len += 3;

    /* use "running status" */
    len += delta_time(event->duration, spec->nticks, spec->fp);
    putc(event->pitch, spec->fp);    
    putc(0, spec->fp);
    len += 2;

    return len;
}

static void write_midi_track (midi_track *track, 
			      midi_spec *spec)
{
    unsigned char bal = 0x7f * (track->channel % 2);
    char tmp[32];
    long len = 0;
    int i, n;

    fwrite(track_hdr, 1, 4, spec->fp);
    write_be_long(0, spec->fp); /* revisited below */

    delta_time_zero(spec->fp);
    len++;
    sprintf(tmp, "Series %d", track->channel + 1);
    n = strlen(tmp);
    len += write_midi_meta(MIDI_TRACK_NAME, spec->fp);
    len += write_var_len(n, spec->fp);
    for (i=0; i<n; i++) {
	putc(tmp[i], spec->fp);
    }
    len += n;

    delta_time_zero(spec->fp);
    len++;
    strcpy(tmp, get_patch_name(track->patch));
    n = strlen(tmp);
    len += write_midi_meta(MIDI_INSTR_NAME, spec->fp);
    len += write_var_len(n, spec->fp);
    for (i=0; i<n; i++) {
	putc(tmp[i], spec->fp);
    }
    len += n;

    /* Stereo separation of series */
    delta_time_zero(spec->fp);
    fputc(0xb0 + track->channel, spec->fp); /* Controller on channel */ 
    fputc(MIDI_PAN_MSB, spec->fp);
    fputc(bal, spec->fp);  /* MSB = 0 or 0x7f */
    delta_time_zero(spec->fp);
    fputc(0xb0 + track->channel, spec->fp);
    fputc(MIDI_PAN_LSB, spec->fp);
    fputc(bal, spec->fp);  /* LSB = 0 or 0x7f */
    len += 8;

    delta_time_zero(spec->fp);
    len++;
    len += write_midi_byte(MIDI_PROGRAM_CHANGE + track->channel, spec->fp);
    len += write_midi_byte(track->patch, spec->fp);

    for (i=0; i<track->n_notes; i++) {
	track->notes[i].channel = track->channel;
	len += write_note_event(&track->notes[i], spec);
    }

    delta_time_zero(spec->fp);
    len++;
    len += write_midi_eot(spec->fp);

    fseek(spec->fp, -(len + 4), SEEK_CUR);
    /* now fix the length */
    write_be_long(len, spec->fp); 
    fseek(spec->fp, len, SEEK_CUR);
}

static void four_four_header (midi_spec *spec)
{
    long min_note_time = 40000;
    long max_note_time = 180000;
    int ticklen, len = 0;
    long pos, pos1;

    write_midi_header(1, spec->ntracks + 1, spec->nticks, spec->fp);

    /* tempo/time track */
    fwrite(track_hdr, 1, 4, spec->fp);
    pos = ftell(spec->fp);
    write_be_long(0, spec->fp); /* revisited below */

#if 0
    /* time signature (4/4 is the MIDI default) */
    delta_time_zero(spec->fp);
    len++;
    len += write_midi_meta(MIDI_TIME_SIG, spec->fp);
    /* size (bytes) */
    len += write_midi_byte(4, spec->fp);
    /* 4/4 time */
    len += write_midi_byte(4, spec->fp);
    len += write_midi_byte(2, spec->fp);
    /* 24 MIDI clocks per metronome click */
    len += write_midi_byte(24, spec->fp);
    /* 8 32nd notes per tick */
    len += write_midi_byte(8, spec->fp);
#endif

    /* tempo */
    delta_time_zero(spec->fp);
    len++;
    len += write_midi_meta(MIDI_TEMPO, spec->fp);
    /* bytes */
    len += write_midi_byte(3, spec->fp);
    if (spec->dset == NULL) {
	/* microseconds per quarter note */
	len += write_be_24(500000, spec->fp);
    } else {
	/* total number of microseconds */
	double nms = spec->nsecs * 1.0e6;
	/* microseconds per quarter note */
	long msq = nms / spec->dset->n;

#if ADEBUG
	fprintf(stderr, "Secs = %d, microsecs = %g, msq = %d\n",
		spec->nsecs, nms, (int) msq);
#endif 

	if (msq > max_note_time) {
	    msq = max_note_time;
	} else if (msq < min_note_time) {
	    msq = min_note_time;
	}

#if ADEBUG
	fprintf(stderr, "After adj, msq = %d\n", (int) msq);
#endif 

	len += write_be_24(msq, spec->fp);
    }

    /* end */
    if (spec->dset == NULL) {
	ticklen = 4 * spec->nticks;
#if ADEBUG
	fprintf(stderr, "dset is NULL, ticklen = %d\n", ticklen);
#endif 
    } else {
	ticklen = spec->dset->n * spec->nticks;
#if ADEBUG
	fprintf(stderr, "ticklen = %d * %d = %d\n", 
		spec->dset->n, spec->nticks, ticklen);
#endif 
    }
    len += write_var_len(ticklen, spec->fp);
    len += write_midi_eot(spec->fp);

    /* length of header */
    pos1 = ftell(spec->fp);
    fseek(spec->fp, pos, SEEK_SET);
    write_be_long(len, spec->fp);
    fseek(spec->fp, pos1, SEEK_SET);
}

static void dataset_init (dataset *dset)
{
    int i;

    dset->points = NULL;
    dset->y2 = NULL;

    dset->n = 0;
    dset->pd = 0;
    dset->series2 = 0;

    dset->intercept = NADBL;
    dset->slope = NADBL;

    for (i=0; i<N_COMMENTS; i++) {
	dset->comments[i] = NULL;
    }
}

static void dataset_free (dataset *dset)
{
    int i;

    free(dset->points);

    if (dset->y2 != NULL) {
	free(dset->y2);
    }
    
    for (i=0; i<N_COMMENTS; i++) {
	if (dset->comments[i] != NULL) {
	    g_free(dset->comments[i]);
	}
    }
}

const char *cent_str (int cent)
{
    switch (cent) {
    case 14:
	return "fourteen";
    case 15:
	return "fifteen";
    case 16:
	return "sixteen";
    case 17:
	return "seventeen";
    case 18:
	return "eighteen";
    case 19:
	return "nineteen";
    case 20:
	return "two thousand";
    default:
	return "unknown";
    }
}

static int make_xrange_comment (const char *line, dataset *dset)
{
    double x1, x2;

    if (sscanf(line, "set xrange [%lf:%lf]", &x1, &x2) != 2) {
	return 0;
    }

    if (dset->pd == 0) {
	dset->comments[XRANGE_COMMENT] = 
	    g_strdup_printf("x range %.4g to %.4g", x1, x2);
    } else {
	char tmp1[64], tmp2[64];
	int ix1 = x1;
	int ix2 = x2;
	int cent, yr;

	cent = x1 / 100;
	yr = ix1 - 100 * cent;
	sprintf(tmp1, "%s %d", cent_str(cent), yr);
    
	cent = x2 / 100;
	yr = ix2 - 100 * cent;
	sprintf(tmp2, "%s %d", cent_str(cent), yr);

	dset->comments[XRANGE_COMMENT] = 
	    g_strdup_printf("%s to %s", tmp1, tmp2);
    }

    return 1;
}

static int make_axis_label_comment (const char *line, dataset *dset)
{
    char tmp[16];
    int ret = 0;

    if (sscanf(line, "set ylabel '%15[^']", tmp) == 1) {
	if (dset->pd > 0) {
	    dset->comments[YLABEL_COMMENT] = g_strdup(tmp);
	} else {
	    dset->comments[YLABEL_COMMENT] =
		g_strdup_printf("y variable %s", tmp);
	}
	ret = 1;
    } else if (sscanf(line, "set xlabel '%15[^']", tmp) == 1) {
	if (dset->pd > 0) {
	    dset->comments[XLABEL_COMMENT] = g_strdup(tmp);
	} else {
	    dset->comments[XLABEL_COMMENT] =
		g_strdup_printf("x variable %s", tmp);
	}
	ret = 1;
    }

    return ret;
}

static void make_ts_comment (const char *line, dataset *dset)
{
    const char *pdstr = NULL;

    sscanf(line, "# timeseries %d", &dset->pd);

    switch (dset->pd) {
    case 1:
	pdstr = "Annual";
	break;
    case 4:
	pdstr = "Quarterly";
	break;
    case 5:
    case 7:
	pdstr = "Daily";
	break;
    case 12:
	pdstr = "Monthly";
	break;
    case 24:
	pdstr = "Hourly";
	break;
    default:
	break;
    }

    if (pdstr != NULL) {
	dset->comments[TS_COMMENT] = 
	    g_strdup_printf("%s time series", pdstr);
    } else {
	dset->comments[TS_COMMENT] = g_strdup("time series");
    }
}

static int get_comment (const char *line, dataset *dset)
{
    int ret = 0;

    if (!strncmp(line, "# timeseries", 12)) {
	make_ts_comment(line, dset);
	ret = 1;
    } else if (!strncmp(line, "set ylabel", 10) ||
	     !strncmp(line, "set xlabel", 10)) {
	make_axis_label_comment(line, dset);
	ret = 1;
    } else if (!strncmp(line, "set xrange", 10)) {
	make_xrange_comment(line, dset);
	ret = 1;
    } 

    return ret;
}

#ifdef PLAY_AUTOFIT_LINE

static int get_fit_params (char *line, dataset *dset)
{
    double a, b;
    char *p = strchr(line, '*');

    if (p != NULL) *p = '\0';

    if (sscanf(line, "%lf + %lf", &a, &b) == 2) {
	dset->intercept = a;
	dset->slope = b;
	dset->series2 = 1;
	return 1;
    }

    return 0;
}

#endif

static int get_data_x_y (const char *line, double *x, double *y)
{
    if (sscanf(line, "%lf %lf", x, y) == 2) return 0;

    if (sscanf(line, "%lf ?", x) == 1) {
	*y = NADBL;
	return 0;
    }

    return 1;
}

static int read_datafile (const char *fname, dataset *dset)
{
#ifdef PLAY_AUTOFIT_LINE
    int fitline = 0;
#endif
    char line[256];
    int i, err = 0;
    int got_e = 0, y2data = 0;
    FILE *fdat;

    dataset_init(dset);

    fdat = gretl_fopen(fname, "r");
    if (fdat == NULL) {
	fprintf(stderr, "Couldn't open '%s'\n", fname);
	return 1;
    } else {
	fprintf(stderr, "Reading %s...\n", fname);
    }

    while (fgets(line, sizeof line, fdat)) {
	tailstrip(line);
	if (get_comment(line, dset)) {
	    continue;
	} else if (!strcmp(line, "e")) {
	    fprintf(stderr, "Got end of data marker\n");
	    got_e++;
	    if (got_e == 2) {
		/* can't handle more than two series! */
		break;
	    }
	} else if (strstr(line, "automatic fitted")) {
#ifdef PLAY_AUTOFIT_LINE
	    fitline = 1;
#endif
	    continue;
	} else if (isdigit((unsigned char) line[0])) {
	    if (strstr(line, "title")) {
#ifdef PLAY_AUTOFIT_LINE
		if (fitline) {
		    get_fit_params(line, dset);
		}
#endif
		continue;
	    }
	    if (!got_e) {
		dset->n += 1;
	    } else if (!y2data) {
		y2data = 1;
	    }
	}
    }

    if (dset->n == 0) {
	fprintf(stderr, "No data in '%s'\n", fname);
	err = 1;
	goto bailout;
    } 

    dset->points = malloc(dset->n * sizeof *dset->points);
    if (dset->points == NULL) {
	err = 1;
	fputs("Out of memory\n", stderr);
	goto bailout;
    }

    if (y2data) {
	dset->y2 = malloc(dset->n * sizeof *dset->y2);
	if (dset->y2 == NULL) {
	    err = 1;
	    fputs("Out of memory\n", stderr);
	    goto bailout;
	}
	dset->series2 = 1;
    }

    rewind(fdat);

    i = got_e = 0;
    while (!err && fgets(line, 256, fdat)) {
	tailstrip(line);
	if (!strcmp(line, "e")) {
	    got_e++;
	    if (got_e == 2) {
		break;
	    } 
	    i = 0;
	} else if (isdigit((unsigned char) line[0])) {
	    double x, y;

	    if (strstr(line, "title")) {
		continue;
	    }

	    if (get_data_x_y(line, &x, &y)) {
		fprintf(stderr, "Couldn't read data on line %d\n", i + 1);
		err = 1;
	    } else {
		if (!got_e) {
		    dset->points[i].x = x;
		    dset->points[i].y = y;
		} else {
		    dset->y2[i] = y;
		}
	    }
	    i++;
	}
    }    

 bailout:

    fclose(fdat);

    if (err) {
	dataset_free(dset);
    }

    return err;
}

static void points_min_max (const datapoint *points, 
			    double *xmin, double *xmax,
			    double *ymin, double *ymax,
			    int n)
{
    int i;

    *xmin = *xmax = points[0].x;
    *ymin = *ymax = points[0].y;

    for (i=1; i<n; i++) {
	if (points[i].x < *xmin) *xmin = points[i].x;
	if (points[i].x > *xmax) *xmax = points[i].x;
	if (!na(points[i].y)) {
	    if (points[i].y < *ymin) *ymin = points[i].y;
	    if (points[i].y > *ymax) *ymax = points[i].y;
	}
    }
}

static void min_max (const double *x, double *min, double *max,
		     int n)
{
    int i;

    *min = *max = x[0];

    for (i=1; i<n; i++) {
	if (!na(x[i])) {
	    if (x[i] < *min) *min = x[i];
	    if (x[i] > *max) *max = x[i];
	}
    }
}

#if defined(HAVE_FLITE)

#if 0
static int save_dataset_comments (const dataset *dset)
{
    int i, j;
    cst_voice *v;
    cst_wave *w, *fullw = NULL;

    flite_init();
    v = register_cmu_us_kal();

    j = 0;
    for (i=0; i<N_COMMENTS; i++) {
	if (dset->comments[i] != NULL) {
	    if (j == 0) {
		fullw = flite_text_to_wave(dset->comments[i], v);
	    } else {
		w = flite_text_to_wave(dset->comments[i], v);
		concat_wave(fullw, w);
		delete_wave(w);
	    }
	    j++;
	}
    }

    cst_wave_save_riff(fullw, "gretl_flite.wav");
    delete_wave(fullw);

    return 0;
}
#endif

static void speak_dataset_comments (const dataset *dset)
{
    int i;
    cst_voice *v;

    flite_init();
    v = register_cmu_us_kal();

    for (i=0; i<N_COMMENTS; i++) {
	if (dset->comments[i] != NULL) {
	    flite_text_to_speech(dset->comments[i], v, "play");
	}
    }
}

#elif defined(WIN32_SAPI)

static WCHAR *wide_string (const char *s)
{
    unsigned char *ret;
    int i, j, n = strlen(s);

    ret = malloc(2 * (n + 1));
    if (ret == NULL) return NULL;

    j = 0;
    for (i=0; i<=n; i++) {
	ret[j++] = s[i];
	ret[j++] = '\0';
    }

    return (WCHAR *) ret;
}

static void speak_dataset_comments (const dataset *dset)
{
    int i;
    ISpVoice *v = NULL;
    HRESULT hr;

    hr = CoInitialize(NULL);
    if (!SUCCEEDED(hr)) {
        fprintf(stderr, "CoInitialize failed\n");
        return;
    }

    hr = CoCreateInstance(&CLSID_SpVoice, 
                          NULL, 
                          CLSCTX_ALL, 
                          &IID_ISpVoice, 
                          (void *) &v); 

    if (SUCCEEDED(hr)) {
	for (i=0; i<N_COMMENTS; i++) {
	    if (dset->comments[i] != NULL) {
		WCHAR *w = wide_string(dset->comments[i]);

		ISpVoice_Speak(v, w, 0, NULL);
		free(w);
	    }
	}
        ISpVoice_Release(v);
    } 

    CoUninitialize();
}

#endif

static void print_dataset_comments (const dataset *dset)
{
    int i;

    for (i=0; i<N_COMMENTS; i++) {
	if (dset->comments[i] != NULL) {
	    printf("Comment[%d]: %s\n", i+1, dset->comments[i]);
	}
    }

#ifdef DO_SPEECH
    speak_dataset_comments(dset);
#endif
}

static void audio_graph_error (const char *msg)
{
#ifdef HAVE_FLITE
    cst_voice *v;
#endif
#ifdef WIN32_SAPI
    ISpVoice *v = NULL;
    HRESULT hr;
#endif    

    fprintf(stderr, "%s\n", msg);

#ifdef HAVE_FLITE
    flite_init();

    v = register_cmu_us_kal();
    flite_text_to_speech(msg, v, "play");
#endif
#ifdef WIN32_SAPI
    hr = CoInitialize(NULL);
    if (!SUCCEEDED(hr)) return;

    hr = CoCreateInstance(&CLSID_SpVoice, 
                          NULL, 
                          CLSCTX_ALL, 
                          &IID_ISpVoice, 
                          (void *) &v);
    if (SUCCEEDED(hr)) {
	wchar_t *w = wide_string(msg);

	ISpVoice_Speak(v, w, 0, NULL);
	free(w);
        ISpVoice_Release(v);
    } 
#endif
}

static int compare_points (const void *a, const void *b)
{
    const datapoint *pa = (const datapoint *) a;
    const datapoint *pb = (const datapoint *) b;
     
    return (pa->x > pb->x) - (pa->x < pb->x);
}

static int play_dataset (midi_spec *spec, midi_track *track,
			 const dataset *dset)
{
    double xmin, xmax, xavg;
    double ymin, ymax;
    double xscale, yscale;
    int i;

    track->notes = malloc(dset->n * sizeof *track->notes); 
    if (track->notes == NULL) {
	fputs("out of memory\n", stderr);
	return 1;
    }

    track->channel = 0;
    track->patch = PC_GRAND; 
    track->n_notes = dset->n;

    if (dset->pd == 0) {
	/* scatter plot: sort data by x value */
	qsort((void *) dset->points, (size_t) dset->n, 
	      sizeof dset->points[0], compare_points);
    }

    points_min_max(dset->points, &xmin, &xmax, &ymin, &ymax, dset->n);
   
    if (dset->y2 != NULL) {
	double y2min, y2max;

	min_max(dset->y2, &y2min, &y2max, dset->n);
	if (y2min < ymin) ymin = y2min;
	if (y2max > ymax) ymax = y2max;
    }

    xavg = (xmax - xmin) / (dset->n - 1);
    /* normalize average x step to quarter note */
    xscale = 1.0 / xavg;

#if ADEBUG
    fprintf(stderr, "xavg = %g, xscale = %g\n", xavg, xscale);
#endif

    yscale = YMAX / (ymax - ymin);

    for (i=0; i<dset->n; i++) {
	double dtx, dux, ypos;

	if (!na(dset->points[i].y)) {	
	    ypos = (dset->points[i].y - ymin) * yscale;
	} else {
	    ypos = NADBL;
	}

	if (i == 0) {
	    dtx = 0.0;
	} else {
	    dtx = xscale * (dset->points[i].x - dset->points[i-1].x);
	}

	if (i == dset->n - 1) {
	    dux = xscale * xavg;
	} else {
	    dux = xscale * (dset->points[i+1].x - dset->points[i].x);
	}

	track->notes[i].dtime = dtx;
	track->notes[i].duration = dux;
	if (!na(ypos)) {
	    track->notes[i].pitch = 36 + (int) (ypos + 0.5);
	    track->notes[i].force = DEFAULT_FORCE;
	} else {
	    track->notes[i].pitch = 36;
	    track->notes[i].force = 0;
	}	    

#if ADEBUG
	fprintf(stderr, "Obs %d: x = %g, y = %g, ypos = %g\n", i, 
		dset->points[i].x, dset->points[i].y, ypos);
	fprintf(stderr, " dtime=%g, duration=%g, pitch=%d\n",
		track->notes[i].dtime, track->notes[i].duration,
		track->notes[i].pitch);
#endif
    }

    write_midi_track(track, spec);

    if (dset->series2) {
	for (i=0; i<dset->n; i++) {
	    double yi, ypos = 0;

	    if (dset->y2 != NULL) {
		yi = dset->y2[i];
	    } else {
		yi = dset->intercept + dset->slope * dset->points[i].x;
	    }

	    if (!na(yi)) {
		ypos = (yi - ymin) * yscale;
		track->notes[i].pitch = 36 + (int) (ypos + 0.5);
		track->notes[i].force = DEFAULT_FORCE;
	    } else {
		track->notes[i].pitch = 36;
		track->notes[i].force = 0;
	    }
#if ADEBUG
	    fprintf(stderr, "Series2, Obs %d: x = %g, y = %g, ypos = %g\n", i, 
		    dset->points[i].x, yi, ypos);
	    fprintf(stderr, " dtime=%g, duration=%g, pitch=%d\n",
		    track->notes[i].dtime, track->notes[i].duration,
		    track->notes[i].pitch);
#endif
	}

	track->channel = 1;
	/* below: was PC_MARIMBA, but that's not in freepats */
	track->patch = PC_CELLO;
	write_midi_track(track, spec);
    }	

    free(track->notes);

    return 0;
}

#if defined(G_OS_WIN32)

static int midi_fork (const char *fname, const char *midiplayer)
{
    int err = 0;

#ifdef _WIN64
    if ((gint64) ShellExecute(NULL, "open", fname, NULL, NULL, SW_SHOW) < 32) {
	err = 1;
    }
#else
    if ((int) ShellExecute(NULL, "open", fname, NULL, NULL, SW_SHOW) < 32) {
	err = 1;
    }
#endif

    return err;
}

#else /* non-Windows version */

static char **get_midi_args (int *argc, const char *fname,
			     const char *midiplayer)
{
    char **argv = NULL;
    const char *midistr;

    if (midiplayer != NULL && *midiplayer != '\0') {
	midistr = midiplayer;
    } else {
	midistr = (const char *) getenv("GRETL_MIDI_PLAYER");
    }
	
    if (midistr == NULL || *midistr == '\0') {
	*argc = 3;

	argv = malloc((*argc + 1) * sizeof *argv);
	if (argv == NULL) return NULL;

	argv[0] = g_strdup(DEFAULT_MIDI_PROG);
	argv[1] = g_strdup(DEFAULT_MIDI_OPT);
	argv[2] = g_strdup(fname);
	argv[3] = NULL;
    } else {
	char *tmp = g_strdup(midistr);
	char *s = tmp;
	int i, ntoks = 1;

	while (*s) {
	    if (*s == ' ') ntoks++;
	    s++;
	}

	*argc = ntoks + 1;

	argv = malloc((*argc + 1) * sizeof *argv); 
	if (argv == NULL) {
	    g_free(tmp);
	    return NULL;
	}

	if (ntoks > 1) {
	    argv[0] = g_strdup(strtok(tmp, " "));
	    for (i=1; i<ntoks; i++) {
		argv[i] = g_strdup(strtok(NULL, " "));
	    }
	} else {
	    argv[0] = g_strdup(tmp);
	    i = 1;
	}
	
	argv[i++] = g_strdup(fname);
	argv[i] = NULL;
	g_free(tmp);
    }

    return argv;
}

static int real_midi_fork (char **argv)
{
    gboolean run;

    run = g_spawn_async(NULL, argv, NULL, G_SPAWN_SEARCH_PATH, 
                        NULL, NULL, NULL, NULL);

    return !run;
}

static int midi_fork (const char *fname, const char *midiplayer)
{
    int i, err, argc;
    char **argv;

    argv = get_midi_args(&argc, fname, midiplayer);
    if (argv == NULL) return 1;

    err = real_midi_fork(argv);

    for (i=0; i<argc; i++) {
	g_free(argv[i]);
    }
    free(argv);  

    return err;
}

#endif /* ! MS Windows */

int midi_play_graph (const char *fname, const char *userdir,
		     const char *midiplayer)
{
    char outname[FILENAME_MAX];
    midi_spec spec;
    midi_track track;
    dataset dset;

    sprintf(outname, "%sgretl.mid", userdir);

    spec.fp = gretl_fopen(outname, "wb");
    if (spec.fp == NULL) {
	fprintf(stderr, "Couldn't write to '%s'\n", outname);
	return 1;
    }

    if (read_datafile(fname, &dset)) {
	audio_graph_error("Error reading data file");
	fclose(spec.fp);
	return 1;
    }

    spec.ntracks = 1 + dset.series2;
    spec.nticks = 96;
    spec.dset = &dset;
    spec.nsecs = 8;
    four_four_header(&spec);
    print_dataset_comments(&dset);
    play_dataset(&spec, &track, &dset);
    dataset_free(&dset);

    fclose(spec.fp);

    midi_fork(outname, midiplayer);

    return 0;
}

#ifdef DO_SPEECH
# include "audioprint.c"
#endif 
