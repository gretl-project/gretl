/* gretl audio graph plugin */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "libgretl.h"

#include "miditypes.h"
#include "midi_utils.h"

#include <glib.h>

#if GLIB_CHECK_VERSION(2,0,0)
# define GLIB2
# include <signal.h>
#endif /* GLIB_CHECK_VERSION */

#undef DEBUG

#define HAVE_FLITE 1

#ifdef HAVE_FLITE
# include <flite/flite.h>
extern cst_voice *register_cmu_us_kal (void);
#endif

const char *track_hdr = "MTrk";

#define N_COMMENTS 3

typedef struct _dataset dataset;
typedef struct _midi_spec midi_spec;
typedef struct _midi_track midi_track;
typedef struct _note_event note_event;

struct _dataset {
    int pd;
    int n;
    double *x;
    double *y;
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
    int dtime;
    int duration;
    unsigned char pitch;
    unsigned char force;
    unsigned char channel;
};

#define delta_time_zero(f) (putc(0x00, f))

#define YMAX   60

static int write_note_event (note_event *event, midi_spec *spec)
{
    int len;

    len = delta_time(event->dtime, spec->nticks, spec->fp);

    putc(MIDI_NOTE_ON + event->channel, spec->fp);
    putc(event->pitch, spec->fp);
    putc(event->force, spec->fp);
    len += 3;

    /* use running status */
    len += delta_time(event->duration, spec->nticks, spec->fp);
    putc(event->pitch, spec->fp);    
    putc(0, spec->fp);
    len += 2;

    return len;
}

static void write_midi_track (midi_track *track, 
			      midi_spec *spec)
{
    char tmp[32];
    long len = 0;
    unsigned char bal = 0x7f * (track->channel % 2);
    int i, n;

    fwrite(track_hdr, 1, 4, spec->fp);
    write_be_long(0, spec->fp); /* revisit below */

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
    write_be_long(len, spec->fp); /* now fix the length */
    fseek(spec->fp, len, SEEK_CUR);
}

static void four_four_header (midi_spec *spec)
{
    int len, ticklen;
    long pos, pos1;
    long min_note_time = 60000;
    long max_note_time = 180000;

    write_midi_header(1, spec->ntracks + 1, spec->nticks, spec->fp);

    /* tempo/time track */
    fwrite(track_hdr, 1, 4, spec->fp);
    pos = ftell(spec->fp);
    write_be_long(0, spec->fp); /* revisit below */

    /* time sig */
    len = write_var_len(0, spec->fp);
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

    /* tempo */
    len += write_var_len(0, spec->fp);
    len += write_midi_meta(MIDI_TEMPO, spec->fp);
    /* bytes */
    len += write_midi_byte(3, spec->fp);
    /* microseconds per quarter note */
    if (spec->dset == NULL) {
	len += write_be_24(500000, spec->fp);
    } else {
	double nms = spec->nsecs * 1.0e6;
	long msq = nms / spec->dset->n;

	if (msq > max_note_time) {
	    msq = max_note_time;
	} else if (msq < min_note_time) {
	    msq = min_note_time;
	}

	len += write_be_24(msq, spec->fp);
    }

    /* end */
    if (spec->dset == NULL) {
	ticklen = 4 * spec->nticks;
    } else {
	ticklen = spec->dset->n * spec->nticks;
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

    dset->x = NULL;
    dset->y = NULL;
    dset->y2 = NULL;

    dset->n = 0;
    dset->pd = 0;

    for (i=0; i<N_COMMENTS; i++) {
	dset->comments[i] = NULL;
    }
}

static void dataset_free (dataset *dset)
{
    int i;

    free(dset->x);
    free(dset->y);

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
	dset->comments[2] = g_strdup_printf("%.4g to %.4g", x1, x2);
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

	dset->comments[2] = g_strdup_printf("%s to %s", tmp1, tmp2);
    }

    return 1;
}

static int make_ylabel_comment (const char *line, dataset *dset)
{
    char tmp[16];
    int ret = 0;

    if (sscanf(line, "set ylabel '%15[^\']", tmp) == 1) {
	dset->comments[1] = g_strdup(tmp);
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
    case 12:
	pdstr = "Monthly";
	break;
    default:
	break;
    }

    if (pdstr != NULL) {
	dset->comments[0] = g_strdup_printf("%s time series", pdstr);
    } else {
	dset->comments[0] = g_strdup("time series");
    }
}

static int get_comment (const char *line, dataset *dset)
{
    int ret = 0;

    if (!strncmp(line, "# timeseries", 12)) {
	make_ts_comment(line, dset);
	ret = 1;
    }

    else if (!strncmp(line, "set ylabel", 10)) {
	make_ylabel_comment(line, dset);
	ret = 1;
    }
    
    else if (!strncmp(line, "set xrange", 10)) {
	make_xrange_comment(line, dset);
	ret = 1;
    }

    return ret;
}

static void tail_strip_line (char *line)
{
    int i, n = strlen(line);

    for (i=n-1; i>0; i--) {
	if (line[i] == '\n' || line[i] == ' ') {
	    line[i] = '\0';
	}
	else break;
    }
}

static int read_datafile (const char *fname, dataset *dset)
{
    char line[256];
    int i, err = 0;
    int got_e = 0, y2data = 0;
    FILE *fdat;

    dataset_init(dset);

    fdat = fopen(fname, "r");
    if (fdat == NULL) {
	fprintf(stderr, "Couldn't open '%s'\n", fname);
	return 1;
    } else {
	fprintf(stderr, "Reading %s...\n", fname);
    }

    while (fgets(line, 256, fdat)) {
	tail_strip_line(line);
	if (get_comment(line, dset)) continue;
	else if (!strcmp(line, "e")) {
	    fprintf(stderr, "Got end of data marker\n");
	    got_e++;
	    if (got_e == 2) {
		/* can't handle more than two series! */
		break;
	    }
	} else if (isdigit((unsigned char) line[0])) {
	    if (strstr(line, "title")) continue;
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

    dset->x = malloc(dset->n * sizeof *dset->x);
    dset->y = malloc(dset->n * sizeof *dset->y);
    if (dset->x == NULL || dset->y == NULL) {
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
    }

    rewind(fdat);

    i = got_e = 0;
    while (!err && fgets(line, 256, fdat)) {
	tail_strip_line(line);
	if (!strcmp(line, "e")) {
	    got_e++;
	    if (got_e == 2) {
		break;
	    } 
	    i = 0;
	} else if (isdigit((unsigned char) line[0])) {
	    double x, y;

	    if (strstr(line, "title")) continue;

	    if (sscanf(line, "%lf %lf", &x, &y) != 2) {
		fprintf(stderr, "Couldn't read data on line %d\n", i + 1);
		err = 1;
	    } else {
		if (!got_e) {
		    dset->x[i] = x;
		    dset->y[i] = y;
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

static void min_max (const double *x, double *min, double *max,
		     int n)
{
    int i;

    *min = *max = x[0];

    for (i=1; i<n; i++) {
	if (x[i] < *min) *min = x[i];
	if (x[i] > *max) *max = x[i];
    }
}

#ifdef HAVE_FLITE
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

# ifdef DEBUG
    fprintf(stderr, "dset->n_comments=%d, done\n", dset->n_comments); 
# endif
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

#ifdef HAVE_FLITE
    speak_dataset_comments(dset);
#endif
}

static void audio_graph_error (const char *msg)
{
#ifdef HAVE_FLITE
    cst_voice *v;
#endif

    fprintf(stderr, "%s\n", msg);

#ifdef HAVE_FLITE
    flite_init();

    v = register_cmu_us_kal();
    flite_text_to_speech(msg, v, "play");
#endif
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
   
    min_max(dset->x, &xmin, &xmax, dset->n);

    min_max(dset->y, &ymin, &ymax, dset->n);
    if (dset->y2 != NULL) {
	double y2min, y2max;

	min_max(dset->y2, &y2min, &y2max, dset->n);
	if (y2min < ymin) ymin = y2min;
	if (y2max > ymax) ymax = y2max;
    }

    xavg = (xmax - xmin) / dset->n;
    /* normalize average x step to quarter note */
    xscale = 1.0 / xavg;

    yscale = YMAX / (ymax - ymin);

    for (i=0; i<dset->n; i++) {
	double dtx, dux, ypos;

	ypos = (dset->y[i] - ymin) * yscale;

	if (i == 0) {
	    dtx = 0.0;
	} else {
	    dtx = xscale * (dset->x[i] - dset->x[i-1]);
	}

	if (i == dset->n - 1) {
	    dux = xscale * xavg;
	} else {
	    dux = xscale * (dset->x[i+1] - dset->x[i]);
	}

	track->notes[i].dtime = dtx;
	track->notes[i].duration = dux;
	track->notes[i].pitch = 36 + (int) (ypos + 0.5);
	track->notes[i].force = 64;

#ifdef DEBUG
	fprintf(stderr, "Obs %d: x = %g, y = %g\n", i, 
		dset->x[i], dset->y[i]);
	fprintf(stderr, " ypos = %g, dtx = %g, dux = %g\n",
		ypos, dtx, dux);
	fprintf(stderr, " dtime=%d, duration=%d, pitch=%d\n",
		track->notes[i].dtime, track->notes[i].duration,
		track->notes[i].pitch);
#endif
    }

    write_midi_track(track, spec);

    if (dset->y2 != NULL) {
	for (i=0; i<dset->n; i++) {
	    double ypos = (dset->y2[i] - ymin) * yscale;

	    track->notes[i].pitch = 36 + (int) (ypos + 0.5);
	}
	track->channel = 1;
	track->patch = PC_MARIMBA; 
	write_midi_track(track, spec);
    }

    free(track->notes);

#ifdef DEBUG
    fprintf(stderr, "freed track->notes, returning\n");
#endif

    return 0;
}

#ifdef GLIB2

static int audio_fork (const char *prog, const char *fname)
{
    gchar *argv[5];
    gboolean run;
    int i;
    
    argv[0] = g_strdup(prog);
    argv[1] = g_strdup("-A");
    argv[2] = g_strdup("600");
    argv[3] = g_strdup(fname);
    argv[4] = NULL;

    run = g_spawn_async(NULL, argv, NULL, G_SPAWN_SEARCH_PATH, 
                        NULL, NULL, NULL, NULL);

    for (i=0; i<5; i++) {
	g_free(argv[i]);
    }

    return !run;
}

#endif /* GLIB2 */

int midi_play_graph (const char *fname, const char *userdir)
{
    char outname[FILENAME_MAX];
#ifndef GLIB2
    char syscmd[1024];
#endif
    midi_spec spec;
    midi_track track;
    dataset dset;

    sprintf(outname, "%sgretl.mid", userdir);

    spec.fp = fopen(outname, "wb");
    if (spec.fp == NULL) {
	fprintf(stderr, "Couldn't write to '%s'\n", outname);
	return 1;
    }

    if (read_datafile(fname, &dset)) {
	audio_graph_error("Error reading data file");
	fclose(spec.fp);
	return 1;
    }

    spec.ntracks = (dset.y2 == NULL)? 1 : 2;
    spec.nticks = 96;
    spec.dset = &dset;
    spec.nsecs = 16;
    four_four_header(&spec);
    print_dataset_comments(&dset);
    play_dataset(&spec, &track, &dset);
    dataset_free(&dset);

    fclose(spec.fp);

#ifdef GLIB2
    audio_fork("timidity", outname);
#else
    sprintf(syscmd, "timidity -A 600 %s", outname);
    gretl_spawn(syscmd);
#endif

#ifdef DEBUG
    fprintf(stderr, "midi_play_graph, returning\n");
#endif

    return 0;
}
