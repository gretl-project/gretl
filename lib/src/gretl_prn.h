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

#ifndef GRETL_PRN_H
#define GRETL_PRN_H

/**
 * PRN:
 *
 * An opaque structure accessed only via gretl_print functions.
 */

typedef struct PRN_ PRN;

typedef enum {
    GRETL_PRINT_STDOUT,
    GRETL_PRINT_STDERR,
    GRETL_PRINT_FILE,
    GRETL_PRINT_BUFFER,
    GRETL_PRINT_TEMPFILE,
    GRETL_PRINT_STREAM,
    GRETL_PRINT_GZFILE
} PrnType;

typedef enum {
    GRETL_FORMAT_TXT       = 1 << 0,
    GRETL_FORMAT_TEX       = 1 << 1,
    GRETL_FORMAT_DOC       = 1 << 2,
    GRETL_FORMAT_RTF       = 1 << 3,
    GRETL_FORMAT_RTF_TXT   = 1 << 4,
    GRETL_FORMAT_EQN       = 1 << 5,
    GRETL_FORMAT_SELECTION = 1 << 6,
    GRETL_FORMAT_CSV       = 1 << 7,
    GRETL_FORMAT_TAB       = 1 << 8,
    GRETL_FORMAT_MODELTAB  = 1 << 9,
    GRETL_FORMAT_LANDSCAPE = 1 << 10,
    GRETL_FORMAT_HAS_MINUS = 1 << 11,
    GRETL_FORMAT_XML       = 1 << 12
} PrnFormat;

/* functions follow */

PRN *gretl_print_new (PrnType ptype, int *err);

void gretl_print_destroy (PRN *prn);

PRN *gretl_print_new_with_filename (const char *fname, int *err);

PRN *gretl_gzip_print_new (const char *fname, int comp_level,
			   int *err);

PRN *gretl_print_new_with_tempfile (int *err);

int gretl_print_has_tempfile (PRN *prn);

const char *gretl_print_get_tempfile_name (PRN *prn);

PRN *gretl_print_new_with_buffer (char *buf);

PRN *gretl_print_new_with_gchar_buffer (gchar *buf);

PRN *gretl_print_new_with_stream (FILE *fp);

void gretl_print_detach_stream (PRN *prn);

int gretl_print_reset_buffer (PRN *prn);

int gretl_print_rename_file (PRN *prn, const char *oldpath, 
			     const char *newpath);

const char *gretl_print_get_buffer (PRN *prn);

const char *gretl_print_get_trimmed_buffer (PRN *prn);

char *gretl_print_steal_buffer (PRN *prn);

int gretl_print_replace_buffer (PRN *prn, char *buf);

void gretl_print_get_size (PRN *prn, int *width, int *height);

int gretl_print_set_save_position (PRN *prn);

void gretl_print_unset_save_position (PRN *prn);

char *gretl_print_get_chunk (PRN *prn);

char *gretl_print_get_chunk_at (PRN *prn, int pos);

int gretl_print_tell (PRN *prn);

void gretl_print_set_format (PRN *prn, PrnFormat format);

void gretl_print_toggle_doc_flag (PRN *prn);

void gretl_print_set_has_minus (PRN *prn);

int gretl_print_has_minus (PRN *prn);

void gretl_print_set_delim (PRN *prn, char delim);

int pprintf (PRN *prn, const char *format, ...);

void pprintf2 (PRN *prn, const char *format, ...);

int pputs (PRN *prn, const char *s);

int pputc (PRN *prn, int c);

void gretl_print_ensure_vspace (PRN *prn);

void gretl_prn_newline (PRN *prn);

void gretl_print_flush_stream (PRN *prn);

void gretl_print_close_stream (PRN *prn);

int printing_to_standard_stream (PRN *prn);

int print_redirection_level (PRN *prn);

const char *print_redirection_filename (PRN *prn);

int print_redirected_at_level (PRN *prn, int level);

int print_start_redirection (PRN *prn, FILE *fp,
			     const char *fname,
			     const char *strvar);

int print_end_redirection (PRN *prn);

int plain_format (PRN *prn);

int rtf_format (PRN *prn);

int rtf_doc_format (PRN *prn);

int tex_format (PRN *prn);

int tex_doc_format (PRN *prn);

int tex_eqn_format (PRN *prn);

int csv_format (PRN *prn);

int prn_format (PRN *prn);

char prn_delim (PRN *prn);

int gretl_print_has_buffer (PRN *prn);

int gretl_print_alloc (PRN *prn, size_t s);

char *gretl_print_read_tempfile (PRN *prn, int *err);

#endif /* GRETL_PRN_H */
