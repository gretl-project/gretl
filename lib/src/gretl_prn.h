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
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 *
 */

#ifndef GRETL_PRN_H
#define GRETL_PRN_H

typedef enum {
    GRETL_PRINT_STDOUT,
    GRETL_PRINT_STDERR,
    GRETL_PRINT_FILE,
    GRETL_PRINT_BUFFER,
    GRETL_PRINT_NULL
} PrnType;

typedef enum {
    GRETL_PRINT_FORMAT_PLAIN,
    GRETL_PRINT_FORMAT_TEX,
    GRETL_PRINT_FORMAT_TEX_DOC,
    GRETL_PRINT_FORMAT_RTF
} PrnFormat;

/* functions follow */

void gretl_print_destroy (PRN *prn);

PRN *gretl_print_new (PrnType ptype);

PRN *gretl_print_new_with_filename (const char *fname);

PRN *gretl_print_new_with_buffer (char *buf);

int gretl_print_reset_buffer (PRN *prn);

const char *gretl_print_get_buffer (PRN *prn);

void gretl_print_set_format (PRN *prn, PrnFormat format);

int pprintf (PRN *prn, const char *template, ...);

int pputs (PRN *prn, const char *s);

int pputc (PRN *prn, int c);

void gretl_print_flush_stream (PRN *prn);

int printing_to_standard_stream (PRN *prn);

int printing_is_redirected (PRN *prn);

int print_start_redirection (PRN *prn, FILE *fp);

int print_end_redirection (PRN *prn);

int plain_format (PRN *prn);

int rtf_format (PRN *prn);

int tex_format (PRN *prn);

int tex_doc_format (PRN *prn);

#endif /* GRETL_PRN_H */
