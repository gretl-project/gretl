/* (C) Copyright 1993,1994 by Carnegie Mellon University
 * All Rights Reserved.
 *
 * Permission to use, copy, modify, distribute, and sell this software
 * and its documentation for any purpose is hereby granted without
 * fee, provided that the above copyright notice appear in all copies
 * and that both that copyright notice and this permission notice
 * appear in supporting documentation, and that the name of Carnegie
 * Mellon University not be used in advertising or publicity
 * pertaining to distribution of the software without specific,
 * written prior permission.  Carnegie Mellon University makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied
 * warranty.
 *
 * CARNEGIE MELLON UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO
 * THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS, IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY BE LIABLE
 * FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN
 * AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING
 * OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
 * SOFTWARE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

#include "mpack.h"
#include "md5.h"

static char basis_64[] =
   "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static void output64chunk (int c1, int c2, int c3, int pads, FILE *outfile)
{
    putc(basis_64[c1>>2], outfile);
    putc(basis_64[((c1 & 0x3)<< 4) | ((c2 & 0xF0) >> 4)], outfile);

    if (pads == 2) {
        putc('=', outfile);
        putc('=', outfile);
    } else if (pads) {
        putc(basis_64[((c2 & 0xF) << 2) | ((c3 & 0xC0) >>6)], outfile);
        putc('=', outfile);
    } else {
        putc(basis_64[((c2 & 0xF) << 2) | ((c3 & 0xC0) >>6)], outfile);
        putc(basis_64[c3 & 0x3F], outfile);
    }
}

static int to64 (FILE *infile, FILE *outfile)
{
    int c1, c2, c3, ct=0, written=0;

    while ((c1 = getc(infile)) != EOF) {
        c2 = getc(infile);
        if (c2 == EOF) {
            output64chunk(c1, 0, 0, 2, outfile);
        } else {
            c3 = getc(infile);
            if (c3 == EOF) {
                output64chunk(c1, c2, 0, 1, outfile);
            } else {
                output64chunk(c1, c2, c3, 0, outfile);
            }
        }
        ct += 4;
        if (ct > 71) {
            putc('\n', outfile);
	    written += 73;
            ct = 0;
        }
    }

    if (ct) {
	putc('\n', outfile);
	ct++;
    }

    return written + ct;
}

static void md5contextTo64 (MD5_CTX *context, char *encodedDigest)
{
    unsigned char digest[18];
    int i;
    char *p;

    MD5Final(digest, context);
    digest[sizeof(digest)-1] = digest[sizeof(digest)-2] = 0;

    p = encodedDigest;

    for (i=0; i < sizeof(digest); i+=3) {
	*p++ = basis_64[digest[i]>>2];
	*p++ = basis_64[((digest[i] & 0x3)<<4) | ((digest[i+1] & 0xF0)>>4)];
	*p++ = basis_64[((digest[i+1] & 0xF)<<2) | ((digest[i+2] & 0xC0)>>6)];
	*p++ = basis_64[digest[i+2] & 0x3F];
    }

    *p-- = '\0';
    *p-- = '=';
    *p-- = '=';
}    

static void md5digest (FILE *infile, char *digest)
{
    MD5_CTX context;
    unsigned char buf[1000];
    int nbytes;
    
    MD5Init(&context);

    while ((nbytes = fread(buf, 1, sizeof(buf), infile))) {
	MD5Update(&context, buf, nbytes);
    }

    rewind(infile);

    md5contextTo64(&context, digest);
}

/* Encode a file into a MIME message */

int mpack_encode (FILE *fpin, const char *fname, const char *note, 
		  const char *subject, const char *recipient, 
		  const char *reply_to, const char *type, 
		  FILE *fpout)
{
    const char *cleanfname, *p;
    char digest[25];

    /* Clean up fname for printing */
    cleanfname = fname;
    if ((p = strrchr(cleanfname, '/')) != NULL) cleanfname = p+1;
    if ((p = strrchr(cleanfname, '\\')) != NULL) cleanfname = p+1;
    if ((p = strrchr(cleanfname, ':')) != NULL) cleanfname = p+1;

    /* Compute MD5 digests */
    md5digest(fpin, digest);

    fprintf(fpout, "Mime-Version: 1.0\n");
    fprintf(fpout, "From: %s\n", reply_to);
    fprintf(fpout, "To: %s\n", recipient);
    fprintf(fpout, "Subject: %s\n", subject);
    fputs("Content-Type: multipart/mixed; boundary=\"-\"\n", fpout);
    fputs("\nThis is a MIME encoded message.\n\n", fpout);

    /* description section */
    if (note != NULL) {
	fputs("---\n\n", fpout);
	fputs(note, fpout);
	fputc('\n', fpout);
    }

    fputs("---\n", fpout);

    fprintf(fpout, "Content-Type: %s; name=\"%s\"\n", type,
	    cleanfname);
    fputs("Content-Transfer-Encoding: base64\n", fpout);
    fprintf(fpout, "Content-Disposition: inline; filename=\"%s\"\n", 
	    cleanfname);

    fprintf(fpout, "Content-MD5: %s\n\n", digest);

    to64(fpin, fpout);

    fputs("\n-----\n", fpout);

    return 0;
}
