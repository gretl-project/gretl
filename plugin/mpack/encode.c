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

static FILE *createnewfile (char *fname)
{
    FILE *ret = NULL;
    int fd;
     
#ifdef O_EXCL
    fd = open(fname, O_RDWR|O_CREAT|O_EXCL, 0644);
#else
    fd = open(fname, O_RDWR|O_CREAT|O_TRUNC, 0644);
#endif

    if (fd != -1) {
	ret = fdopen(fd, "w");
    }

    return ret;
}

/* Encode a file into a MIME message */

int encode (FILE *fpin, const char *fname, const char *note, 
	    const char *subject, const char *recipient, const char *reply_to,
	    const char *type, char *tmpfname)
{
    FILE *fpout;
    const char *cleanfname, *p;
    char digest[25];
    char buf[1024];

    /* Clean up fname for printing */
    cleanfname = fname;
    if ((p = strrchr(cleanfname, '/')) != NULL) cleanfname = p+1;
    if ((p = strrchr(cleanfname, '\\')) != NULL) cleanfname = p+1;
    if ((p = strrchr(cleanfname, ':')) != NULL) cleanfname = p+1;

    /* Compute MD5 digests */
    md5digest(fpin, digest);

    /* Open output file */
    fpout = createnewfile(tmpfname);
    if (!fpout) {
	perror(buf);
	return 1;
    }

    fprintf(fpout, "Mime-Version: 1.0\n");
    fprintf(fpout, "From: %s\n", reply_to);
    /* fprintf(fpout, "Reply-To: %s\n", reply_to); */
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
	
    fclose(fpout);    

    return 0;
}
