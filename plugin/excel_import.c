/*
  gretl importer plugin for MS Excel files, Allin Cottrell, 2002.
  Based on xls2csv (Copyright 1998 David Rysdam)
  as distributed in the "catdoc" package by Vitus Wagner.
  This file is released under the GPL.  Details can be
  found in the file COPYING accompanying this distribution.
*/

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "xltypes.h"
#include "catdoc.h"
#include <stdlib.h>
#include <unistd.h>

/************************************************************************/
/* Displays  help message                                               */
/************************************************************************/

void help (void) {
    printf("Usage:\n xls2csv [-xl] [-s charset] [-d charset] [-c char] [ -q number] files\n");
}

void print_sheet(void);
void free_sheet(void);
void print_value(char *value);
int getshort(char *rec,int offset);
char *format_double(char *rec,int offset,int format_code);
char *gettypename(long rectype);
void printhexy(unsigned char rec[], long reclen);
void process_item (int rectype, int reclen, char *rec); 
void allocate(int row,int col);
char *copy_unicode_string(char **src);
char *convert8to8(char *src,int count);
char *convert16to8(char *src,int count);
void do_table(FILE *input,char *filename);
char *mark_string(char *instr);

short int *source_charset = NULL;
CHARSET target_charset;
SUBSTMAP spec_chars;

/* Defines unicode chars which should be
   replaced by strings before UNICODE->target chatset
   mappigs are applied i.e. TeX special chars like %
*/

SUBSTMAP replacements;
char *input_buffer, *output_buffer;

#define QUOTE_NEVER 0
#define QUOTE_SPACES_ONLY 1
#define QUOTE_ALL_STRINGS 2
#define QUOTE_EVERYTHING 3

int quote_mode = QUOTE_ALL_STRINGS;
char cell_separator = ',';


void do_table (FILE *input,char *filename) {    
    long rectype;
    long reclen;
    int eof_flag=0;
    unsigned char rec[MAX_MS_RECSIZE];
    int itemsread=1;

    while (itemsread) {
	fread(rec, 2, 1, input);
	if (rec[0] != 9 || rec[1] != 8) {
	    itemsread=fread(rec, 126, 1, input);
	} else {
	    fread(rec, 2, 1, input);
	    reclen = getshort(rec, 0);
	    if (reclen == 8 || reclen == 16) {
		itemsread = fread(rec, reclen, 1, input);
		break;
	    } else {
		fprintf(stderr, "%s: Invalid BOF record\n", filename);
	        return;
	    } 
	}
    }    
    if (feof(input)) {
	fprintf(stderr, "%s: No BOF record found\n", filename);
	exit(1);
    }  
   
    while (itemsread) {
	char buffer[2];

	rectype = 0;
	itemsread = fread(buffer, 2, 1, input);
	rectype = getshort(buffer,0);
	if (itemsread == 0)
	    break;
	reclen = 0;

	itemsread = fread(buffer, 2, 1, input);
	reclen = getshort(buffer,0);
	if (reclen && reclen < MAX_MS_RECSIZE && reclen > 0){
	    itemsread = fread(rec, 1, reclen, input);
	    rec[reclen] = '\0';
	}
	if(eof_flag) {
	    if (rectype != BOF) {
		break;
	    }    
	}
	process_item(rectype, reclen, rec);
	if (rectype == MSEOF) {
	    eof_flag = 1;
	} else {
	    eof_flag = 0;	
	}  
    }

    fclose(input);
    return;
}

struct rowdescr {
    int last, end;
    char **cells;
};	

char **sst;
int sstsize = 0, sstnext = 0;
int codepage = 1251; /*default*/
struct rowdescr *rowptr = NULL;
char **saved_reference = NULL;
int startrow = 0, lastrow = 0;

void process_item (int rectype, int reclen, char *rec) {
    switch (rectype) {
    case CODEPAGE: {
	char *source_csname;	
		   
	if (source_charset) break;
	codepage = getshort(rec,0);
	if (codepage != 1200) {
	    char cp[9];

	    sprintf(cp, "cp%d", codepage);
	    check_charset(&source_csname, cp);
	    source_charset = read_charset(source_csname);
	}		 
	break;
    }  
    case FORMAT: 
	break;
    case SST: {
	char *ptr = rec + 8, **outptr;
	int i;

	sstsize = getshort(rec, 4);
	sst=malloc(sstsize * sizeof(char *));
	for (i=0, outptr=sst; i<sstsize && (ptr-rec)<reclen; i++, outptr++) {
	    *outptr = copy_unicode_string(&ptr);
	} 
	if (i < sstsize) {
	    sstnext = i;
	}	  
	break;
    }	
    case CONTINUE: {
	int i = sstnext;
	char **outptr, *ptr = rec;
	if (i == 0) {
	    break;
	}
	for (outptr=sst+i; i<sstsize && (ptr-rec)<reclen; i++, outptr++) {
	    *outptr = copy_unicode_string(&ptr);
	} 
	if (i < sstsize) {
	    sstnext = i;
	}	  
	break;
    }			   
    case LABEL: {
	int row, col, len;
	struct rowdescr *prow;

	saved_reference = NULL;
	row = getshort(rec, 0); 
	col = getshort(rec, 2);
	len = getshort(rec, 6);
	allocate(row, col);
	prow = rowptr + row;
	prow->cells[col] = mark_string(convert8to8(rec+8, len));
	break;
    }     
    case CONSTANT_STRING: {
	struct rowdescr *prow;
	int row = getshort(rec, 0); 
	int col = getshort(rec, 2);
	int string_no = getshort(rec, 6);

	saved_reference = NULL;
	allocate(row, col);
	prow = rowptr + row;
	if (string_no >= sstsize) {
	    fprintf(stderr, "string index too large\n");
	} else if (sst[string_no] != NULL) {	
	    int len = strlen(sst[string_no]);
	    char *outptr;

	    outptr = prow->cells[col] = malloc(len + 2);
	    *(outptr++) = '"';
	    strcpy(outptr, sst[string_no]);
	} else {
	    prow->cells[col] = malloc(1);
	    strcpy(prow->cells[col], "");
	}	
	break;
    }
    case NUMBER: {
	int row, col;
	struct rowdescr *prow;

	saved_reference = NULL;
	row = getshort(rec, 0) - startrow; 
	col = getshort(rec, 2);
	allocate(row, col);
	prow = rowptr + row;
	prow->cells[col] = strdup(format_double(rec, 6, getshort(rec, 4)));
	break;
    }
    case FORMULA: { 
	int row, col;
	struct rowdescr *prow;

	saved_reference = NULL;
	row = getshort(rec, 0) - startrow; 
	col = getshort(rec, 2);
	allocate(row, col);
	prow = rowptr + row;
	if (((unsigned char) rec[12] == 0xFF) && (unsigned char) rec[13] == 0xFF) {
	    /* not a floating point value */
	    if (rec[6] == 1) {
		/*boolean*/
		char buf[2] = "0";

		buf[0] += rec[9];
		prow->cells[col] = strdup(buf);
	    } else if (rec[6] == 2) {
		/*error*/
		char buf[6] = "ERROR";

		prow->cells[col] = strdup(buf);
	    } else if (rec[6] == 0) {
		saved_reference = prow->cells + col;
	    }   
	} else {
	    int format_code = getshort(rec, 4);

	    prow->cells[col] = strdup(format_double(rec, 6, format_code));
	}		 
	break;
    }
    case STRING: {
	int len;

	if (!saved_reference) {
	    fprintf(stderr, "String record without preceeding string formula\n");
	    break;
	}
	len = getshort(rec,0);
	*saved_reference = mark_string(convert8to8(rec + 2, len + 1));
	break;
    }	    
    case BOF: {
	if (rowptr) {
	    fprintf(stderr, "BOF when current sheet is not flushed\n");
	    free_sheet();
	}
	break;
    }	  
    case MSEOF: {
	if (!rowptr) break;
	print_sheet();
	free_sheet();
	break;
    }
    default: {
    }
    }
}  

char *mark_string (char *instr) {
    int len = strlen(instr);
    char *out = malloc(len+2);

    *out = '"';
    strcpy(out+1, instr);
    free(instr);
    return out;

}    
char *copy_unicode_string (char **src) {
    int count = getshort(*src, 0);
    int flags = *(*src + 2);
    int to_skip = 0;
    int charsize = (flags & 0x01)? 2 : 1;
    char *realstart;

    if (flags & 0x04) {
	/* extended string */
	if (flags & 0x08) {
	    /*rich string*/;
	    to_skip= 3 + 4 + (charsize*count) + 2 + 4 * getshort(*src, 3) +
		getshort(*src, 5);
	} else {
	    to_skip= 3 + 4 + (charsize * count) + getshort(*src, 3);
	}
	fprintf(stderr, "Extended string found length=%d charsize=%d %d bytes to skip\n",
		count, charsize, to_skip);
	*src += to_skip;
	return strdup("Extended string");
    } else {
	to_skip = 3 + charsize * count + ((flags & 0x8)? (2+4*getshort(*src,3)) : 0);
    }
    realstart = *src + ((flags & 0x8)? 6 : 3);
    *src += to_skip;
    if (charsize == 1) {
	return convert8to8(realstart, count);
    } else {
	return convert16to8(realstart, count);
    }   
}

char *convert8to8 (char *src,int count) {
    char *dest = malloc(count+1);
    char *s, *d, *c;
    int i, u, len, l;

    len = count;
    *dest = 0; l = 0;
    for (s=src, d=dest, i=0; i<count; i++, s++) {
	if (source_charset) {
	    u = to_unicode(source_charset, (unsigned char) *s);
	} else {
	    u = (unsigned char) *s;
	}	 
	c = convert_char(u);
	l += strlen(c);
	while (l > len) {
	    len += 16;
	    fprintf(stderr, "realloc:count was %d len=%d\n", count, len);
	    dest = realloc(dest, len + 1);
	    d = dest + l;
	}
	strcpy(d, c);
	d += strlen(c);
    }   
    return dest;

}

char *convert16to8 (char *src,int count) {
    char *dest = malloc(count+1);
    char *s, *d, *c;
    int i, u, len, l;

    len = count;
    *dest = 0; l = 0;
    for (s=src, d=dest, i=0; i<count; i++, s+=2) {
	u = getshort(s, 0);
	c = convert_char(u);
	l += strlen(c);
	while (l > len) {
	    len += 16;
	    fprintf(stderr, "realloc:count was %d len=%d\n", count, len);
	    dest = realloc(dest, len+1);
	    d = dest + l;
	}
	strcpy(d, );
	d += trlen(c);
    }      
    return dest;
}

void allocate (int row,int col) {
    int newrow, newcol;

    if (row >= lastrow) {
	newrow = (row/16+1)*16;
	rowptr = realloc(rowptr,newrow*sizeof(struct rowdescr));
	memset(rowptr+lastrow,0,(newrow-lastrow)*sizeof(struct rowdescr));
	lastrow = newrow;
    }
    if (col >= rowptr[row].end) {
	newcol = (col/16+1)*16;
	rowptr[row].cells = realloc(rowptr[row].cells,newcol*sizeof(char *));
	memset(rowptr[row].cells+rowptr[row].end,0,(newcol-rowptr[row].end)
	       *sizeof(char *));
	rowptr[row].end = newcol;
    }  
    if (col>rowptr[row].last) rowptr[row].last = col;
}

char *format_double (char *rec,int offset,int format_code) {	
    static char buffer [128];
    union { 
	char cc[8];
	double d;
    } dconv;
    char *d, *s;
    int i;

# ifdef WORDS_BIGENDIAN     
    for(s=rec+offset+8, d=dconv.cc, i=0; i<8; i++) *(d++)=*(--s);
# else       
    for(s=rec+offset, d=dconv.cc, i=0; i<8; i++) *(d++)=*(s++);
# endif     
    sprintf(buffer, "%.10g", dconv.d);
    return buffer;
}

int getshort(char *rec,int offset) {
    return (signed short int) (*((unsigned char *)(rec+offset)) |
			       ((*((unsigned char *)(rec+offset+1)))<<8));
}	      

void free_sheet (void) {
    int i,j;
    struct rowdescr *row;
    char **col;

    for (row = rowptr, i=0; i<lastrow; i++, row++) {
	if (!row->cells) continue;
	for (col = row->cells, j=0; j<row->end; j++, col++) {
	    if (*col) {
		free(*col);
	    }
	}
	free(row->cells);
    }
    free(rowptr);
    rowptr = NULL;
    lastrow = 0;
}

void print_sheet (void) {
    int i,j;
    struct rowdescr *row;
    char **col;

    lastrow--;
    while (lastrow>0&&!rowptr[lastrow].cells) lastrow--;
    for(i=0,row=rowptr;i<=lastrow;i++,row++) {
	if (row->cells) {
	    for (j=0,col=row->cells;j<=row->last;j++,col++) {
		if (j) fputc(cell_separator,stdout);
		if (*col) {
		    print_value(*col);
		}
	    }
	}
	fputc('\n',stdout);
    }
} 

time_t float2date (double f) { 
    /* Hacked version. Excell stores date as floating point count of days
     * since 1.1.1900. We are substracting value of 1.1.1970 and multiplying
     * by 86400 thus getting seconds from the epoch
     */
    return rint((f-25569.0)*86400); 
}

void print_value (char *value) 
{
    int i,len;
    int quotes=0;
    int is_string=(*value == '"');
    char *ptr = is_string? value+1 : value;

    len = strlen(ptr);
    switch (quote_mode) {
    case QUOTE_NEVER:
	break;
    case QUOTE_SPACES_ONLY:	  
	for (i=0; i<len; i++) {
	    if (isspace(ptr[i]) || ptr[i] == cell_separator || 
		ptr[i] == '"') {
		quotes=1;
		break;
	    }
	}    
	break;
    case QUOTE_ALL_STRINGS:
	quotes = is_string;
	break;
    case QUOTE_EVERYTHING:
	quotes = 1;
	break;     
    }  	  
    if (quotes) {
	fputc('\"', stdout);
	for (i=0; i<len; i++) {
	    if (ptr[i] == '\"') {
		fputc('\"', stdout);
		fputc('\"', stdout);
	    } else {
		fputc(ptr[i], stdout);
	    }
	}   
	fputc('\"', stdout);
    } else {
	fputs(ptr, stdout);
    }
}    

int main (int argc, char *argv[])
{
    FILE *input;
    char *filename;
    short int *tmp_charset;
    char *dest_csname; 
    int c;
    int i;
 
    check_charset(&dest_csname, TARGET_CHARSET); 

    while ((c = getopt(argc, argv, "ls:d:xq:c:")) != -1) {
	switch (c) {
	case 'l':
	    list_charsets(); exit(0);
	case 'x': 
	    unknown_as_hex = 1; break;
	case 's':
	    check_charset(&source_csname, optarg);
	    source_charset = read_charset(source_csname);
	    break;
	case 'd':
	    check_charset(&dest_csname, optarg);
	    break;
	case 'q':
	    {   
		char *errptr;
	    
		quote_mode = strtol(optarg, &errptr, 0);
		if ((errptr && *errptr) || quote_mode < 0 || quote_mode > 3) {
		    fprintf(stderr,
			    "argument of -q should be number from 0 to 3\n");
		    exit(1);
		}
	    }    
	    break;
	case 'c':
	    cell_separator = optarg[0];
	    break;
	default:
	    help();
	    exit(1);
	}	
    }

    /* charset conversion init */
    input_buffer = malloc(FILE_BUFFER);
    tmp_charset = read_charset(dest_csname);
    target_charset = make_reverse_map(tmp_charset);
    free(tmp_charset);

    spec_chars = read_substmap(stradd("ascii", SPEC_EXT));
    if (!spec_chars) {
	fprintf(stderr, "Cannod read substitution map ascii%s\n",
		SPEC_EXT);
	exit(1);
    }  
    replacements = read_substmap(stradd("ascii", REPL_EXT));
    if (!replacements) {
        fprintf(stderr, "Cannod read substitution map ascii%s\n",
		REPL_EXT);
	exit(1);
    }  
    if (optind >= argc) {
	if (isatty(fileno(stdin))) {
	    help();
	    exit(0);
	}    
	do_table(stdin, "STDIN");
        exit (0);
    }	
    for (i=optind; i<argc; i++) {
	filename = argv[i];
	input = fopen(filename, "rb");
	if (!input) {
	    perror(filename);
	    exit(1);
	}   
        do_table(input, filename);
    }	
    return 0;
}
