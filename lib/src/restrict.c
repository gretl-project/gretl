#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define E_ALLOC 4

enum {
    OP_PLUS,
    OP_MINUS
};

typedef struct restriction_ restriction;

struct restriction_ {
    int nterms;
    char *op;
    int *coeff;
    double rhs;
};

int n_restrictions;
restriction **restrictions;

static int count_ops (const char *p)
{
    int n = 0 ;

    while (*p) {
	if (*p == '+' || *p == '-') n++;
	if (*p == '=') break;
	p++;
    }

    return n;
}

static int parse_chunk (const char *s)
{
    int bnum = -1;

    sscanf(s, " b%d", &bnum);

    return bnum;
}

static void restriction_destroy (restriction *r)
{
    free(r->op);
    free(r->coeff);
    free(r);
}

static restriction *restriction_new (int n)
{
    restriction *r;

    r = malloc(sizeof *r);
    if (r == NULL) return NULL;

    r->op = NULL;
    r->coeff = NULL;
	
    r->op = malloc(n * sizeof *r->op);
    r->coeff = malloc(n * sizeof *r->coeff);
    if (r->op == NULL || r->coeff == NULL) {
	restriction_destroy(r);
	return NULL;
    }

    r->nterms = n;
    r->rhs = 0.0;

    return r;
}

static void print_restriction (restriction *r)
{
    int i;

    printf("Number of LHS terms = %d\n", r->nterms);
    for (i=0; i<r->nterms; i++) {
	if (r->op[i] == OP_PLUS) printf(" + ");
	else printf(" - ");
	printf("b%d", r->coeff[i]);
    }
    printf(" = %g\n", r->rhs);
}

static int 
add_to_restriction (restriction *r, char op, int bnum, int i)
{
    r->op[i] = op;
    r->coeff[i] = bnum;

    return 0;
}

int gretl_parse_restriction_line (const char *line)
{
    char op = OP_PLUS;
    const char *p = line;
    restriction *r;
    double x;
    int i, nt, err = 0;

    if (!strncmp(line, "restrict", 8)) p += 8;
    while (isspace((unsigned char) *p)) p++;

    if (*p == '+') p++;
    else if (*p == '-') {
	op = OP_MINUS;
	p++;
    }
	
    nt = 1 + count_ops(p);

    r = restriction_new(nt);
    if (r == NULL) return E_ALLOC;

    for (i=0; i<nt; i++) {
	char chunk[8];
	int len, bnum;

	len = strcspn(p, "+-=");
	if (len > 7) {
	    err = 1;
	    break;
	}
	*chunk = 0;
	strncat(chunk, p, len);
	p += len;

	bnum = parse_chunk(chunk);
	if (bnum < 0) {
	    err = 1;
	    break;
	}

	add_to_restriction(r, op, bnum, i);

	if (*p == '+') {
	    op = OP_PLUS;
	    p++;
	} else if (*p == '-') {
	    op = OP_MINUS;
	    p++;
	}
    }

    if (!err) {
	if (!sscanf(p, " = %lf", &x)) err = 1;
    }

    if (err) {
	restriction_destroy(r);
    } else {
	r->rhs = x;
	print_restriction(r);
    }
    
    return err;
}

#if 0

int gretl_restrict_init (const char *line)
{
    int err = 0;

    if (restrictions != NULL) {
	free_restrictions(restrictions);
    }

    n_restrictions = 0;

    restrictions = malloc(sizeof *restrictions);
    if (restrictions == NULL) return E_ALLOC;

    return err;
}

#endif

int main (void)
{
    const char *foo = "restrict b1 + b2 = 0";
    int err = 0;

    err = gretl_parse_restriction_line(foo);
    if (err) {
	printf("Got %d from gretl_parse_restriction_line\n", err);
    } 

    return 0;
}
