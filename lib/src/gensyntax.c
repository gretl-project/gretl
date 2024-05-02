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

/* parser module for 'genr' and related commands */

#include "genparse.h"
#include "uservar_priv.h"
#include "gretl_string_table.h"
#include "genr_optim.h"

#if GENDEBUG
# define SDEBUG GENDEBUG
# define MDEBUG GENDEBUG
#else
# define SDEBUG 0
# define MDEBUG 0
#endif

static NODE *powterm (parser *p, NODE *l);

#if SDEBUG
static void notify (const char *s, NODE *n, parser *p)
{
    if (n == NULL) {
	fprintf(stderr, "%-8s: returning NULL node, err = %d\n",
		s, p->err);
    } else {
	fprintf(stderr, "%-8s: returning node of type %d (%s) at %p, err = %d\n",
		s, n->t, getsymb(n->t), (void *) n, p->err);
    }
}
#endif

NODE *new_node (int t)
{
    NODE *n = calloc(1, sizeof *n);

#if MDEBUG
    fprintf(stderr, "new_node: allocated node of type %d (%s) at %p\n",
	    t, getsymb(t), (void *) n);
#endif

    if (n != NULL) {
	n->t = t;
	n->vnum = NO_VNUM;
	n->L = n->M = n->R = NULL;
	n->parent = NULL;
    }

    return n;
}

NODE *newempty (void)
{
    NODE *n = new_node(EMPTY);

    if (n != NULL) {
	n->v.idnum = 0;
    }

    return n;
}

static NODE *newref (parser *p, int t)
{
    NODE *n = new_node(t);

    if (n != NULL) {
	if (t == SERIES) {
	    n->vnum = p->idnum;
	    n->v.xvec = p->dset->Z[n->vnum];
	    n->vname = p->idstr;
	} else if (t == NUM || t == NUM_P || t == NUM_M) {
	    user_var *u = p->data;

	    n->vname = p->idstr;
	    n->v.xval = *(double *) u->ptr;
	    n->uv = u;
	} else if (t == MAT || t == LIST || t == BUNDLE ||
		   t == ARRAY || t == STR) {
	    user_var *u = p->data;

	    n->vname = p->idstr;
	    n->v.ptr = u->ptr;
	    n->uv = u;
	} else if (t == OSL) {
	    user_var *u = p->data;

	    n->vname = p->idstr;
	    n->v.ptr = u->ptr;
	    n->uv = u;
	    n->flags |= LHT_NODE;
	} else if (t == UFUN) {
	    n->vname = p->idstr; /* function name */
	    n->v.ptr = p->data;  /* pointer to function */
        } else if (t == RFUN) {
            n->vname = p->idstr; /* function name */
	} else if (t == DBUNDLE) {
	    n->v.idnum = p->idnum;
	} else if (t == MMEMB) {
	    /* for case of named model */
	    n->vname = p->idstr;
	    n->v.idnum = -1;
	    n->t = DBUNDLE; /* roll into DBUNDLE case */
	} else if (t == UNDEF) {
	    n->vname = p->idstr;
	} else if (t == UOBJ || t == WLIST) {
	    n->v.str = p->idstr;
	} else {
	    n->v.idnum = p->idnum;
	}
    }

    return n;
}

/* node storing an anonymous string value */

static NODE *newstr (char *s)
{
    NODE *n = NULL;

    if (s == NULL) {
	fprintf(stderr, "newstr: input is NULL\n");
    } else {
	n = new_node(STR);
	if (n != NULL) {
	    n->v.str = s;
	    n->flags = TMP_NODE;
	}
    }

    return n;
}

/* node storing an anonymous floating point value */

NODE *newdbl (double x)
{
    NODE *n = new_node(NUM);

    if (n != NULL) {
	n->v.xval = x;
    }

    return n;
}

/* node for binary operator or two-argument function */

NODE *newb2 (int t, NODE *l, NODE *r)
{
    NODE *n = new_node(t);

    if (n != NULL) {
	n->L = l;
	n->R = r;
    }

    return n;
}

/* node for ternary operator */

static NODE *newb3 (int t, NODE *l)
{
    NODE *n = new_node(t);

    if (n != NULL) {
	n->L = l;
    }

    return n;
}

/* node for unknown number of subnodes */

static NODE *newbn (int t)
{
    NODE *n = new_node(t);

    if (n != NULL) {
	n->v.bn.n_nodes = 0;
	n->v.bn.n = NULL;
    }

    return n;
}

NODE *bncopy (NODE *t, int *err)
{
    NODE *n = new_node(t->t);

    if (n != NULL) {
	int i, k = t->v.bn.n_nodes;
	NODE **nn = malloc(k * sizeof *nn);

	if (nn == NULL) {
	    *err = E_ALLOC;
	} else {
	    n->v.bn.n_nodes = k;
	    n->v.bn.n = nn;
	    for (i=0; i<k; i++) {
		n->v.bn.n[i] = t->v.bn.n[i];
	    }
	}
    }

    return n;
}

static int push_bn_node (NODE *t, NODE *n)
{
    NODE **nn;
    int k = t->v.bn.n_nodes;

    nn = realloc(t->v.bn.n, (k + 1) * sizeof *nn);
    if (nn == NULL) {
	return E_ALLOC;
    }

    t->v.bn.n = nn;
    t->v.bn.n[k] = n;
    t->v.bn.n_nodes += 1;
    n->parent = t;

#if SDEBUG
    fprintf(stderr, "push_bn_node: n_nodes now = %d, "
	    "added node of type %d (vnum = %d)\n",
	    t->v.bn.n_nodes, n->t, n->vnum);
#endif

    return 0;
}

static void expected_symbol_error (int c, parser *p, int badc)
{
    parser_print_input(p);

    if (badc) {
	gretl_errmsg_sprintf(_("Expected '%c' but found '%c'\n"), c,
			     badc);
    } else {
	const char *found = getsymb_full(p->sym, p);

	if (found == NULL || *found == '\0') {
	    gretl_errmsg_sprintf(_("Expected '%c' but formula ended\n"), c);
	} else {
	    gretl_errmsg_sprintf(_("Expected '%c' but found '%s'\n"), c,
				 found);
	}
    }

    p->err = E_PARSE;
}

static void unmatched_symbol_error (int c, parser *p)
{
    if (c == '(' || c == ')') {
	fprintf(stderr, "*** gensyntax: unmatched_symbol_error '%c'\n", c);
    } else {
	fprintf(stderr, "*** gensyntax: unmatched_symbol_error (%c)\n", c);
    }
    parser_print_input(p);
    pprintf(p->prn, _("Unmatched '%c'\n"), c);
    p->err = E_PARSE;
}

static NODE *base (parser *p, NODE *up)
{
    NODE *t = NULL;

    if (p->err) {
	return NULL;
    }

#if SDEBUG
    fprintf(stderr, "base(): on input sym = %d ('%s'), ch = '%c'\n",
	    p->sym, getsymb_full(p->sym, p), p->ch? p->ch : '0');
#endif

    switch (p->sym) {
    case CNUM:
	t = newdbl(p->xval);
	lex(p);
	break;
    case CSTR:
	t = newstr(p->idstr);
	lex(p);
	break;
    case NUM:
    case SERIES:
    case MAT:
    case BUNDLE:
    case ARRAY:
    case STR:
    case UOBJ:
    case CON:
    case DVAR:
    case MVAR:
    case LIST:
    case WLIST:
    case DBUNDLE:
    case NUM_P:
    case NUM_M:
    case OSL:
    case UNDEF:
	t = newref(p, p->sym);
	lex(p);
	break;
    case DUM:
	if (p->idnum == DUM_NULL || p->idnum == DUM_EMPTY) {
	    t = newempty();
	} else {
	    t = newref(p, p->sym);
	}
	lex(p);
	break;
    case EMPTY:
    case NULLARG:
	t = newempty();
	if (p->sym == NULLARG) {
	    t->vname = p->idstr;
	}
	lex(p);
	break;
    case B_RANGE:
    case P_DOT:
	lex(p);
	t = expr(p);
	break;
    case B_JOIN:
	/* list joiner with empty LHS */
	t = newempty();
	break;
    case G_LPR: /* left paren '(' */
	lex(p);
	t = expr(p);
	if (p->sym == G_RPR) {
	    lex(p);
	} else if (p->err == 0) {
	    expected_symbol_error(')', p, 0);
	}
	break;
    case G_LBR: /* left bracket '[' */
	if (up == NULL) {
	    goto deferr;
	}
	if (up->t == OBS) {
	    t = obs_node(p);
	}
	if (p->sym == G_RBR) {
	    lex(p);
	} else if (p->err == 0) {
	    expected_symbol_error(']', p, 0);
	}
	break;
    default:
    deferr:
	context_error(0, p, "base");
	break;
    }

#if SDEBUG
    notify("base", t, p);
    fprintf(stderr, "on exit from base, p->sym = %s (p->err = %d)\n",
	    getsymb(p->sym), p->err);
#endif

    return t;
}

/* Special for the unary '&' operator: the operand must
   be a named series or "user_var", with an additional
   restriction on the type.
*/

static NODE *u_addr_base (parser *p)
{
    NODE *t = base(p, NULL);

    if (t != NULL && !p->err) {
	switch (t->t) {
	case SERIES:
	    if (t->v.xvec == NULL) {
		p->err = E_TYPES;
	    }
	    break;
	case NUM:
	case MAT:
	case BUNDLE:
	case ARRAY:
	case STR:
	case LIST:
	    if (t->uv == NULL) {
		p->err = E_TYPES;
	    }
	    break;
	case OSL:
	case UNDEF:
	    break;
	case EMPTY:
	    if (t->vname == NULL) {
		p->err = E_TYPES;
	    }
	    break;
	default:
	    p->err = E_TYPES;
	}

	if (p->err) {
	    pputs(p->prn, _("Wrong type of operand for unary '&'"));
	    pputc(p->prn, '\n');
	}
    }

    return t;
}

static void unwrap_string_arg (parser *p)
{
    int n = strlen(p->idstr);

    if (p->idstr[n-1] == '"') {
	p->idstr[n-1] = '\0';
    } else {
	unmatched_symbol_error('"', p);
    }
}

#define varargs_func(f) (f == F_PRINTF || f == F_SPRINTF || \
			 f == F_SSCANF)

#define testc(c) (c == '(' || c == '.' || c == '[')

/* Grab a string argument. Note: we have a mechanism in genlex.c for
   retrieving arguments that take the form of quoted string literals
   or names of string variables.  The special use of this function is
   to grab a literal string without requiring the user to wrap it in
   quotes; we use it only where we know the only acceptable argument
   is a string.  This function assumes that the argument in question
   is the last (or only) argument to a function: we look for a closing
   parenthesis and flag an error if we don't find one.

   Depending on the context we may or may not want to "consume" the
   trailing right paren that follows the string of interest; that's
   controlled by the second argument.

   Added 2013-08-25: The above is all very well, but it breaks the
   case where a function that returns a string is given as an
   argument, rather than a plain string variable or string literal --
   the function call was getting passed as a string literal. A block
   is now inserted to test for this case.
*/

static NODE *get_final_string_arg (parser *p, NODE *t, int sym,
				   int eat_last)
{
    const char *src = NULL;
    int n, wrapped = 0;
    int strvar = 0;

    while (p->ch == ' ') {
	parser_getc(p);
    }

    if (!varargs_func(sym) && sym != F_EXISTS) {
	/* Check for a nested function call (2013-08-25) or
	   bundle/array member (2015-09-25). Further fixes
	   applied July 2018.
	*/
	src = p->point - 1;
	n = gretl_namechar_spn(src);
	if (n > 0 && n < VNAMELEN && testc(src[n])) {
	    NODE *ret = NULL;
	    char tmp[VNAMELEN];
	    char c = src[n];

	    *tmp = '\0';
	    strncat(tmp, src, n);
	    src = NULL;
	    if (c == '(') {
		if (function_lookup(tmp) || get_user_function_by_name(tmp)) {
		    /* watch out: this is VERY finicky! */
		    /* was: ret = base(p, t); */
		    lex(p);
		    ret = expr(p);
		}
	    } else if (gretl_is_bundle(tmp)) {
		lex(p);
		ret = expr(p);
	    } else if (c == '[' && get_array_by_name(tmp)) {
		lex(p);
		ret = expr(p);
	    }
	    if (ret != NULL) {
		if (eat_last) {
		    /* consume trailing right paren */
		    /* note: was parser_getc(p); */
		    lex(p);
		}
		return ret;
	    }
	}
    }

    if (p->ch == '"') {
	wrapped = 1;
	if (!varargs_func(sym)) {
	    parser_getc(p);
	}
    }

    if (p->ch == ')') {
	if (wrapped) {
	    p->err = E_PARSE;
	    return NULL;
	}
	/* handle empty arg */
	p->idstr = gretl_strdup("");
    } else {
	char p0 = '(', p1 = ')';
	int i, paren = 1, quoted = wrapped;
	int close = -1, started = 0;
	const char *s = p->point;

	if (p->ch == p0) paren++;

	/* find length of string to closing paren */
	i = 0;
	while (*s) {
	    if (!quoted) {
		if (*s == p0) {
		    paren++;
		} else if (*s == p1) {
		    paren--;
		}
	    }
	    if (paren == 0) {
		close = i;
		break;
	    }
	    if (*s == '"') {
		quoted = !quoted;
	    }
	    s++;
	    i++;
	}

	if (close < 0) {
	    unmatched_symbol_error(p0, p);
	    return NULL;
	}

	for (i=0; i<=close; i++) {
	    if (!started && !isspace(p->ch)) {
		src = p->point - 1;
		p->idstr = gretl_strndup(src, close - i + 1);
		started = 1;
	    }
	    parser_getc(p);
	}
    }

    if (p->idstr == NULL) {
	p->err = E_ALLOC;
	return NULL;
    }

    tailstrip(p->idstr);
    if (eat_last) {
	/* consume trailing right paren */
	parser_getc(p);
    }
    lex(p);

    if (!varargs_func(sym)) {
	if (wrapped) {
	    unwrap_string_arg(p);
	} else if (sym != F_ISDISCR &&
		   sym != F_TYPEOF &&
		   sym != F_EXISTS) {
	    /* not quoted: give priority to string variables
	       unless we need the _names_ of string variables
	       rather then their content
	    */
	    p->data = get_user_var_of_type_by_name(p->idstr,
						   GRETL_TYPE_STRING);
	    if (p->data != NULL) {
		strvar = 1;
	    } else {
		char *s = get_built_in_string_by_name(p->idstr);

		if (s != NULL) {
		    free(p->idstr);
		    p->idstr = gretl_strdup(s);
		}
	    }
	}
    }

#if SDEBUG
    fprintf(stderr, "get_final_string_arg: '%s' (strvar=%d)\n",
	    p->idstr, strvar);
#endif

    if (p->err) {
	return NULL;
    } else if (strvar == 1) {
	return newref(p, STR);
    } else {
	return newstr(p->idstr);
    }
}

/* Here we're allowing for the possibility that an "fncall"
   argument is given in the form of a string variable. If the
   relevant portion of the statement contains any quotation
   marks that cannot be the case.
*/

static void try_for_string_var (parser *p, int close, int anyquote)
{
    int err = 0;

    if (!anyquote) {
	char *sname = gretl_strndup(p->point - 1, close + 1);
	char *sval = get_string_by_name(sname);

	if (sval != NULL) {
	    p->idstr = gretl_strdup(sval);
	} else {
	    err = 1;
	}
	free(sname);
    } else {
	err = 1;
    }

    if (err) {
	unmatched_symbol_error(')', p);
    }
}

enum {
    RIGHT_STR  = 1 << 0,
    MID_STR    = 1 << 1
};

/* get_literal_string_arg() is mostly used to retrieve a string
   that defines a function call, e.g. in BFGSmax(). It can also
   be used to trap a non-quoted string that functions as an
   "enum" value, as in the "vcat" flag for mpireduce().
   In the former case @opt will be 0, in the latter MID_STR.
*/

static NODE *get_literal_string_arg (parser *p, int opt)
{
    int close = 0;

    while (p->ch == ' ') {
	parser_getc(p);
    }

    if (p->ch == '"') {
	/* handle the case where the argument is quoted */
	close = parser_char_index(p, '"');
	if (close < 0) {
	    unmatched_symbol_error('"', p);
	    return NULL;
	}
	p->idstr = gretl_strndup(p->point, close);
	close++;
    } else if (opt == MID_STR) {
	/* handle the "enum"-type unquoted case */
	int i1 = parser_char_index(p, ',');
	int i2 = parser_char_index(p, ')');

	if (i2 < 0) {
	    unmatched_symbol_error('(', p);
	    return NULL;
	}
	close = (i1 >= 0 && i1 < i2)? i1 : i2;
	p->idstr = gretl_strndup(p->point - 1, close + 1);
    } else {
	/* handle the function-call case: the terminator of
	   the argument will be a bare comma (unquoted, not
	   in parentheses) or a right paren which matches
	   the paren that opens the argument list of
	   the function call.
	*/
	const char *s = p->point;
	int gotparen = 0;
	int quoted = 0;
	int anyquote = 0;
	int paren = 0;
	int i = 0;

	if (p->ch == '(') {
	    /* an "anonymous function"? */
	    gotparen = paren = 1;
	}

	while (*s) {
	    if (*s == '"') {
		anyquote = 1;
		quoted = !quoted;
	    }
	    if (!quoted) {
		if (*s == '(') {
		    gotparen++;
		    paren++;
		} else if (*s == ')') {
		    if (paren == 0) {
			paren--;
			close = i;
			break;
		    }
		    paren--;
		}
		if (paren == 0) {
		    if (gotparen) {
			/* include right paren */
			close = i + 1;
			break;
		    } else if (*s == ',') {
			/* don't include comma */
			close = i;
			break;
		    }
		}
	    }
	    s++;
	    i++;
	}

	if (paren > 0) {
	    unmatched_symbol_error('(', p);
	} else if (paren < 0) {
	    try_for_string_var(p, close, anyquote);
	} else if (quoted) {
	    unmatched_symbol_error('"', p);
	} else {
	    p->idstr = gretl_strndup(p->point - 1, close + 1);
	}
    }

    if (p->idstr == NULL) {
	if (!p->err) {
	    p->err = E_ALLOC;
	}
    } else {
	tailstrip(p->idstr);
	parser_advance(p, close);
	lex(p);
#if SDEBUG
	fprintf(stderr, "get_literal_string_arg: '%s'\n", p->idstr);
#endif
    }

    return (p->idstr == NULL)? NULL : newstr(p->idstr);
}

static NODE *get_bundle_member_name (parser *p, int dollarize)
{
    NODE *ret = NULL;

#if SDEBUG
    fprintf(stderr, "bundle_member_name: sym='%s', ch='%c', point='%s'\n",
	    getsymb(p->sym), p->ch, p->point);
#endif

    if (p->ch == '[') {
	/* using bracketed key notation */
	parser_getc(p); /* eat opening '[' */
	if (!p->err) {
	    lex(p); /* get next character */
	    ret = expr(p);
	}
	if (!p->err) {
	    if (p->sym != G_RBR) {
		unmatched_symbol_error('[', p);
	    } else {
		lex(p); /* eat closing ']' */
	    }
	}
    } else if (p->ch == '.') {
	/* using bundle dot notation */
	int i, n;

	if (dollarize && *p->point == '$') {
	    /* allow leading '$' for objects under named models */
	    n = 1 + gretl_namechar_spn(p->point + 1);
	} else {
	    /* otherwise should be standard valid identifier */
	    n = gretl_namechar_spn(p->point);
	}
	if (n == 0 || n >= VNAMELEN) {
	    p->err = E_PARSE;
	} else {
	    p->idstr = gretl_strndup(p->point, n);
	    if (p->idstr == NULL) {
		p->err = E_ALLOC;
	    } else {
		for (i=0; i<=n; i++) {
		    parser_getc(p);
		}
		lex(p);
		ret = newstr(p->idstr);
	    }
	}
    } else {
	fprintf(stderr, "HERE, get_bundle_member_name, p->ch = '%c'\n", p->ch);
	p->err = E_PARSE;
    }

    return ret;
}

/* Find the name of the (putative) series to the right of '.'
   in an expression of the form <listname>.<series-name>,
   and return it in a string node.
*/

static NODE *series_name_node (parser *p)
{
    int n = gretl_namechar_spn(p->point);
    NODE *ret = NULL;

    if (n == 0 || n >= VNAMELEN) {
	p->err = E_PARSE;
    } else {
	p->idstr = gretl_strndup(p->point, n);
	if (p->idstr == NULL) {
	    p->err = E_ALLOC;
	} else {
	    int i;

	    for (i=0; i<=n; i++) {
		parser_getc(p);
	    }
	    lex(p);
	    ret = newstr(p->idstr);
	}
    }

    return ret;
}

static void get_matrix_def (NODE *t, parser *p)
{
    NODE *n;
    char cexp = 0;

#if SDEBUG
    fprintf(stderr, "get_matrix_def, p->sym = %d ('%s')\n",
	    p->sym, getsymb(p->sym));
#endif

    if (p->sym == G_LCB) {
	lex(p);
	if (!p->err && p->sym == G_RCB) {
	    /* empty matrix def, OK */
	    lex(p);
	    return;
	}
	while (p->ch != 0 && !p->err) {
	    n = expr(p);
	    if (p->err) {
		break;
	    } else {
		p->err = push_bn_node(t, n);
	    }
	    if (p->sym == P_COM) {
		/* comma: on to the next column */
		lex(p);
	    } else if (p->sym == P_SEMI) {
		/* semicolon: on to the next row */
		n = newempty();
		p->err = push_bn_node(t, n);
		lex(p);
	    } else if (p->sym == G_RCB) {
		/* right curly bracket: reached the end */
		break;
	    } else {
		/* something that doesn't belong here */
		cexp = '}';
		break;
	    }
	}
    } else {
	cexp = '{';
    }

#if SDEBUG
    fprintf(stderr, " after processing, sym = %d, err = %d\n",
	    p->sym, p->err);
#endif

    if (!p->err && cexp == 0) {
	if (p->sym == G_RCB) {
	    lex(p);
	} else {
	    unmatched_symbol_error('{', p);
	}
    }

    if (cexp && !p->err) {
	expected_symbol_error(cexp, p, 0);
    }
}

#define set_slice_on(p) (p->flags |= P_SLICING)
#define set_slice_off(p) (p->flags &= ~P_SLICING)

#define set_lag_parse_on(p) (p->flags |= P_LAGPRSE)
#define set_lag_parse_off(p) (p->flags &= ~P_LAGPRSE)

static void get_slice_parts (NODE *t, parser *p)
{
    int slice_upsym = p->upsym;
    char cexp = 0;

#if SDEBUG
    fprintf(stderr, "get_slice_parts, p->sym = %d (%s), upsym %s\n",
	    p->sym, getsymb(p->sym), getsymb(p->upsym));
#endif

    set_slice_on(p);

    if (p->sym == G_LBR) {
	lex(p);
	if (p->sym == P_COM) {
	    /* empty row spec, OK */
	    t->L = newempty();
	} else {
	    t->L = expr(p);
	}
	if (p->sym == P_COL) {
	    /* second part of colon-separated range */
	    t->L = newb2(SUBSL, t->L, NULL);
	    lex(p);
	    if (p->sym == P_COM || p->sym == G_RBR) {
		/* reached end: second part implicitly empty */
		t->L->R = newempty();
	    } else {
		t->L->R = expr(p);
	    }
	}
	if (p->sym == G_RBR) {
	    /* no comma, no second arg string: may be OK */
	    t->R = NULL;
	    lex(p);
	    set_slice_off(p);
	    return;
	}
	if (p->sym == P_COM) {
	    p->upsym = slice_upsym; /* 2021-09-19 */
	    lex(p);
	    if (p->sym == G_RBR) {
		/* empty column spec, OK */
		t->R = newempty();
	    } else {
		t->R = expr(p);
	    }
	    if (p->sym == P_COL) {
		/* second part of colon-separated range */
		t->R = newb2(SUBSL, t->R, NULL);
		lex(p);
		if (p->sym == G_RBR) {
		    /* reached end: second part implicitly empty */
		    t->R->R = newempty();
		} else {
		    t->R->R = expr(p);
		}
	    }
	    if (p->sym == G_RBR) {
		lex(p);
	    } else {
		cexp = ']';
	    }
	} else {
	    cexp = ',';
	}
    } else {
	cexp = '[';
    }

    if (cexp && !p->err) {
	expected_symbol_error(cexp, p, 0);
    }

    set_slice_off(p);
}

/* contexts where empty arguments are never acceptable */
#define no_empty(f) (f == F_DEFBUNDLE || f == F_DEFARRAY || f == F_DEFLIST)

static void attach_child (NODE *parent, NODE *child, int f,
			  int np, int i, parser *p)
{
    if (p->err) {
	return;
    }

#if SDEBUG
    fprintf(stderr, "attach_child: np=%d, i=%d, type '%s'\n", np, i,
	    getsymb(child->t));
#endif

    /* catch erroneous cases */
    if (np > 0 && i == np) {
	gretl_errmsg_sprintf("%s: %s", getsymb_full(f, p),
			     _("too many arguments"));
	p->err = E_INVARG;
    } else if (child->t == EMPTY && no_empty(f)) {
	p->err = E_PARSE;
    }

    if (p->err) {
	free_tree(child, p, 0);
	return;
    }

    if (np == 1) {
	/* 1-place node */
	parent->L = child;
	child->parent = parent;
    } else if (np == 2) {
	/* 2-place node */
	if (i == 0) {
	    parent->L = child;
	} else {
	    parent->R = child;
	}
	child->parent = parent;
    } else if (np == 3) {
	/* 3-place node */
	if (i == 0) {
	    parent->L = child;
	} else if (i == 1) {
	    parent->M = child;
	} else {
	    parent->R = child;
	}
	child->parent = parent;
    } else if (child != NULL) {
	/* n-place node */
	p->err = push_bn_node(parent, child);
	child->parent = parent;
    }
}

static void pad_parent (NODE *parent, int np, int i, parser *p)
{
    int j;

    for (j=i; j<np; j++) {
	attach_child(parent, newempty(), 0, np, j, p);
    }
}

struct argrecord {
    int f;
    int vec[4];
};

/* Here we record the positions of function-call arguments to
   built-in functions that take such arguments. These require
   special handling: they should be passed forward as string
   literals even if they are not wrapped in double-quotes.
*/

struct argrecord fncall_argrec[] = {
    {F_BFGSMAX,  {0, 1, 1, 0}},
    {F_NRMAX,    {0, 1, 1, 1}},
    {F_FDJAC,    {0, 1, 0, 0}},
    {F_SIMANN,   {0, 1, 0, 0}},
    {F_BFGSCMAX, {0, 0, 1, 1}},
    {F_NMMAX,    {0, 1, 0, 0}},
    {F_GSSMAX,   {0, 1, 0, 0}},
    {F_NUMHESS,  {0, 1, 0, 0}},
    {F_FZERO,    {1, 0, 0, 0}},
};

static const int *get_callargs (int f)
{
    int i, n = G_N_ELEMENTS(fncall_argrec);

    for (i=0; i<n; i++) {
	if (f == fncall_argrec[i].f) {
	    return fncall_argrec[i].vec;
	}
    }

    if (fncall_func(f)) {
	/* inconsistent with info in genparse.h */
	fprintf(stderr, "ERROR: function %s has an 'fncall' arg"
		" but no argrecord\n", getsymb(f));
    }

    return NULL;
}

static int next_arg_is_string (int i, const int *callargs, int k,
			       int opt)
{
    if (callargs && callargs[i]) {
	return 1;
    }
    if ((opt & MID_STR) && i > 0 && i < k-1) {
	return 1;
    }
    if ((opt & RIGHT_STR) && i == k-1) {
	return 1;
    }
    return 0;
}

/* Get up to @np comma-separated arguments (possibly optional).
   However, if np < 0 this is a signal to get as many arguments as we
   can find (the number unknown in advance).
*/

static void get_args (NODE *t, parser *p, int f, int np,
		      int opt, int *next)
{
    NODE *child;
    const int *callargs;
    char cexp = 0;
    int i = 0;

#if SDEBUG
    fprintf(stderr, "get_args: f=%s, np=%d, ch='%c', point='%s'\n",
	    getsymb(f), np, p->ch, p->point);
#endif

    if (p->sym != G_LPR) {
	expected_symbol_error('(', p, 0);
	return;
    }

    if (f == F_STACK) {
	/* for stack(), the first arg should be a list */
	p->flags |= P_LISTDEF;
    }

    callargs = get_callargs(f);
    if (callargs == NULL || callargs[0] == 0) {
	/* nothing special about the first arg */
	lex(p);
    }

    if (f == F_STACK) {
	/* revert to regular parsing */
	p->flags &= ~P_LISTDEF;
    }

    while (p->sym != G_RPR && p->ch != '\0' && !p->err) {
	/* get the next argument */
	if (p->sym == P_COM) {
	    /* implies an empty argument slot */
	    child = newempty();
	} else if (i < 4 && callargs && callargs[i]) {
	    /* a function-call argument: don't insist on quotation */
	    child = get_literal_string_arg(p, 0);
	} else if (i > 0 && i < np - 1 && (opt & MID_STR)) {
	    child = get_literal_string_arg(p, opt);
	} else if (i == np - 1 && (opt & RIGHT_STR)) {
	    child = get_final_string_arg(p, t, f, 0);
	} else {
	    child = expr(p);
	}

	if (!p->err) {
	    attach_child(t, child, f, np, i++, p);
	}

	if (p->err || p->sym == G_RPR) {
	    break;
	} else if (p->sym == P_COM) {
	    /* turn off flag for accepting string as first arg */
	    p->flags &= ~P_GETSTR;
	    if (next_arg_is_string(i, callargs, np, opt)) {
		/* the next arg needs special handling */
		p->sym = EMPTY;
	    } else {
		lex(p);
		if (p->sym == G_RPR) {
		    /* trailing empty argument slot */
		    attach_child(t, newempty(), f, np, i++, p);
		}
	    }
	} else {
	    /* either arg-separating comma or closing paren was expected */
	    cexp = (i < np)? ',' : ')';
	    break;
	}
    }

    if (!p->err) {
	if (cexp) {
	    expected_symbol_error(cexp, p, 0);
	} else if (p->sym == G_RPR) {
	    if (i < np) {
		pad_parent(t, np, i, p);
	    }
	    lex(p);
	    if (p->sym == G_LBR) {
		*next = '[';
	    }
	} else {
	    expected_symbol_error(')', p, 0);
	}
    }

#if SDEBUG
    fprintf(stderr, "get_args: returning with p->err=%d\n", p->err);
#endif
}

static void get_assertion (NODE *t, parser *p)
{
    char *str, *s = strrchr(p->point, ')');
    int n = s - p->point + 1;
    int next = 0;

    str = gretl_strndup(p->point - 1, n);

    get_args(t, p, F_ASSERT, 1, 0, &next);
    if (!p->err) {
	str = g_strchomp(g_strchug(str));
	t->R = newstr(str);
    }
}

/* For defining a bundle via _(): get one or more comma-separated
   terms: each must take the form key=value, or just key if the
   key and the name of a pre-defined object are one and the same.
*/

static void get_bundle_pairs (NODE *t, parser *p, int *next)
{
    NODE *child;
    const char *src;
    int n, j, i = 0;

#if SDEBUG
    fprintf(stderr, "get_bundle_pairs: ch='%c', point='%s'\n",
	    p->ch, p->point);
#endif

    /* allow for an empty definition */
    while (p->ch == ' ') parser_getc(p);
    if (p->ch == ')') {
	parser_getc(p);
	lex(p);
	return;
    }

    while (p->ch && !p->err) {
	/* first get an unquoted key */
	while (p->ch == ' ') parser_getc(p);
	src = p->point - 1;
	n = gretl_namechar_spn(src);
	if (n == 0) {
	    p->err = E_PARSE;
	    break;
	}
	for (j=0; j<n; j++) {
	    parser_getc(p);
	}
	p->idstr = gretl_strndup(src, n);
	child = newstr(p->idstr);
	attach_child(t, child, 0, -1, i++, p);
	while (p->ch == ' ') parser_getc(p);
	if (p->ch == '=') {
	    /* parse the following expression */
	    parser_getc(p);
	    while (p->ch == ' ') parser_getc(p);
	    lex(p);
	    child = expr(p);
	} else {
	    /* backtrack to get named object */
	    parser_advance(p, -(p->point - src));
	    lex(p);
	    child = base(p, NULL);
	}
	if (!p->err) {
	    attach_child(t, child, 0, -1, i++, p);
	}
	if (p->sym == P_COM) {
	    ; /* OK */
	} else if (p->sym == G_RPR) {
	    lex(p);
	    break;
	} else {
	    gretl_errmsg_set(_("Expected ',' or ')'"));
	    p->err = E_PARSE;
	}
    }
}

static NODE *powterm (parser *p, NODE *l)
{
    int sym = p->sym;
    int opt = OPT_NONE;
    int next = 0;
    NODE *t = NULL;

    if (p->err) {
	return NULL;
    }

#if SDEBUG
    fprintf(stderr, "powterm, starting: p->sym = %d ('%s'), p->ch = '%c' (%d)\n",
	    p->sym, getsymb(p->sym), p->ch? p->ch : '0', p->ch);
#endif

    if (l != NULL) {
	/* powterm recursion: swallowing prior node @l */
	if (sym == BMEMB || sym == DBMEMB) {
	    t = newb2(sym, l, NULL);
	    if (t != NULL) {
		parser_ungetc(p);
		t->R = get_bundle_member_name(p, 0);
	    }
	} else if (sym == G_LBR) {
	    /* "OSL": we're being somewhat agnostic here regarding
	       the type of object of which we're taking a slice;
	       That will be sorted out at eval() time.
	    */
	    t = newb2(OSL, l, NULL);
	    if (t != NULL) {
		t->R = new_node(SLRAW);
		if (t->R != NULL) {
		    get_slice_parts(t->R, p);
		}
	    }
	} else if (sym == G_LPR) {
	    /* (bundle member or listvar) plus lag spec? */
	    if (l->t == BMEMB || l->t == LISTVAR) {
		t = newb2(LAG, l, NULL);
		if (t != NULL) {
		    t->R = base(p, NULL);
		}
	    } else {
		context_error('(', p, "powterm");
	    }
	}
	goto maybe_recurse;
    }

    if (string_last_func(sym)) {
	opt |= RIGHT_STR;
    }
    if (string_mid_func(sym)) {
	opt |= MID_STR;
    }

    if (sym == U_ADDR) {
        t = new_node(sym);
        if (t != NULL) {
            lex(p);
	    t->L = u_addr_base(p);
        }
    } else if (func2_symb(sym)) {
	int unset = 0;

	if (sym == F_GENSERIES) {
	    set_doing_genseries(1);
	    unset = 1;
	}
	t = new_node(sym);
	if (t != NULL) {
	    lex(p);
	    if (sym == F_ASSERT) {
		get_assertion(t, p);
	    } else {
		get_args(t, p, sym, 2, opt, &next);
	    }
	}
	if (unset) {
	    set_doing_genseries(0);
	}
    } else if (func3_symb(sym)) {
	t = newb3(sym, NULL);
	if (t != NULL) {
	    lex(p);
	    if (str0_func(sym)) {
		p->flags |= P_GETSTR;
	    }
	    get_args(t, p, sym, 3, opt, &next);
	}
    } else if (str0_func(sym)) {
	t = new_node(sym);
	if (t != NULL) {
	    lex(p);
	    t->L = newbn(FARGS);
	    if (t != NULL) {
		p->flags |= P_GETSTR;
		get_args(t->L, p, sym, -1, opt, &next);
	    }
	}
    } else if (string_arg_func(sym)) {
	t = new_node(sym);
	if (t != NULL) {
	    lex(p);
	    t->L = get_final_string_arg(p, t, sym, 1);
	}
    } else if (func1_symb(sym)) {
	if (undef_arg_ok(sym)) {
	    set_parsing_query(1);
	}
	t = new_node(sym);
	if (t != NULL) {
	    if (sym < FP_MAX) {
		t->v.ptr = p->data; /* attach function pointer */
	    }
	    lex(p);
	    get_args(t, p, sym, 1, opt, &next);
	}
	if (undef_arg_ok(sym)) {
	    set_parsing_query(0);
	}
    } else if (sym == LAG || sym == OBS) {
	if (sym == LAG) {
	    set_lag_parse_on(p);
	}
	t = new_node(sym);
	if (t != NULL) {
	    t->L = newref(p, p->upsym);
	    lex(p);
	    t->R = base(p, t);
	}
	if (sym == LAG) {
	    set_lag_parse_off(p);
	}
    } else if (sym == OSL) {
	t = new_node(sym);
	if (t != NULL) {
	    t->L = newref(p, p->upsym);
	    t->R = new_node(SLRAW);
	    if (t->R != NULL) {
		lex(p);
		get_slice_parts(t->R, p);
	    }
	}
    } else if (sym == DMSTR) {
	t = new_node(sym);
	if (t != NULL) {
	    t->L = newref(p, MVAR);
	    lex(p);
	    t->R = get_final_string_arg(p, t, sym, 1);
	}
    } else if (sym == BMEMB) {
	if (p->data == NULL) {
	    p->err = E_PARSE;
	    return NULL;
	}
	t = new_node(sym);
	if (t != NULL) {
	    t->L = newref(p, BUNDLE);
	    t->R = get_bundle_member_name(p, 0);
	}
    } else if (sym == DBMEMB) {
	t = new_node(sym);
	if (t != NULL) {
	    t->L = newref(p, DBUNDLE);
	    t->R = get_bundle_member_name(p, 0);
	}
    } else if (sym == MMEMB) {
	t = new_node(DBMEMB);
	if (t != NULL) {
	    t->L = newref(p, MMEMB);
	    t->R = get_bundle_member_name(p, 1);
	}
    } else if (sym == LISTVAR) {
	t = new_node(LISTVAR);
	if (t != NULL) {
	    t->L = newref(p, LIST);
	    t->R = series_name_node(p);
	}
    } else if (sym == G_LPR) {
	/* parenthesized expression */
	t = base(p, NULL);
    } else if (sym == G_LCB) {
	/* explicit matrix definition */
	t = new_node(F_DEFMAT);
	if (t != NULL) {
	    t->L = newbn(FARGS);
	    if (t->L != NULL) {
		get_matrix_def(t->L, p);
	    }
	}
    } else if (sym == UFUN || sym == RFUN) {
        t = newref(p, sym);
        if (t != NULL) {
            lex(p);
            t->R = newbn(FARGS);
	    if (t->R != NULL) {
		get_args(t->R, p, sym, -1, opt, &next);
	    }
        }
    } else if (sym == F_DEFARGS) {
	t = new_node(sym);
	if (t != NULL) {
	    lex(p);
	    t->L = newbn(FARGS);
	    if (t->L != NULL) {
		get_bundle_pairs(t->L, p, &next);
	    }
	}
    } else if (funcn_symb(sym)) {
	t = new_node(sym);
	if (t != NULL) {
	    lex(p);
	    t->L = newbn(FARGS);
	    if (t->L != NULL) {
		int k = -1;

		if (sym == F_NRMAX ||
		    sym == F_BFGSCMAX ||
		    sym == F_MOVAVG) {
		    k = 4;
		}
		get_args(t->L, p, sym, k, opt, &next);
	    }
	}
    } else {
	t = base(p, NULL);
    }

    if (als_func(sym) && (p->flags & P_ALIASED)) {
	/* transfer flag to newly created node */
	t->flags |= ALS_NODE;
	p->flags ^= P_ALIASED;
    }

 maybe_recurse:

    if (!p->err) {
	if (t != NULL && next == '[') {
	    /* support func(args)[slice], etc. */
	    t = newb2(OSL, t, NULL);
	    if (t != NULL) {
		t->R = new_node(SLRAW);
		if (t->R != NULL) {
		    get_slice_parts(t->R, p);
		}
	    }
	} else if (t != NULL) {
	    if (p->sym == BMEMB || p->sym == DBMEMB) {
		/* these types can recurse */
		t = powterm(p, t);
	    } else if (p->sym == G_LBR || p->sym == G_LPR) {
		t = powterm(p, t);
	    }
	}
    }

#if SDEBUG
    notify("powterm", t, p);
#endif

    return t;
}

static int got_hex_val (parser *p)
{
    int c = *p->point;

    return p->ch == '0' && (c == 'x' || c == 'X');
}

/* Test for whether the ' symbol must represent the unary
   transposition operator rather than binary transpose-multiply,
   based on the following symbol, @t.
*/

#define must_be_unary(t) (t == EOT || \
                          t < OP_MAX || \
			  t == P_COM || \
                          t == G_RPR || \
	                  t == G_RBR ||	\
                          t == G_RCB || \
			  t == P_COL)

static NODE *factor (parser *p)
{
    int sym = p->sym == B_SUB ? U_NEG :
	p->sym == B_ADD ? U_POS : p->sym;
    NODE *t;

#if SDEBUG
    fprintf(stderr, "factor: starting...\n");
#endif

    if (sym == U_NEG && got_hex_val(p)) {
	gretl_errmsg_set("hexadecimal values must be unsigned");
	p->err = E_TYPES;
    }
    if (p->err) {
	return NULL;
    }

    if (unary_op(sym) && sym != U_ADDR) {
        t = new_node(sym);
        if (t != NULL) {
            lex(p);
	    t->L = factor(p);
        }
    } else {
	t = powterm(p, NULL);
	if (t != NULL) {
	    int upsym = p->sym;

	    while (!p->err && (p->sym == B_POW ||
			       p->sym == B_DOTPOW ||
			       p->sym == B_TRMUL)) {
		t = newb2(p->sym, t, NULL);
		if (t != NULL) {
		    lex(p);
		    if (upsym == B_TRMUL && must_be_unary(p->sym)) {
			/* dummy RHS for unary transpose */
#if SDEBUG
			fprintf(stderr, "factor: B_TRMUL must in fact be unary\n");
#endif
			t->R = newempty();
		    } else {
			t->R = factor(p);
		    }
		}
	    }
	}
    }

#if SDEBUG
    notify("factor", t, p);
#endif

    return t;
}

static NODE *term (parser *p)
{
    NODE *t;

#if SDEBUG
    fprintf(stderr, "term: starting...\n");
#endif

    if (p->err || (t = factor(p)) == NULL) {
	return NULL;
    }

    while (!p->err && (p->sym == B_MUL || p->sym == B_DIV ||
		       p->sym == B_LDIV || p->sym == B_MOD ||
		       p->sym == B_DOTMULT || p->sym == B_DOTDIV ||
		       p->sym == B_KRON)) {
	t = newb2(p->sym, t, NULL);
	if (t != NULL) {
	    lex(p);
	    t->R = factor(p);
	}
    }

#if SDEBUG
    notify("term", t, p);
#endif

    return t;
}

static NODE *expr4 (parser *p)
{
    NODE *t;

#if SDEBUG
    fprintf(stderr, "expr4: starting...\n");
#endif

    if (p->err || (t = term(p)) == NULL) {
	return NULL;
    }

    while (!p->err && (p->sym == B_ADD || p->sym == B_SUB ||
		       p->sym == B_DOTADD || p->sym == B_DOTSUB ||
		       p->sym == B_HCAT || p->sym == B_VCAT ||
		       p->sym == B_LCAT || p->sym == B_RANGE ||
		       p->sym == B_ELLIP || p->sym == B_JOIN)) {
	t = newb2(p->sym, t, NULL);
	if (t != NULL) {
	    lex(p);
	    t->R = term(p);
	}
    }

#if SDEBUG
    notify("expr4", t, p);
#endif

    return t;
}

static NODE *expr3 (parser *p)
{
    NODE *t;

#if SDEBUG
    fprintf(stderr, "expr3: starting...\n");
#endif

    if (p->err || (t = expr4(p)) == NULL) {
	return NULL;
    }

    while (!p->err && (p->sym == B_GT || p->sym == B_LT ||
		       p->sym == B_DOTGT || p->sym == B_DOTLT ||
		       p->sym == B_GTE || p->sym == B_LTE ||
		       p->sym == B_DOTGTE || p->sym == B_DOTLTE)) {
	t = newb2(p->sym, t, NULL);
	if (t != NULL) {
	    lex(p);
	    t->R = expr4(p);
	}
    }

#if SDEBUG
    notify("expr3", t, p);
#endif

    return t;
}

static NODE *expr2 (parser *p)
{
    NODE *t;

#if SDEBUG
    fprintf(stderr, "expr2: starting...\n");
#endif

    if (p->err || (t = expr3(p)) == NULL) {
	return NULL;
    }

    while (!p->err && (p->sym == B_EQ || p->sym == B_NEQ ||
		       p->sym == B_DOTEQ || p->sym == B_DOTNEQ)) {
	t = newb2(p->sym, t, NULL);
	if (t != NULL) {
	    lex(p);
	    t->R = expr3(p);
	}
    }

#if SDEBUG
    notify("expr2", t, p);
#endif

    return t;
}

static NODE *expr1 (parser *p)
{
    NODE *t;

#if SDEBUG
    fprintf(stderr, "expr1: starting...\n");
#endif

    if (p->err || (t = expr2(p)) == NULL) {
	return NULL;
    }

    while (!p->err && p->sym == B_AND) {
	t = newb2(p->sym, t, NULL);
	if (t != NULL) {
	    p->flags |= P_AND;
	    lex(p);
	    t->R = expr2(p);
	    p->flags &= ~P_AND;
	}
    }

#if SDEBUG
    notify("expr1", t, p);
#endif

    return t;
}

static NODE *expr0 (parser *p)
{
    NODE *t;

#if SDEBUG
    fprintf(stderr, "expr0: starting...\n");
#endif

    if (p->err || (t = expr1(p)) == NULL) {
	return NULL;
    }

    while (!p->err && p->sym == B_OR) {
	t = newb2(p->sym, t, NULL);
	if (t != NULL) {
	    p->flags |= P_OR;
	    lex(p);
	    t->R = expr1(p);
	    p->flags &= ~P_OR;
	}
    }

#if SDEBUG
    notify("expr0", t, p);
#endif

    return t;
}

NODE *expr (parser *p)
{
    NODE *t;

#if SDEBUG
    fprintf(stderr, "expr: starting...\n");
#endif

    if (p->err || (t = expr0(p)) == NULL) {
	return NULL;
    }

    while (!p->err && p->sym == QUERY) {
	t = newb3(p->sym, t);
	if (t != NULL) {
	    set_parsing_query(1);
	    lex(p);
	    if (!p->err) {
		t->M = expr(p);
		if (p->sym == P_COL) {
		    lex(p);
		    if (!p->err) {
			t->R = expr(p);
		    }
		} else {
		    expected_symbol_error(':', p, 0);
		}
	    }
	    set_parsing_query(0);
	}
    }

#if SDEBUG
    notify("expr", t, p);
#endif

    return t;
}
