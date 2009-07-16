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

#if GENDEBUG
# define SDEBUG 1
# define MDEBUG 1
#else
# define SDEBUG 0
# define MDEBUG 0
#endif

#define set_transpose(n) (n->flags |= TRANSP_NODE)

static NODE *powterm (parser *p);

#if SDEBUG
static void notify (const char *s, NODE *t, parser *p)
{
    fprintf(stderr, "%-8s: returning node at %p (type %d), err = %d\n", 
	    s, (void *) t, (t != NULL)? t->t : 0, p->err);
}
#endif

static NODE *newempty (int t)
{  
    NODE *n = malloc(sizeof *n);

#if MDEBUG
    fprintf(stderr, "newempty: allocated node at %p\n", (void *) n);
#endif

    if (n != NULL) {
	n->t = t;
	n->v.idnum = 0;
	n->flags = 0;
	n->vnum = NO_VNUM;
    }

    return n;
}

static NODE *newref (parser *p, int t)
{  
    NODE *n = malloc(sizeof *n);

#if MDEBUG
    fprintf(stderr, "newref: allocated node at %p (type %d)\n", 
	    (void *) n, t);
#endif

    if (n != NULL) {
	if (t == USERIES) {
	    n->vnum = p->idnum;
	    n->t = VEC;
	    n->v.xvec = (*p->Z)[n->vnum];
	} else {
	    n->t = t;
	    if (t == USCALAR || t == UMAT || t == UOBJ || t == LIST) {
		n->v.str = p->idstr;
	    } else {
		n->v.idnum = p->idnum;
	    }
	    n->vnum = NO_VNUM;
	}
	n->flags = 0;
    }

    return n;
}

static NODE *newstr (parser *p, int t)
{  
    NODE *n;

    if (p->idstr == NULL) {
	fprintf(stderr, "newstr: input is NULL\n");
	return NULL;
    }

    n = malloc(sizeof *n);

#if MDEBUG
    fprintf(stderr, "newstr: allocated node at %p (s = '%s')\n", 
	    (void *) n, p->idstr);
#endif

    if (n != NULL) {
	n->t = t;
	n->v.str = p->idstr;
	n->flags = 0;
	n->vnum = NO_VNUM;
    }

    return n;
}

/* node storing a floating point value */

NODE *newdbl (double x)
{  
    NODE *n = malloc(sizeof *n);

#if MDEBUG
    fprintf(stderr, "newdbl: allocated node at %p (x = %g)\n", 
	    (void *) n, x);
#endif

    if (n != NULL) {
	n->t = NUM;
	n->v.xval = x;
	n->flags = 0;
	n->vnum = NO_VNUM;
    }

    return n;
}

/* node for unary operator, or single-argument function */

static NODE *newb1 (int t, NODE *b)
{  
    NODE *n = malloc(sizeof *n);

#if MDEBUG
    fprintf(stderr, "newb1:  allocated node at %p (type = %d)\n", 
	    (void *) n, t);
#endif

    if (n != NULL) {
	n->t = t;
	n->v.b1.b = b;
	n->flags = 0;
	n->vnum = NO_VNUM;
    }

    return n;
}

/* node for binary operator or two-argument function */

static NODE *newb2 (int t, NODE *l, NODE *r)
{  
    NODE *n = malloc(sizeof *n);

#if MDEBUG
    fprintf(stderr, "newb2:  allocated node at %p (type = %d)\n", 
	    (void *) n, t);
#endif

    if (n != NULL) {
	n->t = t;
	n->v.b2.l = l;
	n->v.b2.r = r;
	n->flags = 0;
	n->vnum = NO_VNUM;
    }

    return n;
}

/* node for ternary operator */

static NODE *newb3 (int t, NODE *l)
{  
    NODE *n = malloc(sizeof *n);

#if MDEBUG
    fprintf(stderr, "newb3:  allocated node at %p (type = %d)\n", 
	    (void *) n, t);
#endif

    if (n != NULL) {
	n->t = t;
	n->v.b3.l = l;
	n->v.b3.m = NULL;
	n->v.b3.r = NULL;
	n->flags = 0;
	n->vnum = NO_VNUM;
    }

    return n;
}

/* node for unknown number of subnodes */

static NODE *newbn (int t)
{  
    NODE *n = malloc(sizeof *n);

#if MDEBUG
    fprintf(stderr, "newbn:  allocated node at %p (type = %d)\n", 
	    (void *) n, t);
#endif

    if (n != NULL) {
	n->t = t;
	n->v.bn.n_nodes = 0;
	n->v.bn.n = NULL;
	n->flags = 0;
	n->vnum = NO_VNUM;
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

#if SDEBUG
    fprintf(stderr, "push_bn_node: n_nodes now = %d, "
	    "added node of type %d (vnum = %d)\n", 
	    t->v.bn.n_nodes, n->t, n->vnum);
#endif

    return 0;
}

static void expected_symbol_error (int c, parser *p)
{
    const char *found = getsymb(p->sym, p);

    parser_print_input(p);

    if (*found == '\0') {
	pprintf(p->prn, _("Expected '%c' but formula ended\n"), c);
    } else {
	pprintf(p->prn, _("Expected '%c' but found '%s'\n"), c, 
		found);
    }

    if (!strcmp(found, "&")) {
	pputs(p->prn, "(for logical AND, please use \"&&\")\n");
    } else if (!strcmp(found, "|")) {
	pputs(p->prn, "(for logical OR, please use \"||\")\n");
    }

    p->err = 1;
}

static void unmatched_symbol_error (int c, parser *p)
{
    parser_print_input(p);
    pprintf(p->prn, _("Unmatched '%c'\n"), c);
    p->err = E_PARSE;
}

#define matrix_ref_node(p) (p->sym == UMAT || (p->sym == MVAR && \
                            model_data_matrix(p->idnum)))

/* try to recognize a unary ', indicating that we're taking
   the transpose of the preceding matrix */

static int unary_apost (parser *p) 
{
    if (p->ch == '\'') {
	/* peek ahead */
	int c = parser_next_nonspace_char(p);

	if (isalpha(c) || c == '$' || c == '(') {
	    /* next symbol is variable or expression: the
	       apostrophe must be the binary operator, 
	       B_TRMUL
	    */
	    return 0;
	} else {
	    return 1;
	}
    } 

    return 0;
}

static NODE *base (parser *p, NODE *up)
{ 
    NODE *t = NULL;

    if (p->err) {
	return NULL;
    }

#if SDEBUG
    fprintf(stderr, "base(): on input sym = %d ('%s'), ch = '%c'\n", 
	    p->sym, getsymb(p->sym, p), p->ch? p->ch : '0');
#endif

    switch (p->sym) {
    case NUM: 
	t = newdbl(p->xval);
	lex(p);
	break;
    case USCALAR: 
    case USERIES:
    case UMAT:
    case UOBJ:
    case CON: 
    case DVAR:
    case MVAR:
    case OBS:
    case LIST:
	t = newref(p, p->sym);
	if (matrix_ref_node(p) && unary_apost(p)) {
	    set_transpose(t);
	    parser_getc(p);
	}
	lex(p);
	break;
    case DUM:
	if (p->idnum == DUM_NULL) {
	    t = newempty(EMPTY);
	} else {
	    t = newref(p, p->sym);
	}
	lex(p);
	break;
    case B_SUB:
    case B_ADD:
    case B_RANGE:
    case P_DOT:
	lex(p);
	t = expr(p);
	break;
    case G_LPR: /* left paren '(' */
	lex(p);
	t = expr(p);
	if (p->sym == G_RPR) {
	    if (up != NULL && up->t != LAG && unary_apost(p)) {
		set_transpose(up);
		parser_getc(p);
	    } 
	    lex(p);
	} else if (p->err == 0) {
	    expected_symbol_error(')', p);
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
	    if (up->t == MSL || up->t == DMSL) {
		if (p->ch == '\'') {
		    set_transpose(up);
		    parser_getc(p);
		}
	    }
	    lex(p);
	} else if (p->err == 0) {
	    expected_symbol_error(']', p);
	}
	break;
    default: 
    deferr:
	context_error(0, p);
	break;
    }

#if SDEBUG
    notify("base", t, p);
#endif    

    return t;
}

/* construct a node that contains a reference to a data series, based
   on the name of a list (already captured in p->idstr) and the name
   of a variable that is supposed to be in the list, which is separated
   from the listname by '.'.
*/

static NODE *listvar_node (parser *p)
{
    int *list = get_list_by_name(p->idstr);
    NODE *ret = NULL;
    char vname[VNAMELEN] = {0};
    int i, n, v, ok = 0;

    if (list == NULL) {
	p->err = E_UNKVAR;
	return NULL;
    }

    n = gretl_namechar_spn(p->point);
    if (n >= VNAMELEN) {
	/* too long -- can't be a valid varname */
	p->err = E_UNKVAR;
	return NULL;
    }

    for (i=0; i<n; i++) {
	parser_getc(p);
	vname[i] = p->ch;
    }

    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (!strcmp(vname, p->dinfo->varname[v])) {
	    /* found the variable */
	    free(p->idstr);
	    p->idstr = NULL;
	    p->idnum = v;
	    ret = newref(p, USERIES);
	    if (ret == NULL) {
		p->err = E_ALLOC;
	    } else {
		ok = 1;
	    }
	    break;
	}
    }

    if (!ok && !p->err) {
	p->err = E_UNKVAR;
    }

    if (!p->err) {
	parser_getc(p);
	lex(p);
    }

    return ret;    
}

/* Grab a string argument.  Note: we have a mechanism in genlex.c for
   retrieving arguments that take the form of quoted string literals
   or names of string variables.  The special use of this function is
   to grab a literal string without requiring the user to wrap it in
   quotes; we use it only where we know the only acceptable argument
   is a string.  This function assumes that the argument in question
   is the last (or only) argument to a function: we look for a closing
   parenthesis and flag an error if we don't find one.
*/

static NODE *get_string_arg (parser *p)
{
    int wrapped = 0;

    while (p->ch == ' ') {
	parser_getc(p);
    }

    if (p->ch == '"') {
	wrapped = 1;
	parser_getc(p);
    }

    if (p->ch == ')') {
	/* allow empty arg string "()" */
	p->idstr = gretl_strdup("");
    } else {
	int i, paren = 1, quoted = wrapped;
	int close = -1, started = 0;
	const char *s = p->point;

	/* find length of string to closing paren */
	i = 0;
	while (*s) {
	    if (!quoted && *s == ')') paren--;
	    if (paren == 0) {
		close = i;
		break;
	    }
	    if (*s == '"') {
		quoted = !quoted;
	    } else if (!quoted && *s == '(') {
		paren++;
	    }
	    s++;
	    i++;
	}

	if (close < 0) {
	    unmatched_symbol_error('(', p);
	    return NULL;
	}

	for (i=0; i<=close; i++) {
	    if (!started && !isspace(p->ch)) {
		p->idstr = gretl_strndup(p->point - 1, close - i + 1);
		started = 1;
	    }
	    parser_getc(p);
	}
    }

    if (p->idstr == NULL) {
	p->err = E_ALLOC;
	return NULL;
    }

    parser_getc(p);
    lex(p);

    tailstrip(p->idstr);

    if (wrapped) {
	int n = strlen(p->idstr);

	if (p->idstr[n-1] == '"') {
	    p->idstr[n-1] = '\0';
	} else {
	    unmatched_symbol_error('"', p);
	}
    }

#if SDEBUG
    fprintf(stderr, "get_string_arg: '%s'\n", p->idstr);
#endif

    return newstr(p, STR);
}

enum {
    RIGHT_STR = 1
};

static void get_matrix_def (NODE *t, parser *p, int *sub)
{
    NODE *n;
    char cexp = 0;

#if SDEBUG
    fprintf(stderr, "get_matrix_def, p->sym = %d\n", p->sym);
#endif    

    if (p->sym == G_LCB) {
	lex(p);
	while (p->ch != 0 && !p->err) {
	    n = expr(p);
	    if (p->err) {
		break;
	    } else {
		p->err = push_bn_node(t, n);
	    }
	    if (p->sym == P_COM) {
		lex(p);
	    } else if (p->sym == P_SEMI) {
		n = newempty(EMPTY);
		p->err = push_bn_node(t, n);
		lex(p);
	    } else if (p->sym == G_RCB) {
		if (p->ch == '\'') {
		    set_transpose(t);
		    parser_getc(p);
		} else if (p->ch == '[') {
		    parser_ungetc(p);
		    *sub = 1;
		}
		break;
	    }
	}
    } else {
	cexp = '{';
    }

    if (cexp == 0) {
	if (p->sym == G_RCB) {
	    lex(p);
	} else {
	    unmatched_symbol_error('{', p);
	}
    }
	    
    if (cexp && !p->err) {
	expected_symbol_error(cexp, p);
    }
}	

static void get_slice_parts (NODE *t, parser *p)
{
    char cexp = 0;

#if SDEBUG
    fprintf(stderr, "get_slice_parts, p->sym = %d\n", p->sym);
#endif  

    set_matrix_slice_on();

    if (p->sym == G_LBR) {
	lex(p);
	if (p->sym == P_COM) {
	    /* empty row spec, OK */
	    t->v.b2.l = newempty(EMPTY);
	} else {
	    t->v.b2.l = expr(p);
	}
	if (p->sym == P_COL) {
	    /* second part of colon-separated range */
	    t->v.b2.l = newb2(SUBSL, t->v.b2.l, NULL);
	    lex(p);
	    if (p->sym == P_COM || p->sym == G_RBR) {
		/* reached end: second part implicitly empty */
		t->v.b2.l->v.b2.r = newempty(EMPTY);
	    } else {
		t->v.b2.l->v.b2.r = expr(p);
	    }
	}
	if (p->sym == G_RBR) {
	    /* no comma, no second arg string: may be OK */
	    t->v.b2.r = newempty(ABSENT);
	    lex(p);
	    set_matrix_slice_off();
	    return;
	}
	if (p->sym == P_COM) {
	    lex(p);
	    if (p->sym == G_RBR) {
		/* empty column spec, OK */
		t->v.b2.r = newempty(EMPTY);
	    } else {
		t->v.b2.r = expr(p);
	    }
	    if (p->sym == P_COL) {
		/* second part of colon-separated range */
		t->v.b2.r = newb2(SUBSL, t->v.b2.r, NULL);
		lex(p);
		if (p->sym == G_RBR) {
		    /* reached end: second part implicitly empty */
		    t->v.b2.r->v.b2.r = newempty(EMPTY);
		} else {
		    t->v.b2.r->v.b2.r = expr(p);
		}
	    }
	    if (p->sym == G_RBR) {
		if (p->ch == '\'') {
		    set_transpose(t); /* ?? */
		    parser_getc(p);
		}
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
	expected_symbol_error(cexp, p);
    }

    set_matrix_slice_off();
}

static void attach_child (NODE *parent, NODE *child, int k, int i,
			  parser *p)
{
#if SDEBUG
    fprintf(stderr, "attach_child: i=%d, type = %d\n", i, child->t);
#endif

    if (p->err) {
	return;
    }

    if (k == 2) {
	/* 2-place node */
	if (i == 0) {
	    parent->v.b2.l = child;
	} else {
	    parent->v.b2.r = child;
	}
    } else if (k == 3) {
	/* 3-place node */
	if (i == 0) {
	    parent->v.b3.l = child;
	} else if (i == 1) {
	    parent->v.b3.m = child;
	} else {
	    parent->v.b3.r = child;
	}
    } else {
	/* n-place node */
	p->err = push_bn_node(parent, child);
    }
}

static void pad_parent (NODE *parent, int k, int i, parser *p)
{
    NODE *n;
    int j;

    for (j=i; j<k; j++) {
	n = newempty(EMPTY);
	attach_child(parent, n, k, j, p);
    }
}

/* Get up to k comma-separated arguments (possibly optional).
   However, if k < 0 this is a signal for get as many arguments as we
   can find (the number is unknown in advance).
*/

static void get_args (NODE *t, parser *p, int k, int opt, int *next)
{
    NODE *child;
    char cexp = 0;
    int i = 0;

#if SDEBUG
    fprintf(stderr, "get_args: k = %d...\n", k);
#endif    

    if (p->sym != G_LPR) {
	expected_symbol_error('(', p);
	return;
    }	

    lex(p);

    while (((k > 0 && i < k) || p->ch) && !p->err) {
	if (p->sym == G_RPR) {
	    break;
	}
	child = expr(p);
	attach_child(t, child, k, i++, p);
	if (p->err) {
	    break;
	}
	if (p->sym == G_RPR) {
	    break;
	} else if (p->sym == P_COM) {
	    /* turn off flag for accepting string as first arg */
	    p->flags &= ~P_GETSTR;
	    /* and handle final string arg if relevant */
	    if (i == k - 1 && opt == RIGHT_STR) {
		child = get_string_arg(p);
		attach_child(t, child, k, i++, p);
	    } else {
		lex(p);
	    }
	} else {
	    cexp = ',';
	}
    }

    if (!p->err) {
	if (cexp) {
	    expected_symbol_error(cexp, p);
	} else if (p->sym == G_RPR) {
	    if (i < k) {
		pad_parent(t, k, i, p);
	    } 
	    /* handle trailing transpose */
	    if (p->ch == '\'') {
		set_transpose(t);
		parser_getc(p);
	    }
	    lex(p);
	    if (p->sym == G_LBR) {
		*next = G_LBR;
	    }
	}
    }	
}

static void get_ovar_ref (NODE *t, parser *p)
{
    if (p->ch != '.' || parser_charpos(p, '$') != 0) {
	p->err = E_PARSE;
	return;
    }

    p->idnum = 0;

    /* handle the '.' */
    lex(p);

    /* get the following '$' name */
    lex(p);

    if (p->idnum == 0) {
	p->err = E_PARSE;
    } else if (p->sym == DMSL || p->sym == DMSTR) {
	/* followed by '[' or '(' matrix subspec? */
	t->v.b2.r = powterm(p);
    } else {
	t->v.b2.r = newref(p, p->sym);
	lex(p);
    }
}

static NODE *powterm (parser *p)
{ 
    /* watch out for unary operators */
    int sym = p->sym == B_SUB ? U_NEG : 
	p->sym == B_ADD ? U_POS : p->sym;
    int opt = OPT_NONE;
    int next = 0;
    NODE *t;

    if (p->err) {
	return NULL;
    }

#if SDEBUG
    fprintf(stderr, "powterm: p->sym = %d, p->ch = '%c' (%d)\n",
	    p->sym, p->ch? p->ch : '0', p->ch);
#endif

    if (string_last_func(sym)) {
	opt = RIGHT_STR;
    } 

    if (unary_op(sym)) {
	if (p->ch == 0) {
	    context_error(0, p);
	    return NULL;
	}
        t = newb1(sym, NULL);
        if (t != NULL) {
            lex(p);
            t->v.b1.b = powterm(p);
        }
    } else if (func2_symb(sym)) {
	t = newb2(sym, NULL, NULL);
	if (t != NULL) {
	    lex(p);
	    get_args(t, p, 2, opt, &next);
	}
    } else if (func3_symb(sym)) {
	t = newb3(sym, NULL);
	if (t != NULL) {
	    lex(p);
	    if (string0_func(sym)) {
		p->flags |= P_GETSTR;
	    }
	    get_args(t, p, 3, opt, &next);
	}
    } else if (string0_func(sym)) {
	t = newb1(sym, NULL);
	if (t != NULL) {
	    lex(p);
	    t->v.b1.b = newbn(FARGS);
	    if (t != NULL) {
		p->flags |= P_GETSTR;
		get_args(t->v.b1.b, p, -1, opt, &next);
	    }
	}
    } else if (sym == F_URCPVAL) {
	t = newb1(sym, NULL);
	if (t != NULL) {
	    lex(p);
	    t->v.b1.b = newbn(FARGS);
	    if (t != NULL) {
		get_args(t->v.b1.b, p, -1, opt, &next);
	    }
	}	
    } else if (string_arg_func(sym)) {
	t = newb1(sym, NULL);
	if (t != NULL) {
	    lex(p);
	    t->v.b1.b = get_string_arg(p);
	}	
    } else if (func1_symb(sym)) {
	t = newb1(sym, NULL);
	if (t != NULL) {
	    lex(p);
	    t->v.b1.b = base(p, t);
	    if (p->sym == G_LBR) {
		next = G_LBR;
	    }
	}
    } else if (sym == LAG || sym == OBS) {
	t = newb2(sym, NULL, NULL);
	if (t != NULL) {
	    t->v.b2.l = newref(p, USERIES); 
	    lex(p);
	    t->v.b2.r = base(p, t);
	}
    } else if (sym == MSL || sym == DMSL) {
	t = newb2(sym, NULL, NULL);
	if (t != NULL) {
	    if (p->sym == MSL) {
		t->v.b2.l = newstr(p, STR);
	    } else {
		t->v.b2.l = newref(p, MVAR);
	    }
	    t->v.b2.r = newb2(MSL2, NULL, NULL);
	    if (t->v.b2.r != NULL) {
		lex(p);
		get_slice_parts(t->v.b2.r, p);
	    }
	}
    } else if (sym == DMSTR) {
	t = newb2(sym, NULL, NULL);
	if (t != NULL) {
	    t->v.b2.l = newref(p, MVAR);
	    lex(p);
	    t->v.b2.r = get_string_arg(p);
	}
    } else if (sym == OVAR) {
	t = newb2(sym, NULL, NULL);
	if (t != NULL) {
	    t->v.b2.l = newstr(p, STR);
	    get_ovar_ref(t, p);
	}
    } else if (sym == LISTVAR) {
	t = listvar_node(p);
	if (t != NULL && (p->sym == G_LPR || p->sym == G_LBR)) {
	    /* list.series node may be "inflected" as lag or obs */
	    p->sym = (p->sym == G_LPR)? LAG : OBS;
	    t = newb2(p->sym, t, NULL);
	    if (t != NULL) {
		parser_ungetc(p);
		lex(p);
		t->v.b2.r = base(p, t);
	    }
	}
    } else if (sym == G_LPR) {
	/* dummy root for parenthesized expressions, to facilitate
	   taking the transpose of matrix stuff, e.g. (A*B)' */
	t = newb1(EROOT, NULL);
	if (t != NULL) {
	    t->v.b1.b = base(p, t);
	}
    } else if (sym == G_LCB) {
	/* explicit matrix definition, possibly followed by
	   a "subslice" specification */
	int sub = 0;

	t = newbn(MDEF);
	if (t != NULL) {
	    get_matrix_def(t, p, &sub);
	    if (sub) {
		t = newb2(MSL, t, NULL);
		if (t != NULL) {
		    t->v.b2.r = newb2(MSL2, NULL, NULL);
		    if (t->v.b2.r != NULL) {
			lex(p);
			get_slice_parts(t->v.b2.r, p);
		    }
		}		
	    }
	}
    } else if (sym == UFUN || sym == RFUN) {
	t = newb2(sym, NULL, NULL);
	if (t != NULL) {
	    t->v.b2.l = newstr(p, STR);
	    lex(p);
	    t->v.b2.r = newbn(FARGS);
	    if (t != NULL) {
		get_args(t->v.b2.r, p, -1, opt, &next);
	    }
	}
    } else if (funcn_symb(sym)) {
	t = newb1(sym, NULL);
	if (t != NULL) {
	    lex(p);
	    t->v.b1.b = newbn(FARGS);
	    if (t != NULL) {
		get_args(t->v.b1.b, p, -1, opt, &next);
	    }
	}
    } else if (sym == STR || sym == VSTR) {
	t = newstr(p, sym);
	lex(p);
    } else {
	t = base(p, NULL);
    }

    if (next == G_LBR) {
	/* support func(args)[slice] */
	t = newb2(MSL, t, NULL);
	if (t != NULL) {
	    t->v.b2.r = newb2(MSL2, NULL, NULL);
	    if (t->v.b2.r != NULL) {
		get_slice_parts(t->v.b2.r, p);
	    }
	}	
    }

#if SDEBUG
    notify("powterm", t, p);
#endif

    return t;
}

/* convert B_POW to right associativity: that is,
   
       pow           pow
      L   R         L   R
      |   |   ->    |   |
     pow  c         a  pow
     | |               | |
     a b               b c

*/

static void convert_pow_term (NODE *n)
{
    NODE *L = n->v.b2.l;
    NODE *a = L->v.b2.l;
    NODE *b = L->v.b2.r;
    NODE *c = n->v.b2.r;

    n->v.b2.l = a;
    n->v.b2.r = L;
    L->v.b2.l = b;
    L->v.b2.r = c;
}

#define pow_sym(t)     (t == B_POW || t == B_DOTPOW)
#define pow_pow_node(n) (pow_sym(n->t) && pow_sym(n->v.b2.l->t))

static NODE *factor (parser *p)
{  
    int sym = p->sym == B_SUB ? U_NEG : 
	p->sym == B_ADD ? U_POS : p->sym;
    NODE *t;

    if (p->err) {
	return NULL;
    }

    if (unary_op(sym)) {
	if (p->ch == 0) {
	    context_error(0, p);
	    return NULL;
	}
        t = newb1(sym, NULL);
        if (t != NULL) {
            lex(p);
            t->v.b1.b = factor(p);
        }
    } else {
	t = powterm(p);
	if (t != NULL) {
	    if (p->sym == B_TRMUL && p->ch == 0) {
		/* can't really be TRMUL at end of input */
		set_transpose(t);
	    } else {
		while (!p->err && (p->sym == B_POW || 
				   p->sym == B_DOTPOW ||
				   p->sym == B_TRMUL)) {
		    t = newb2(p->sym, t, NULL);
		    if (t != NULL) {
			lex(p);
			t->v.b2.r = powterm(p);
		    }
		    if (!p->err && pow_pow_node(t)) {
			/* make exponentiation associate rightward */
			convert_pow_term(t);
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

    if (p->err || (t = factor(p)) == NULL) {
	return NULL;
    }

    while (!p->err && (p->sym == B_MUL || p->sym == B_DIV || 
		       p->sym == B_MOD || p->sym == B_DOTMULT || 
		       p->sym == B_DOTDIV || p->sym == B_KRON)) {
	t = newb2(p->sym, t, NULL);
	if (t != NULL) {
	    lex(p);
	    t->v.b2.r = factor(p);
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

    if (p->err || (t = term(p)) == NULL) {
	return NULL;
    }

    while (!p->err && (p->sym == B_ADD || p->sym == B_SUB || 
		       p->sym == B_DOTADD || p->sym == B_DOTSUB ||
		       p->sym == B_HCAT || p->sym == B_VCAT ||
		       p->sym == B_LCAT || p->sym == B_RANGE)) {
	t = newb2(p->sym, t, NULL);
	if (t != NULL) {
	    lex(p);
	    t->v.b2.r = term(p);
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

    if (p->err || (t = expr4(p)) == NULL) {
	return NULL;
    }

    while (!p->err && (p->sym == B_GT || p->sym == B_LT || 
		       p->sym == B_DOTGT || p->sym == B_DOTLT || 
		       p->sym == B_GTE || p->sym == B_LTE)) {
	t = newb2(p->sym, t, NULL);
	if (t != NULL) {
	    lex(p);
	    t->v.b2.r = expr4(p);
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

    if (p->err || (t = expr3(p)) == NULL) {
	return NULL;
    }

    while (!p->err && (p->sym == B_EQ || p->sym == B_NEQ ||
		       p->sym == B_DOTEQ)) {
	t = newb2(p->sym, t, NULL);
	if (t != NULL) {
	    lex(p);
	    t->v.b2.r = expr3(p);
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

    if (p->err || (t = expr2(p)) == NULL) {
	return NULL;
    }

    while (!p->err && p->sym == B_AND) {
	t = newb2(p->sym, t, NULL);
	if (t != NULL) {
	    lex(p);
	    t->v.b2.r = expr2(p);
	}
    }

#if SDEBUG
    notify("expr1", t, p);
#endif

    return t;
}

NODE *expr0 (parser *p)
{  
    NODE *t;

    if (p->err || (t = expr1(p)) == NULL) {
	return NULL;
    }

    while (!p->err && p->sym == B_OR) {
	t = newb2(p->sym, t, NULL);
	if (t != NULL) {
	    lex(p);
	    t->v.b2.r = expr1(p);
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

    if (p->err || (t = expr0(p)) == NULL) {
	return NULL;
    }

    while (!p->err && p->sym == QUERY) {
	t = newb3(p->sym, t);
	if (t != NULL) {
	    lex(p);
	    t->v.b3.m = expr(p);
	    if (p->sym == P_COL) {
		lex(p);
		t->v.b3.r = expr(p);
	    } else {
		expected_symbol_error(':', p);
	    }
	}
    }

#if SDEBUG
    notify("expr", t, p);
#endif

    return t;
}

/* for use when we need to evaluate a sub-matrix specification
   on the left-hand side of a genr formula */

NODE *msl_node_direct (parser *p)
{
    NODE *t;

    t = newb2(MSL2, NULL, NULL);
    if (t != NULL) {
	lex(p);
	get_slice_parts(t, p);
    }

    return t;
}
