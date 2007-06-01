/*
 *   Copyright (c) by Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
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

enum {
    STR_COPY,
    STR_STEAL
};

static NODE *powterm (parser *p);

#if SDEBUG
static void notify (const char *s, NODE *t, parser *p)
{
    fprintf(stderr, "%-8s: returning node at %p, err = %d\n", 
	    s, (void *) t, p->err);
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
	n->ext = 0;
	n->aux = 0;
	n->tmp = 0;
    }

    return n;
}

static NODE *newref (parser *p)
{  
    NODE *n = malloc(sizeof *n);

#if MDEBUG
    fprintf(stderr, "newref: allocated node at %p (type = %d, idnum = %d)\n", 
	    (void *) n, p->sym, p->idnum);
#endif

    if (n != NULL) {
	n->t = p->sym;
	if (n->t == UMAT || n->t == UOBJ ||
	    n->t == LOOPIDX || n->t == LIST) {
	    n->v.str = p->idstr;
	} else {
	    n->v.idnum = p->idnum;
	}
	n->ext = 0;
	n->aux = 0;
	n->tmp = 0;
    }

    return n;
}

static NODE *newstr (char *s, int ext, int flag)
{  
    NODE *n;

    if (s == NULL) {
	fprintf(stderr, "newstr: input is NULL\n");
	return NULL;
    }

    n = malloc(sizeof *n);

#if MDEBUG
    fprintf(stderr, "newstr: allocated node at %p (s = '%s')\n", 
	    (void *) n, s);
#endif

    if (n != NULL) {
	n->t = STR;
	if (flag == STR_COPY) {
	    n->v.str = gretl_strdup(s);
	} else {
	    n->v.str = s;
	}
	n->tmp = 0;
	n->ext = ext;
	n->aux = 0;
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
	n->tmp = 0;
	n->ext = 0;
	n->aux = 0;
    }

    return n;
}

/* node for unary operator, or single-argument function */

static NODE *newb1 (int t, NODE *b, int ext)
{  
    NODE *n = malloc(sizeof *n);

#if MDEBUG
    fprintf(stderr, "newb1:  allocated node at %p (type = %d)\n", 
	    (void *) n, t);
#endif

    if (n != NULL) {
	n->t = t;
	n->v.b1.b = b;
	n->tmp = 0;
	n->ext = ext;
	n->aux = 0;
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
	n->tmp = 0;
	n->ext = 0;
	n->aux = 0;

	if (n->t == AST2) {
	    n->t = B_POW;
	    n->ext = 1; /* record possible ambiguity */
	}
    }

    return n;
}

/* node for ternary operator */

static NODE *newb3 (int t, NODE *l, NODE *m, NODE *r)
{  
    NODE *n = malloc(sizeof *n);

#if MDEBUG
    fprintf(stderr, "newb3:  allocated node at %p (type = %d)\n", 
	    (void *) n, t);
#endif

    if (n != NULL) {
	n->t = t;
	n->v.b3.l = l;
	n->v.b3.m = m;
	n->v.b3.r = r;
	n->tmp = 0;
	n->ext = 0;
	n->aux = 0;
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
	n->tmp = 0;
	n->ext = 0;
	n->aux = 0;
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
	    "added node of type %d\n", t->v.bn.n_nodes,
	    n->t);
#endif

    return 0;
}

static void expected_symbol_error (int c, parser *p)
{
    parser_print_input(p);
    pprintf(p->prn, _("Expected '%c' but found '%s'\n"), c, 
	    getsymb(p->sym, p));
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
	/* peek at the next non-space character */
	int c = parser_next_char(p);

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
    case UVAR: 
    case UMAT:
    case UOBJ:
    case CON: 
    case DVAR:
    case MVAR:
    case OBS:
    case LIST:
    case LOOPIDX:
	t = newref(p);
	if (matrix_ref_node(p) && unary_apost(p)) {
	    t->ext = TRANSP;
	    parser_getc(p);
	}
	lex(p);
	break;
    case DUM:
	if (p->idnum == DUM_NULL) {
	    t = newempty(EMPTY);
	} else {
	    t = newref(p);
	}
	lex(p);
	break;
    case B_SUB:
    case B_ADD:
    case DOT:
	lex(p);
	t = expr(p);
	break;
    case LPR: /* left paren '(' */
	lex(p);
	t = expr(p);
	if (p->sym == RPR) {
	    if (up != NULL && up->t != LAG && unary_apost(p)) {
		up->ext = TRANSP; 
		parser_getc(p);
	    } 
	    lex(p);
	} else if (p->err == 0) {
	    expected_symbol_error(')', p);
	}
	break;
    case LBR: /* left bracket '[' */
	if (up == NULL) {
	    goto deferr;
	}
	if (up->t == OBS) {
	    t = obs_node(p);
	}
	if (p->sym == RBR) {
	    if (up->t == MSL || up->t == DMSL) {
		if (p->ch == '\'') {
		    up->ext = TRANSP;
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

/* grab parenthesized string argument, for some special
   functions */

static NODE *get_string_arg (parser *p)
{
    char str[GENSTRLEN] = {0};
    int i, close = -1;

    if (p->ch != ')') {
	/* allow for empty arg string "()" */
	int j, started = 0;

	close = parser_charpos(p, ')');

	if (close < 0 || close > GENSTRLEN - 2) {
	    p->err = E_PARSE;
	    if (close > 0) {
		pprintf(p->prn, _("String is too long (%d versus %d max)\n"),
			close, GENSTRLEN);
	    } else {
		unmatched_symbol_error('(', p);
	    }
	    return NULL;
	}

	j = 0;
	for (i=0; i<=close; i++) {
	    if (!started && !isspace(p->ch)) {
		started = 1;
	    }
	    if (started) {
		str[j++] = p->ch;
	    }
	    parser_getc(p);
	}
    }

    parser_getc(p);
    lex(p);

    return newstr(tailstrip(str), 0, STR_COPY);
}

enum {
    BOTH_OPT = 1,
    RIGHT_OPT
};

/* gather an unknown number of comma-separated arguments
   for a (possibly user-defined) function */

static void get_multi_args (NODE *t, parser *p)
{
    NODE *n;
    char cexp = 0;

#if SDEBUG
    fprintf(stderr, "get_multi_args, p->sym = %d\n", p->sym);
#endif    

    if (p->sym == LPR) {
	lex(p);
	while (p->ch != 0 && !p->err) {
	    n = expr(p);
	    if (p->err) {
		break;
	    } else {
		p->err = push_bn_node(t, n);
	    }
	    if (p->sym == COM) {
		/* in case acceptance of plain strings was set,
		   turn it off after the first arg */
		p->getstr = 0;
		lex(p);
	    } else if (p->sym == RPR) {
		break;
	    }
	}
    } else {
	cexp = '(';
    }

    if (cexp == 0 && !p->err) {
	if (p->sym == RPR) {
	    lex(p);
	} else {
	    unmatched_symbol_error('(', p);
	}
    }
	    
    if (cexp && p->err == 0) {
	expected_symbol_error(cexp, p);
    }
}

static void get_matrix_def (NODE *t, parser *p, int *sub)
{
    NODE *n;
    char cexp = 0;

#if SDEBUG
    fprintf(stderr, "get_matrix_def, p->sym = %d\n", p->sym);
#endif    

    if (p->sym == LCB) {
	lex(p);
	while (p->ch != 0 && !p->err) {
	    n = expr(p);
	    if (p->err) {
		break;
	    } else {
		p->err = push_bn_node(t, n);
	    }
	    if (p->sym == COM) {
		lex(p);
	    } else if (p->sym == SEMI) {
		n = newempty(EMPTY);
		p->err = push_bn_node(t, n);
		lex(p);
	    } else if (p->sym == RCB) {
		if (p->ch == '\'') {
		    t->ext = TRANSP;
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
	if (p->sym == RCB) {
	    lex(p);
	} else {
	    unmatched_symbol_error('{', p);
	}
    }
	    
    if (cexp && p->err == 0) {
	expected_symbol_error(cexp, p);
    }
}	

static void get_slice_parts (NODE *t, parser *p)
{
    char cexp = 0;

#if SDEBUG
    fprintf(stderr, "get_slice_parts, p->sym = %d\n", p->sym);
#endif    

    if (p->sym == LBR) {
	lex(p);
	if (p->sym == COM) {
	    /* empty row spec, OK */
	    t->v.b2.l = newempty(EMPTY);
	} else {
	    t->v.b2.l = expr(p);
	}
	if (p->sym == COL) {
	    /* need second part of colon-separated range */
	    t->v.b2.l = newb2(SUBSL, t->v.b2.l, NULL);
	    lex(p);
	    t->v.b2.l->v.b2.r = expr(p);
	}
	if (p->sym == RBR) {
	    /* co comma, no second arg string: may be OK */
	    t->v.b2.r = newempty(ABSENT);
	    lex(p);
	    return;
	}
	if (p->sym == COM) {
	    lex(p);
	    if (p->sym == RBR) {
		/* empty column spec, OK */
		t->v.b2.r = newempty(EMPTY);
	    } else {
		t->v.b2.r = expr(p);
	    }
	    if (p->sym == COL) {
		/* need second part of colon-separated range */
		t->v.b2.r = newb2(SUBSL, t->v.b2.r, NULL);
		lex(p);
		t->v.b2.r->v.b2.r = expr(p);
	    }
	    if (p->sym == RBR) {
		if (p->ch == '\'') {
		    t->ext = TRANSP; /* ?? */
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
	    
    if (cexp && p->err == 0) {
	expected_symbol_error(cexp, p);
    }
}

/* get up to two comma-separated arguments 
   (possibly optional) */

static void get_args (NODE *t, parser *p, int opt)
{
    char cexp = 0;

#if SDEBUG
    fprintf(stderr, "get_args...\n");
#endif    

    if (p->sym == LPR) {
	lex(p);
	if (p->sym == RPR && opt == BOTH_OPT) {
	    /* no args, but it's OK */
	    t->v.b2.l = newempty(EMPTY);
	    t->v.b2.r = newempty(EMPTY);
	    lex(p);
	    return;
	}
	t->v.b2.l = expr(p);
	if (p->sym == RPR && opt == RIGHT_OPT) {
	    /* no second arg, but it's OK */
	    t->v.b2.r = newempty(EMPTY);
	    lex(p);
	    return;
	}	    
	if (p->sym == COM) {
	    lex(p);
	    t->v.b2.r = expr(p);
	    if (p->sym == RPR) {
		if (p->ch == '\'') {
		    t->ext = TRANSP; 
		    parser_getc(p);
		}
		lex(p);
	    } else {
		cexp = ')';
	    }
	} else {
	    cexp = ',';
	}
    } else {
	cexp = '(';
    }
	    
    if (cexp && p->err == 0) {
	expected_symbol_error(cexp, p);
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
    } else if (p->sym == DMSL) {
	/* followed by '[' slice mechanism? */
	t->v.b2.r = powterm(p);
    } else {
	t->v.b2.r = newref(p);
	lex(p);
    }
}

#define idnum_to_ext(t) (t == LAG || t == OBS || t == MVAR || \
                         t == DMSL || t == DMSTR)

static NODE *powterm (parser *p)
{  
    int ext = (idnum_to_ext(p->sym)) ? p->idnum : 0;
    int opt = 0;
    NODE *t;

    if (p->err) {
	return NULL;
    }

#if SDEBUG
    fprintf(stderr, "powterm: p->sym = %d, p->ch = %c\n",
	    p->sym, p->ch? p->ch : '0');
#endif

    if (p->sym == UNIFORM || p->sym == NORMAL) {
	opt = BOTH_OPT;
    }

    if (func2_symb(p->sym)) {
	t = newb2(p->sym, NULL, NULL);
	if (t != NULL) {
	    lex(p);
	    get_args(t, p, opt);
	}
    } else if (string0_func(p->sym)) {
	t = newb1(p->sym, NULL, 0);
	if (t != NULL) {
	    lex(p);
	    t->v.b1.b = newbn(FARGS);
	    if (t != NULL) {
		p->getstr = 1;
		get_multi_args(t->v.b1.b, p);
	    }
	}	
    } else if (string_arg_func(p->sym)) {
	t = newb1(p->sym, NULL, ext);
	if (t != NULL) {
	    lex(p);
	    t->v.b1.b = get_string_arg(p);
	}	
    } else if (func_symb(p->sym)) {
	/* includes LAG, OBS */
	t = newb1(p->sym, NULL, ext);
	if (t != NULL) {
	    lex(p);
	    t->v.b1.b = base(p, t);
	}
    } else if (p->sym == MSL || p->sym == DMSL) {
	t = newb2(p->sym, NULL, NULL);
	if (t != NULL) {
	    t->v.b2.l = newstr(p->idstr, ext, STR_STEAL);
	    t->v.b2.r = newb2(MSL2, NULL, NULL);
	    if (t->v.b2.r != NULL) {
		lex(p);
		get_slice_parts(t->v.b2.r, p);
	    }
	}
    } else if (p->sym == DMSTR) {
	t = newb2(p->sym, NULL, NULL);
	if (t != NULL) {
	    t->v.b2.l = newstr(p->idstr, ext, STR_STEAL);
	    lex(p);
	    t->v.b2.r = get_string_arg(p);
	}
    } else if (p->sym == OVAR) {
	t = newb2(p->sym, NULL, NULL);
	if (t != NULL) {
	    t->v.b2.l = newstr(p->idstr, 0, STR_STEAL);
	    get_ovar_ref(t, p);
	}
    } else if (p->sym == LPR) {
	/* dummy root for parenthesized expressions, to facilitate
	   taking the transpose of matrix stuff, e.g. (A*B)' */
	t = newb1(EROOT, NULL, 0);
	if (t != NULL) {
	    t->v.b1.b = base(p, t);
	}
    } else if (p->sym == LCB) {
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
    } else if (p->sym == UFUN) {
	t = newb2(p->sym, NULL, NULL);
	if (t != NULL) {
	    t->v.b2.l = newstr(p->idstr, 0, STR_STEAL);
	    lex(p);
	    t->v.b2.r = newbn(FARGS);
	    if (t != NULL) {
		get_multi_args(t->v.b2.r, p);
	    }
	}
    } else if (funcn_symb(p->sym)) {
	t = newb1(p->sym, NULL, 0);
	if (t != NULL) {
	    lex(p);
	    t->v.b1.b = newbn(FARGS);
	    if (t != NULL) {
		get_multi_args(t->v.b1.b, p);
	    }
	}
    } else if (p->sym == STR) {
	t = newstr(p->idstr, 0, STR_STEAL);
	lex(p);
    } else {
	t = base(p, NULL);
    }

#if SDEBUG
    notify("powterm", t, p);
#endif

    return t;
}

#if 0
/* chunk that may be needed, modified */
if (unary_op(sym)) {
    if (p->ch == 0) {
	/* no input left: provoke error */
	p->sym = sym;
	t = base(p, NULL);
    } else {
	/* continue as usual */
    }
}
#endif


static NODE *factor (parser *p)
{  
    int sym = p->sym == B_SUB ? U_NEG : 
	p->sym == B_ADD ? U_POS : 
	p->sym == B_AND ? U_ADDR :
	p->sym;
    NODE *t;

    if (p->err) {
	return NULL;
    }

    if (unary_op(sym)) {
	if (p->ch == 0) {
	    context_error(0, p);
	    return NULL;
	}
        t = newb1(sym, NULL, 0);
        if (t != NULL) {
            lex(p);
            t->v.b1.b = factor(p);
        }
    } else {
	t = powterm(p);
	if (t != NULL) {
	    while (!p->err && (p->sym == B_POW || 
			       p->sym == AST2 ||
			       p->sym == DOTPOW ||
			       p->sym == B_TRMUL)) {
		t = newb2(p->sym, t, NULL);
		if (t != NULL) {
		    lex(p);
		    t->v.b2.r = powterm(p);
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
		       p->sym == B_MOD || p->sym == DOTMULT || 
		       p->sym == DOTDIV || p->sym == KRON)) {
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
		       p->sym == DOTADD || p->sym == DOTSUB ||
		       p->sym == MCCAT || p->sym == MRCAT)) {
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
		       p->sym == DOTGT || p->sym == DOTLT || 
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
		       p->sym == DOTEQ)) {
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
	t = newb3(p->sym, t, NULL, NULL);
	if (t != NULL) {
	    lex(p);
	    t->v.b3.m = expr0(p);
	    if (p->sym == COL) {
		lex(p);
		t->v.b3.r = expr0(p);
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
