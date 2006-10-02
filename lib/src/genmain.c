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

/* driver module for 'genr' and related commands */

#include "genparse.h"
#include "libset.h"
#include "gretl_func.h"

#if GENDEBUG
# define GDEBUG 1
#else
# define GDEBUG 0
#endif

static void gen_write_message (const parser *p, int oldv, PRN *prn)
{
    if (prn == NULL || !gretl_messages_on()) {
	return;
    }

    if (p->targ == NUM) {
	if (var_is_series(p->dinfo, p->lh.v)) {
	    pprintf(prn, "Modified series %s (ID %d)",
		    p->lh.name, p->lh.v);
	} else if (p->flags & P_DECL) {
	    pprintf(prn, "Added scalar %s (ID %d)",
		    p->lh.name, p->lh.v);
	} else {
	    double x = (*p->Z)[p->lh.v][p->lh.obs];

	    if (p->lh.v < oldv) {
		pprintf(prn, "Replaced scalar %s (ID %d)",
			p->lh.name, p->lh.v);
	    } else {
		pprintf(prn, "Generated scalar %s (ID %d)",
			p->lh.name, p->lh.v);
	    }
	    if (na(x)) {
		pputs(prn, " = NA");
	    } else {
		pprintf(prn, " = %g", x);
	    }
	}
    } else if (p->targ == VEC) {
	if (p->flags & P_DECL) {
	    pprintf(prn, "Added series %s (ID %d)",
		    p->lh.name, p->lh.v);
	} else if (p->lh.v < oldv) {
	    pprintf(prn, "Replaced series %s (ID %d)",
		    p->lh.name, p->lh.v);
	} else if (p->lh.v < oldv) {
	    pprintf(prn, "Replaced series %s (ID %d)",
		    p->lh.name, p->lh.v);
	} else {
	    pprintf(prn, "Generated series %s (ID %d)",
		    p->lh.name, p->lh.v);
	}
    } else if (p->targ == MAT) {
	if (p->flags & P_DECL) {
	    pprintf(prn, "Added matrix %s\n", p->lh.name);
	} else if (p->lh.m0 != NULL && p->lh.substr != NULL && 
		   *p->lh.substr !='\0') {
	    pprintf(prn, "Modified matrix %s\n", p->lh.name);
	} else if (p->lh.m0 != NULL) {
	    pprintf(prn, "Replaced matrix %s\n", p->lh.name);
	} else {
	    pprintf(prn, "Generated matrix %s\n", p->lh.name);
	}
    }

    pputc(prn, '\n');
}

static void gen_write_label (parser *p, int oldv)
{
    char tmp[MAXLABEL];
    const char *src;
    size_t len = 0;

    /* don't touch the label if we generated a single
       observation in a series */
    if (p->lh.substr != NULL) {
	return;
    }

    *tmp = '\0';

    if (p->lh.v < oldv) {
	int m = get_model_count();

	if (m > 0) {
	    sprintf(tmp, _("Replaced after model %d: "), m);
	    len = strlen(tmp);
	}
    }

    if (*p->lh.label != '\0' && dollar_node(p->tree)) {
	src = p->lh.label;
    } else {
	src = p->rhs;
    }

    if (strlen(src) > MAXLABEL - 1 - len) {
	strncat(tmp, src, MAXLABEL - 4 - len);
	strcat(tmp, "...");
    } else {
	strncat(tmp, src, MAXLABEL - 1);
    }

    strcpy(VARLABEL(p->dinfo, p->lh.v), tmp);
    p->dinfo->varinfo[p->lh.v]->flags |= VAR_GENERATED;

#if GDEBUG
    fprintf(stderr, "varlabel: '%s'\n", VARLABEL(p->dinfo, p->lh.v));
#endif
}

/**
 * function_from_string:
 * @s: the string to look up.
 *
 * Returns: 1 if there is a funtion corresponding
 * to the name @s, or 0 if there is no such function.
 */

int function_from_string (const char *s)
{
    char word[VNAMELEN];
    const char *p;

    *word = 0;

    p = strchr(s, '(');
    if (p != NULL && p - s <= 8) {
	strncat(word, s, p - s);
    } else {
	strncat(word, s, 8);
    }

    if (function_lookup(word)) {
	return 1;
    }

    /* aliases */
    if (!strcmp(word, "ln")) {
	return 1;
    }

    /* user-defined functions */
    if (gretl_is_user_function(s)) {
	return 1;
    }

    return 0;
}

/* "reserved word" data */

static const char *res1[] = {
    "const",
    "CONST",
    "pi",
    "NA",
    "null"
};

static const char *res2[] = {
    "t"
};

static const char *res3[] = {
    "obs",
    "scalar",
    "series",
    "matrix"
};

/**
 * gretl_reserved_word:
 * @str: string to be tested.
 *
 * Returns non-zero if @str is a reserved word that cannot figure as the
 * name of a user-defined variable, otherwise 0.
 */

int gretl_reserved_word (const char *str)
{
    const char *ruses[] = {
	N_("constant"),
	N_("plotting variable"),
	N_("internal variable"),
	N_("math function")
    };
    static int n1 = sizeof res1 / sizeof res1[0];
    static int n2 = sizeof res2 / sizeof res2[0];
    static int n3 = sizeof res3 / sizeof res3[0];
    int i, ret = 0;

    for (i=0; i<n1 && !ret; i++) {
	if (!strcmp(str, res1[i])) {
	    ret = 1;
	}
    }

    for (i=0; i<n2 && !ret; i++) {
	if (!strcmp(str, res2[i])) {
	    ret = 2;
	}
    }

    for (i=0; i<n3 && !ret; i++) {
	if (!strcmp(str, res3[i])) {
	    ret = 3;
	}
    }

    if (!ret) {
	if (function_from_string(str)) {
	    ret = 4;
	}
    }

    if (ret > 0) {
	sprintf(gretl_errmsg, _("'%s' refers to a %s and may not be used as a "
			    "variable name"), str, _(ruses[ret - 1])); 
    }
 
    return ret;
}

#define GEN_LEVEL_DEBUG 0

/**
 * varindex:
 * @pdinfo: data information struct.
 * @varname: name of variable to test.
 *
 * Returns: the ID number of the variable whose name is given,
 * or the next available ID number if there is no variable of
 * that name.
 */

int varindex (const DATAINFO *pdinfo, const char *varname)
{
    const char *s = varname;
    int fsd = 0;
    int i, ret = pdinfo->v;

    if (s == NULL || *s == 0) {
	return ret;
    }

    if (isdigit(*s)) {
	return ret;
    }

    if (!strcmp(s, "const")) {
	return 0;
    }

    fsd = gretl_function_depth();

#if GEN_LEVEL_DEBUG
    fprintf(stderr, "varindex for '%s': fsd = %d\n", s, fsd);
#endif

    if (fsd > 0) {
	/* inside a function: see only vars at that level */
	for (i=1; i<pdinfo->v; i++) { 
	    if (STACK_LEVEL(pdinfo, i) == fsd && 
		!strcmp(pdinfo->varname[i], s)) {
		ret = i;
		break;
	    }
	}
    } else {
	/* see all vars */
	for (i=1; i<pdinfo->v; i++) { 
	    if (!strcmp(pdinfo->varname[i], s)) { 
		ret = i;
		break;
	    }
	}
    }

#if GEN_LEVEL_DEBUG
    fprintf(stderr, "varindex: '%s': returning %d (pdinfo->v = %d)\n", 
	    s, ret, pdinfo->v);
#endif 

    return ret;
}

static int gen_special (const char *s, double ***pZ, 
			DATAINFO *pdinfo, PRN *prn,
			parser *p)
{
    int orig_v = pdinfo->v;
    int msg = 0;
    int err = 0;

    if (!strcmp(s, "dummy")) {
	int di0 = dummy(pZ, pdinfo, 0);

	if (di0 == 0) {
	    err = 1;
	} else if (di0 == orig_v) {
	    pputs(prn, _("Periodic dummy variables generated.\n"));
	} else {
	    pputs(prn, _("Periodic dummy variables already present.\n"));
	}
    } else if (!strcmp(s, "timedum")) {
	err = panel_dummies(pZ, pdinfo, OPT_T);
	if (!err) {
	    pputs(prn, _("Panel dummy variables generated.\n"));
	}
    } else if (!strcmp(s, "unitdum")) {
	err = panel_dummies(pZ, pdinfo, OPT_NONE);
	if (!err) {
	    pputs(prn, _("Panel dummy variables generated.\n"));
	}
    } else if (!strcmp(s, "time")) {
	err = gen_time(pZ, pdinfo, 1);
	msg = 1;
    } else if (!strcmp(s, "index")) {
	err = gen_time(pZ, pdinfo, 0);
	msg = 1;
    } else if (!strcmp(s, "unit")) {
	err = gen_unit(pZ, pdinfo);
	msg = 1;
    } else if (!strcmp(s, "weekday")) {
	err = gen_wkday(pZ, pdinfo);
	msg = 1;
    }	

    if (!err && msg) {
	strcpy(p->lh.name, s);
	p->lh.v = varindex(pdinfo, s);
	p->Z = pZ;
	p->dinfo = pdinfo;
	p->targ = VEC;
	p->flags = 0;
	p->err = p->warn = 0;
	p->prn = prn;
	gen_write_message(p, orig_v, prn);
    }	    

    return err;
}

static int gen_special_call (const char *s)
{
    if (strncmp(s, "genr ", 5)) {
	return 0;
    }

    s += 5;

    if (!strcmp(s, "dummy") || 
	!strcmp(s, "timedum") || 
	!strcmp(s, "unitdum") || 
	!strcmp(s, "time") || 
	!strcmp(s, "index") || 
	!strcmp(s, "unit") || 
	!strcmp(s, "weekday")) {
	return 1;
    }

    return 0;
}

#define gen_verbose(f) (!(f & P_PRINT) && \
                        !(f & P_DISCARD) && \
                        !(f & P_PRIVATE) && \
                        !(f & P_UFUN))

int generate (const char *line, double ***pZ, DATAINFO *pdinfo,
	      gretlopt opt, PRN *prn)
{
    int oldv = pdinfo->v;
    int flags = 0;
    parser p;

    if (opt & OPT_P) {
	flags |= P_PRIVATE;
    }

    if (opt & OPT_U) {
	flags |= P_UFUN;
    }

#if GDEBUG
    fprintf(stderr, "\n*** generate: line = '%s'\n", line);
#endif

    if (gen_special_call(line)) {
	return gen_special(line + 5, pZ, pdinfo, prn, &p);
    }

    realgen(line, &p, pZ, pdinfo, prn, flags);

    if (!(opt & OPT_U)) {
	gen_save_or_print(&p, prn);
    }

    if (!p.err && gen_verbose(p.flags)) {
	gen_write_label(&p, oldv);
	if (!(opt & OPT_Q)) {
	    gen_write_message(&p, oldv, prn);
	}
    }

    gen_cleanup(&p);

#if GDEBUG
    fprintf(stderr, "generate: returning %d\n", p.err);
#endif

    return p.err;
}

/* simply retrieve a scalar result */

double generate_scalar (const char *s, double ***pZ, 
			DATAINFO *pdinfo, int *err)
{
    parser p;
    double x = NADBL;

    *err = realgen(s, &p, pZ, pdinfo, NULL, P_SCALAR | P_PRIVATE);

    if (!*err) {
	if (p.ret->t == MAT) {
	    x = p.ret->v.m->val[0];
	} else {
	    x = p.ret->v.xval;
	}
    }

    gen_cleanup(&p);

    return x;
}

/* retrieve and print a variable from "within" a saved
   object
*/

int print_object_var (const char *oname, const char *param,
		      double ***pZ, DATAINFO *pdinfo,
		      PRN *prn)
{
    char line[MAXLEN];
    parser p;

    sprintf(line, "%s.%s", oname, param);
    realgen(line, &p, pZ, pdinfo, prn, P_DISCARD);
    gen_save_or_print(&p, prn);
    gen_cleanup(&p);

    return p.err;
}

/* create a parsed tree that can be evaluated later, 
   probably multiple times */

parser *genr_compile (const char *s, double ***pZ, DATAINFO *pdinfo, 
		      int *err)
{
    parser *p = malloc(sizeof *p);

#if GDEBUG
    fprintf(stderr, "\n*** genr_compile: s = '%s'\n", s);
#endif

    if (p == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    *err = realgen(s, p, pZ, pdinfo, NULL, P_COMPILE | P_PRIVATE);

    return p;
}

/* run a previously compiled generator */

int execute_genr (parser *p, double ***pZ, DATAINFO *pdinfo)
{
    PRN *prn = NULL;

#if GDEBUG
    fprintf(stderr, "\n*** execute_genr: p = %p\n s='%s'\n", 
	    (void *) p, p->input);
    gretl_print_new(GRETL_PRINT_STDERR);
#endif

    realgen(NULL, p, pZ, pdinfo, prn, P_EXEC);
    gen_save_or_print(p, prn);
    gen_cleanup(p);

#if GDEBUG
    gretl_print_destroy(prn);
#endif

    return p->err;
}

/* destroy a previously compiled generator */

void destroy_genr (parser *p)
{
#if GDEBUG
    fprintf(stderr, "\n*** destroy_genr: p = %p\n", (void *) p);
#endif

    p->flags = 0;
    gen_cleanup(p);
    free(p);
}

int genr_get_varnum (const parser *p)
{
    return p->lh.v;
}
