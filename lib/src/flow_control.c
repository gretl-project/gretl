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

/* if-then stuff - conditional execution */

#include "libgretl.h"
#include "cmd_private.h"
#include "flow_control.h"

#define IFDEBUG 0

enum {
    SET_FALSE,
    SET_TRUE,
    SET_ELSE,
    SET_ELIF,
    SET_ENDIF,
    IS_FALSE,
    IS_TRUE,
    UNINDENT,
    GETINDENT,
    RELAX,
    IFRESET
};

enum {
    TOK_IF = 1,
    TOK_ELIF,
    TOK_ELSE,
    TOK_ENDIF
};

#if 0 /* not yet */

static int inline_if (const char *s)
{
    int ret = 0;

    if (strchr(s, ';') != NULL) {
	int quoted = 0;

	while (*s) {
	    if (*s == '"') {
		quoted = !quoted;
	    } else if (!quoted && *s == ';') {
		ret = 1;
		break;
	    }
	    s++;
	}
    }

    return ret;
}

#endif

/* if_eval: evaluate an "if" condition by generating a scalar
   (integer) representing the truth or falsity of the condition.
   The condition is expressed in the string @s. If a loop is
   being executed currently the @ptr argument may be non-NULL:
   this will happen if the conditional is known NOT to be
   subject to string substitution, in which case it can be
   "compiled" and reused.
*/

static int if_eval (int ci, const char *s, DATASET *dset,
		    PRN *prn, void *ptr, int *err)
{
    GENERATOR *ifgen = NULL;
    double val = NADBL;
    int ret = -1;

#if IFDEBUG
    if (s != NULL) {
	fprintf(stderr, "if_eval: s = '%s'\n", s);
    }
#endif

    if (ptr != NULL) {
	/* We're being called from a loop, with the implicit
	   request that the if-condition be "compiled" (if
	   that's not already done) and subsequently executed
	   without having to be evaluated from scratch.
	*/
	ifgen = *(GENERATOR **) ptr;
	if (ifgen == NULL && s != NULL) {
	    /* Generator not compiled yet: do it now. The
	       flag OPT_P indicates that we're generating
	       a "private" scalar.
	    */
	    GENERATOR **pgen = (GENERATOR **) ptr;

	    *pgen = ifgen = genr_compile(s, dset, GRETL_TYPE_BOOL,
					 OPT_P | OPT_N, NULL, err);
#if IFDEBUG
	    fprintf(stderr, "if_eval: genr_compile (%s), err = %d\n",
		    s, *err);
#endif
	}
    }

    if (ifgen != NULL) {
	val = evaluate_if_cond(ifgen, dset, prn, err);
#if IFDEBUG
	fprintf(stderr, "ran @ifgen: val %g, err %d\n", val, *err);
#endif
    } else if (s == NULL) {
	*err = E_DATA;
    } else {
	*err = 0;
	val = generate_boolean(s, dset, prn, err);
    }

#if IFDEBUG > 1
    fprintf(stderr, "if_eval: generate returned %d\n", *err);
#endif

    if (*err) {
	gretl_errmsg_append(_("error evaluating 'if'"), 0);
    } else if (na(val)) {
	*err = 1;
	gretl_errmsg_append(_("indeterminate condition for 'if'"), 0);
    } else {
	ret = (int) val;
    }

    if (*err && s != NULL && *s != '\0') {
	gchar *cond = g_strdup_printf("> %s %s", gretl_command_word(ci), s);

	gretl_errmsg_append(cond, 0);
	g_free(cond);
    }

#if IFDEBUG
    fprintf(stderr, "if_eval: returning %d\n", ret);
#endif

    return ret;
}

#if IFDEBUG > 1
static const char *ifstr (int c)
{
    if (c == SET_FALSE) return "SET_FALSE";
    if (c == SET_TRUE)  return "SET_TRUE";
    if (c == SET_ELSE)  return "SET_ELSE";
    if (c == SET_ELIF)  return "SET_ELIF";
    if (c == SET_ENDIF) return "SET_ENDIF";
    if (c == IS_FALSE)  return "IS_FALSE";
    if (c == IS_TRUE)   return "IS_TRUE";
    if (c == UNINDENT)  return "UNINDENT";
    if (c == GETINDENT) return "GETINDENT";
    if (c == RELAX)     return "RELAX";
    if (c == IFRESET)   return "RESET";
    return "UNKNOWN";
}
#endif

static void unmatched_message (int code)
{
    gretl_errmsg_sprintf(_("Unmatched \"%s\""),
			 (code == SET_ELSE)? "else" :
			 (code == SET_ELIF)? "elif" : "endif");
}

#define IF_DEPTH 1024

/* Note: the @got_T boolean array below is used to record,
   within an "if ... endif" block, whether any true condition
   has been encountered. This is relevant in the case of
   blocks containing "elif" clauses (as well as, possibly,
   an "else" clause): by reference to @got_T we can ensure
   that at most one branch is followed. Otherwise we'd be in
   danger of following all branches for which the "if" part of
   an "elif" condition turns out to be true, so disregarding
   the "else" implicit in "elif".

   A simple example:

   x = 2
   if x == 1
      print "x = 1"
   elif x == 2
      print "x = 2"
   elif x < 3
      print "x < 3"
   endif

   Use of @got_T in this case ensures that the branch with
   condition x < 3 is not followed.
*/

static int ifstate (int code, int val, int *err)
{
    static unsigned char T[IF_DEPTH];
    static unsigned char tok[IF_DEPTH];
    static unsigned char got_T[IF_DEPTH];
    static unsigned short indent;
    int i, ret = 0;

#if IFDEBUG > 1
    if (code != IS_FALSE) {
	fprintf(stderr, "ifstate: code = %s\n", ifstr(code));
    }
#endif

    if (code == IS_FALSE || code == IS_TRUE) {
	for (i=1; i<=indent; i++) {
	    if (T[i] == 0) {
		ret = 1; /* blocked */
		break;
	    }
	}
	return code == IS_TRUE ? !ret : ret;
    } else if (code == RELAX) {
	indent = 0;
    } else if (code == IFRESET) {
	indent = val;
    } else if (code == GETINDENT) {
	ret = indent;
    } else if (code == UNINDENT) {
	ret = --indent;
    } else if (code == SET_FALSE || code == SET_TRUE) {
	indent++;
	if (indent >= IF_DEPTH) {
	    gretl_errmsg_sprintf("IF depth (%d) exceeded", IF_DEPTH);
	    *err = E_DATA;
	} else {
	    T[indent] = got_T[indent] = (code == SET_TRUE);
	    tok[indent] = TOK_IF;
	}
    } else if (code == SET_ELIF || code == SET_ELSE) {
	if (tok[indent] != TOK_IF && tok[indent] != TOK_ELIF) {
	    unmatched_message(code);
	    *err = E_PARSE;
	} else {
	    tok[indent] = (code == SET_ELSE)? TOK_ELSE : TOK_ELIF;
	    if (T[indent]) {
		T[indent] = 0;
	    } else if (!got_T[indent]) {
		T[indent] = 1;
	    }
	}
    } else if (code == SET_ENDIF) {
	if (tok[indent] != TOK_IF &&
	    tok[indent] != TOK_ELIF &&
	    tok[indent] != TOK_ELSE) {
	    unmatched_message(code);
	    *err = E_PARSE;
	} else {
	    tok[indent] = TOK_ENDIF;
	    got_T[indent] = 0;
	    indent--;
	}
    }

#if IFDEBUG > 1
    fprintf(stderr, "ifstate: returning %d (indent %d, err %d)\n",
	    ret, indent, (err == NULL)? 0 : *err);
#endif

    return ret;
}

static int set_if_state (int code)
{
    int err = 0;

    ifstate(code, 0, &err);
    return err;
}

static int get_if_state (int code)
{
    int err = 0;

    return ifstate(code, 0, &err);
}

void gretl_if_state_clear (void)
{
#if IFDEBUG
    fprintf(stderr, "gretl_if_state_clear called\n");
#endif
    ifstate(RELAX, 0, NULL);
}

int gretl_if_state_finalize (void)
{
    int ret, err = 0;

    ret = ifstate(IS_TRUE, 0, NULL);

    if (!ret) {
	ifstate(RELAX, 0, NULL);
	err = E_PARSE;
    }

    return err;
}

int gretl_if_state_record (void)
{
    return ifstate(GETINDENT, 0, NULL);
}

void gretl_if_state_reset (int indent)
{
    ifstate(IFRESET, indent, NULL);
}

int gretl_if_state_false (void)
{
    return get_if_state(IS_FALSE);
}

int gretl_if_state_check (int indent0)
{
    int indent = ifstate(GETINDENT, 0, NULL);
    int err = 0;

    if (indent != indent0) {
	gretl_errmsg_sprintf(_("Unmatched \"%s\""), "if");
	ifstate(RELAX, 0, NULL);
	err = E_PARSE;
    }

    return err;
}

/* flow_control: if the ci (command index) member of @cmd
   is something other than one of the flow control symbols
   IF, ELSE, ELIF or ENDIF, this function simply returns
   1 if execution if blocked by a false IF condition or 0
   otherwise.

   If we get one of the flow control symbols we operate
   on the program's "if state", pushing a term onto, or
   popping a term off, the existing stack. And in this case
   we always return 1, which indicates to the machinery in
   interact.c that execution of the current command is
   completed.

   We need the @dset argument in case we have to evaluate a
   new IF condition.
*/

int flow_control (ExecState *s, DATASET *dset, void *ptr)
{
    CMD *cmd = s->cmd;
    int ci = cmd->ci;
    int blocked = get_if_state(IS_FALSE);
    int ok, err = 0;

    if (ci != IF && ci != ELSE && ci != ELIF && ci != ENDIF) {
	return blocked;
    }

    if (ci == IF) {
	if (blocked) {
	    /* just increase the "indent" level */
	    err = set_if_state(SET_FALSE);
	} else {
	    /* actually evaluate the condition */
	    ok = if_eval(ci, cmd->vstart, dset, s->prn, ptr, &err);
	    if (!err) {
		err = set_if_state(ok? SET_TRUE : SET_FALSE);
	    }
	}
    } else if (ci == ENDIF) {
	err = set_if_state(SET_ENDIF);
    } else if (ci == ELIF) {
	err = set_if_state(SET_ELIF);
	if (!err && get_if_state(IS_TRUE)) {
	    set_if_state(UNINDENT);
	    ok = if_eval(ci, cmd->vstart, dset, s->prn, ptr, &err);
	    if (!err) {
		err = set_if_state(ok? SET_TRUE : SET_FALSE);
	    }
	}
    } else if (ci == ELSE) {
	err = set_if_state(SET_ELSE);
    }

    if (err) {
	set_if_state(RELAX);
	cmd->err = err;
    }

    return 1;
}
