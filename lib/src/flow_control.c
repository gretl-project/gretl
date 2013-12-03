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

/* if_eval: evaluate an "if" condition by generating a scalar
   (integer) representing the truth or falsity of the condition.
   The condition is expressed in the string @s. If a loop is
   being executed currently the @ptr argument may be non-NULL:
   this will happen if the conditional is known NOT to be 
   subject to string substitution, in which case it can be
   "compiled" and reused.
*/

static int if_eval (const char *s, DATASET *dset, void *ptr, int *err)
{
    GENERATOR *ifgen = NULL;
    double val = NADBL;
    int ret = -1;

#if IFDEBUG
    fprintf(stderr, "if_eval: s = '%s'\n", s);
#endif

    while (*s == ' ') s++;

    if (ptr != NULL) {
	/* We're being called from a loop, with the implicit
	   request that the if-condition be "compiled" (if
	   that's not already done) and subsequently executed
	   without having to be evaluated from scratch.
	*/
	ifgen = *(GENERATOR **) ptr;

	if (ifgen == NULL) {
	    /* Generator not compiled yet: do it now. The
	       flags OPT_P and OPT_S indicate that we're
	       generating a "private" scalar.
	    */
	    GENERATOR **pgen = (GENERATOR **) ptr;

	    *pgen = ifgen = genr_compile(s, dset, OPT_P | OPT_S, err);
	}
    }

    if (ifgen != NULL) {
	val = evaluate_if_cond(ifgen, dset, err);
    } else {
	*err = 0;
	val = generate_scalar(s, dset, err);
    }

#if IFDEBUG
    fprintf(stderr, "if_eval: generate returned %d\n", *err);
#endif

    if (*err) {
	gretl_errmsg_set(_("error evaluating 'if'"));
    } else if (na(val)) {
	*err = 1;
	gretl_errmsg_set(_("indeterminate condition for 'if'"));
    } else {
	ret = (int) val;
    }

#if IFDEBUG
    fprintf(stderr, "if_eval: returning %d\n", ret);
#endif

    return ret;
}

#if IFDEBUG
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
			 (code == SET_ELIF)? "elif": "endif");
}

#define IF_DEPTH 1024

static int ifstate (int code, int val, int *err)
{
    static unsigned short T[IF_DEPTH];
    static unsigned short got_if[IF_DEPTH];
    static unsigned short got_else[IF_DEPTH];
    static unsigned short got_T[IF_DEPTH];
    static unsigned short indent;
    int i, ret = 0;

#if IFDEBUG
    fprintf(stderr, "ifstate: code = %s\n", ifstr(code));
#endif

    if (code == RELAX) {
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
	    got_if[indent] = 1;
	    got_else[indent] = 0;
	}
    } else if (code == SET_ELSE || code == SET_ELIF) {
	if (got_else[indent] || !got_if[indent]) {
	    unmatched_message(code);
	    *err = E_PARSE;
	} else {
	    got_else[indent] = (code == SET_ELSE);
	    if (T[indent]) {
		T[indent] = 0;
	    } else if (!got_T[indent]) {
		T[indent] = 1;
	    }
	}
    } else if (code == SET_ENDIF) {
	if (!got_if[indent] || indent == 0) {
	    unmatched_message(code);
	    *err = E_PARSE;
	} else {
	    got_if[indent] = 0;
	    got_else[indent] = 0;
	    got_T[indent] = 0;
	    indent--;
	}
    } else if (code == IS_FALSE || code == IS_TRUE) {
	for (i=1; i<=indent; i++) {
	    if (T[i] == 0) {
		ret = 1;
		break;
	    }
	}
	if (code == IS_TRUE) {
	    ret = !ret;
	}
    } 

#if IFDEBUG
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

static int trailing_junk_error (const char *s)
{
    char junk[16] = {0};
 
    s += strspn(s, " \t");
    sscanf(s, "%15[^ ]", junk);
    gretl_errmsg_sprintf(_("field '%s' in command is invalid"), junk);
    return E_PARSE;
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

   We need the @line (command line) and @dset arguments
   in case we have to evaluate a new IF condition.
*/

int flow_control (const char *line, DATASET *dset, CMD *cmd,
		  void *ptr)
{
    int ci = cmd->ci;
    int blocked, ok, err = 0;

    blocked = get_if_state(IS_FALSE);

    if (ci != IF && ci != ELSE && ci != ELIF && ci != ENDIF) {
	return blocked;
    }

    if (ci == IF) {
	if (blocked) {
	    err = set_if_state(SET_FALSE);
	} else {
	    ok = if_eval(line + 2, dset, ptr, &err);
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
	    ok = if_eval(line + 4, dset, ptr, &err);
	    if (!err) {
		err = set_if_state(ok? SET_TRUE : SET_FALSE);
	    }
	}
    } else if (ci == ELSE) {
	if (!string_is_blank(line + 4)) {
	    err = trailing_junk_error(line + 4);
	} else {
	    err = set_if_state(SET_ELSE);
	}
    }

    if (err) {
	set_if_state(RELAX);
	cmd->err = err;
    } 

    return 1;
}

