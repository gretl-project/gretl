/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef WIN32
# include "config.h"
#endif

#include <stdio.h> 
#include <stdlib.h>
#include <string.h>

#include "gretl_commands.h"

#include <readline/readline.h>
#include <readline/history.h> 

extern char *xmalloc();

char *dupstr (const char *s) 
{
    char *r = xmalloc(strlen(s) + 1);

    strcpy(r, s);
    return r;
}

/* ........................................................... */

static char *command_generator (char *text, int state)
/* Generator function for command completion.  STATE lets us know whether
   to start from scratch; without any state (i.e. STATE == 0), then we
   start at the top of the list. */
{
    static int list_index;
    const char *cword;

    /* If this is a new word to complete, initialize now. */
    if (!state) {
	list_index = 0;
    }

    /* Return the next name which partially matches from the command list. */
    cword = gretl_command_complete_next(text, list_index++);
    if (cword != NULL) {
	return dupstr(cword);
    }
    
    return NULL;
}

/* ........................................................... */

static char **gretl_completion (char *text, int start, int end)
/* Attempt to complete on the contents of TEXT.  START and END bound the
   region of rl_line_buffer that contains the word to complete.  TEXT is
   the word to complete.  We can use the entire contents of rl_line_buffer
   in case we want to do some simple parsing.  Return the array of matches,
   or NULL if there aren't any. */
{
    char **matches = NULL;

    /* If this word is at the start of the line, then it is a command
       to complete.  Otherwise it is the name of a variable.
    */
#ifdef NEW_READLINE 
    if (start == 0) {
	matches = 
	    rl_completion_matches(text, 
				  (rl_compentry_func_t *) command_generator);
    }
#else
    if (start == 0) {
	matches = completion_matches(text, command_generator);
    }
#endif

    return matches;
}

/* ........................................................... */

char *rl_gets (char **line_read, const char *prompt)
     /* Read a string, and return a pointer to it.  Returns NULL on EOF. */
{
    /* If the buffer has already been allocated, return the memory
       to the free pool. */
    if (*line_read) {
	free(*line_read);
	*line_read = NULL;
    }
     
    /* Get a line from the user. */
    *line_read = readline(prompt);

    /* If the line has any text in it, save it on the history. */
    if (*line_read && **line_read) {
	add_history(*line_read);
    }
     
    return *line_read;
}

/* ........................................................... */

void initialize_readline (void) 
/* Tell the GNU Readline library how to complete.  We want to try to complete
   on command names if this is the first word in the line, or on variable
   names if not. */
{
    /* Allow conditional parsing of the ~/.inputrc file. */
    rl_readline_name = "gretl";

    /* Tell the completer that we want a crack first. */
    rl_attempted_completion_function = (CPPFunction *) gretl_completion;
}

