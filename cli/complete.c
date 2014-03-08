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

#ifdef WIN32
# include "winconfig.h"
#else
# include "config.h"
#endif

#include <stdio.h> 
#include <stdlib.h>
#include <string.h>

#include "gretl_commands.h"

#include <readline/readline.h>
#include <readline/history.h> 

char *dupstr (const char *s) 
{
    char *r = malloc(strlen(s) + 1);

    if (r != NULL) {
	strcpy(r, s);
    }

    return r;
}

/* Generator function for command completion.  STATE lets us know whether
   to start from scratch; without any state (i.e. STATE == 0), then we
   start at the top of the list. */

static char *command_generator (char *text, int state)
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

/* Attempt to complete on the contents of TEXT.  START and END bound the
   region of rl_line_buffer that contains the word to complete.  TEXT is
   the word to complete.  We can use the entire contents of rl_line_buffer
   in case we want to do some simple parsing.  Return the array of matches,
   or NULL if there aren't any. */

static char **gretl_completion (char *text, int start, int end)
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

/* Read a string, and return a pointer to it.  Returns NULL on EOF. */

char *rl_gets (char **line_read, const char *prompt)
{
    /* If the buffer has already been allocated, return the memory
       to the free pool. */
    if (*line_read != NULL) {
	free(*line_read);
	*line_read = NULL;
    }
     
    /* Get a line from the user. */
    *line_read = readline(prompt);

    /* If the line has any text in it, save it on the history. */
    if (*line_read != NULL && **line_read != '\0') {
	add_history(*line_read);
    }
     
    return *line_read;
}

/* Tell the GNU Readline library how to complete.  We want to try to complete
   on command names if this is the first word in the line, or on variable
   names if not. */

void initialize_readline (void) 
{
    /* Allow conditional parsing of the ~/.inputrc file. */
    rl_readline_name = "gretl";

    /* Tell the completer that we want a crack first. */
    rl_attempted_completion_function = (rl_completion_func_t *) gretl_completion;
}

