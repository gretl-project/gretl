/* us2a4.c -- a simple filter designed to input a gretl 
   TeX driver file with US letter-size geometry and output
   a version with geometry suitable for A4 paper.

   Allin Cottrell, September 2011.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define A4GEOM "\\usepackage[a4paper,body={6.07in,9.7in}," \
               "top=.8in,left=1.1in]{geometry}"

int main (void)
{
    char buf[2048];
    int geom_done = 0;

    while (fgets(buf, sizeof buf, stdin)) {
	if (!geom_done && !strncmp(buf, "\\usepackage[letter", 18)) {
	    puts(A4GEOM);
	    geom_done = 1;
	} else {
	    fputs(buf, stdout);
	}
    }

    return 0;
}
