/* longname plugin for gretl -- allow win98 and higher to
   access the win32 function GetLongPathName() */

#include <windows.h>

__declspec(dllexport) 
void real_unmangle (const char *dosname, char *longname, 
		    int maxlen, int *err)
{
    *err = GetLongPathName(dosname, longname, maxlen);

    if (*err > 0 && *err <= maxlen) *err = 0;
    else *err = 1;
}



