/* longname plugin for gretl -- allow win98 and higher to
   access the win32 function GetLongPathName() */

#undef _WIN32_WINDOWS
#define _WIN32_WINDOWS 0x0410
#include <windows.h>

__declspec(dllexport) 
void real_unmangle (const char *dosname, char *longname, 
		    int maxlen, int *err)
{
    if (dosname == NULL || longname == NULL) {
	*err = 1;
	return;
    }

    *err = GetLongPathName(dosname, longname, maxlen);

    if (*err > 0 && *err <= maxlen) *err = 0;
    else *err = 1;
}



