/* longname plugin for gretl -- allow win98 and higher to
   access the win32 function GetLongPathName() */

int real_unmangle (const char *dosname, char *longname, int maxlen)
{
    int err = GetLongPathName(dosname, longname, maxlen);

    if (err > 0 && err <= maxlen) err = 0;
    else err = 1;
    return err;
}


