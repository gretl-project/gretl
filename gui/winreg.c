#include <windows.h>
#include <stdio.h>


int set_reg_val (char *keyname, char *keyval)
{
    int error = 0;
    HKEY regkey;

    if (RegCreateKeyEx(
		       HKEY_CLASSES_ROOT,           /* handle to open key */
		       "gretl",                     /* subkey name */
		       0,                           /* reserved */
		       NULL,                        /* class string */
		       REG_OPTION_NON_VOLATILE,     /* special options */
		       KEY_ALL_ACCESS,              /* desired security access */
		       NULL,                        /* inheritance */
		       &regkey,                     /* key handle */
		       NULL                         /* disposition value buffer */
		       ) != ERROR_SUCCESS) {
	fprintf(stderr, "couldn't open registry\n");
	return 1;
    }

    if (RegSetValueEx(
		  regkey,
		  keyname,
		  0,
		  REG_SZ,
		  keyval,
		  strlen(keyval) + 1) != ERROR_SUCCESS) {
	fprintf(stderr, "couldn't set registry entry\n");
	error = 1;
    }
		  
    RegCloseKey(regkey);

    return error;
}

int read_reg_val (char *keyname)
{
    char getval[255];
    unsigned long datalen = 255;
    int winerr, error = 0;
    HKEY regkey;
    LPVOID lpMsgBuf;

    if (RegOpenKeyEx(
		     HKEY_CLASSES_ROOT,           /* handle to open key */
		     "gretl",                     /* subkey name */
		     0,                           /* reserved */
		     KEY_READ,                    /* access mask */
		     &regkey                      /* key handle */
		     ) != ERROR_SUCCESS) {
	fprintf(stderr, "couldn't open registry\n");
	return 1;
    }

    winerr = RegQueryValueEx(
			     regkey,
			     keyname,
			     NULL,
			     NULL,
			     getval,
			     &datalen
			     );
    if (winerr != ERROR_SUCCESS) {
	fprintf(stderr, "couldn't read registry entry (error %d)\n", winerr);
	error = 1;
	FormatMessage( 
		      FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		      FORMAT_MESSAGE_FROM_SYSTEM | 
		      FORMAT_MESSAGE_IGNORE_INSERTS,
		      NULL,
		      GetLastError(),
		      MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), 
		      (LPTSTR) &lpMsgBuf,
		      0,
		      NULL 
		      );
	MessageBox(NULL, (LPCTSTR)lpMsgBuf, "Error",  MB_OK | MB_ICONINFORMATION);
	LocalFree(lpMsgBuf);
    } else {
	fprintf(stderr, "value for '%s' is '%s'\n", keyname, getval);
    }
		  
    RegCloseKey(regkey);

    return error;
}

int main (void)
{
    char *gptpath = "gptpath";
    char *gptpath_val = "c:\\foo\\bar\\gpt.exe";

    set_reg_val(gptpath, gptpath_val);
    read_reg_val(gptpath);

    return 0;
}
