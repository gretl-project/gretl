/*++
            Serves as an intermediate stub Win32 console application to
            avoid a hanging pipe when redirecting 16-bit console based
            programs (including MS-DOS console based programs and batch
            files) on Window 95, Windows 98, and Windows Me.
            This program is to be launched with redirected standard
            handles. It will launch the command line specified 16-bit
            console based application in the same console, forwarding
            it's own redirected standard handles to the 16-bit child.

--*/    

#include <windows.h>

int main (int argc, char *argv[])
{
   BOOL                bRet = FALSE;
   STARTUPINFO         si   = {0};
   PROCESS_INFORMATION pi   = {0};
   char cmdline[256];
   char *workdir;
   int i;

   if (argc < 3) return 1;

   *cmdline = 0;
   for (i=1; i<argc-1; i++) {
       strcat(cmdline, argv[i]);
       if (i < argc-2) strcat(cmdline, " ");
   }

   workdir = argv[argc - 1];

   si.cb = sizeof si;
   si.dwFlags = STARTF_USESHOWWINDOW | STARTF_USESTDHANDLES;
   si.wShowWindow = SW_HIDE;
   si.dwFlags    = STARTF_USESTDHANDLES;
   si.hStdInput  = GetStdHandle (STD_INPUT_HANDLE);
   si.hStdOutput = GetStdHandle (STD_OUTPUT_HANDLE);
   si.hStdError  = GetStdHandle (STD_ERROR_HANDLE);

   bRet = CreateProcess (NULL, cmdline,
                         NULL, NULL,
                         TRUE, 0,
                         NULL, workdir,
                         &si, &pi
                         );
   if (bRet) {
      WaitForSingleObject (pi.hProcess, INFINITE);
      CloseHandle (pi.hProcess);
      CloseHandle (pi.hThread);
      return 0;
   }

   return 1;
}

