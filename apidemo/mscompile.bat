@ECHO off
REM This program is intended exclusively for native MS Windows compilation.
REM You need the libgretl stack for MS SDK, which can be downloaded from:
REM https://sourceforge.net/projects/gretl/files/mscompile/
REM or
REM https://sourceforge.net/projects/libgretl-for-ms-sdk/

REM architecture we're building for: x86_64 or aarch64
SET "ARCH=x86_64"

REM gretl installation path
SET "PREFIX=\"C:\\\\"Program Files\\gretl"\""

REM libgretl stack path
SET "STACK_PATH=..\..\libgretl-for-ms-sdk-code\trunk\stack"

REM #############################################################
REM standard clang.exe
REM #############################################################
SET "CC=clang.exe"
SET "CFLAGS=-march=native -mtune=native -O2 -DPREFIX=%PREFIX%"
SET "GRETLIN=-I%STACK_PATH%\include -I%STACK_PATH%\include\glib-2.0 -I%STACK_PATH%\lib\%ARCH%\glib-2.0\include"
SET "GRETLLD=-L%STACK_PATH%\lib\%ARCH% -l"gretl-1.0" -lxml2 -lfftw3 -l"glib-2.0" -lintl -l"gobject-2.0""

REM simple_client
%CC% %CFLAGS% %GRETLIN% %GRETLLD% simple_client.c -o simple_client.exe

REM arma_example
%CC% %CFLAGS% %GRETLIN% %GRETLLD% arma_example.c -o arma_example.exe

REM nls_example
%CC% %CFLAGS% %GRETLIN% %GRETLLD% nls_example.c -o nls_example.exe

REM #############################################################
REM clang-cl.exe (MSVC syntax)
REM #############################################################
SET "CC=clang-cl.exe"
SET "CFLAGS=/O2 /D PREFIX=%PREFIX%"
SET "GRETLIN=/I%STACK_PATH%\include /I%STACK_PATH%\include\glib-2.0 /I%STACK_PATH%\lib\%ARCH%\glib-2.0\include"
SET "GRETLLD=/link /LIBPATH:%STACK_PATH%\lib\%ARCH% gretl-1.0.lib xml2.lib fftw3.lib glib-2.0.lib intl.lib gobject-2.0.lib"

REM simple_client
%CC% simple_client.c /Fe:simple_client_cl.exe %CFLAGS% %GRETLIN% %GRETLLD%

REM arma_example
%CC% arma_example.c /Fe:arma_example_cl.exe %CFLAGS% %GRETLIN% %GRETLLD%

REM nls_example
%CC% nls_example.c /Fe:nls_example_cl.exe %CFLAGS% %GRETLIN% %GRETLLD%

ECHO.
ECHO Remember to set up paths to libgretl binaries, for example:
ECHO SET PATH=%%PATH%%;C:\Program Files\gretl
ECHO SET "GRETL_PLUGIN_PATH=C:\Program Files\gretl\plugins"
ECHO SET "GRETL_DATA_PATH=C:\Program Files\gretl\data"
ECHO.

ECHO ON