@ECHO OFF

SET "STACK_PATH=..\stack"

clang++.exe .\gretlapi.cpp -o .\gretlapi-clang-cpp.exe^
 -I%STACK_PATH%\include -I%STACK_PATH%\include\glib-2.0 -I%STACK_PATH%\lib\glib-2.0\include^
 -L%STACK_PATH%\lib^
 -l"gretl-1.0" -lxml2 -lfftw3 -l"glib-2.0" -lintl -l"gobject-2.0"

REM using a driver program for clang that attempts to be compatible with MSVC's cl.exe.
clang-cl.exe .\gretlapi.cpp -o .\gretlapi-clang-cl-cpp.exe^
 /I%STACK_PATH%\include /I%STACK_PATH%\include\glib-2.0 /I%STACK_PATH%\lib\glib-2.0\include^
 /link /LIBPATH:"%STACK_PATH%\lib"^
 gretl-1.0.lib xml2.lib fftw3.lib glib-2.0.lib intl.lib gobject-2.0.lib
