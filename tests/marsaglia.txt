The *.c files in this directory were created by the
program f2c, which translates Fortran code to C code.
The primary purpose of this set of files is to provide
means to create  six executable files:
 diehard, makewhat, diequick, getrnor, asc2bin, meld, asc2bin
from the six C files
 diehard.c, makewhat.c, diequick.c, getrnor.c, meld.c, asc2bin.c
thus getting the full DIEHARD battery of tests.

Please note that the Fortran source files for DIEHARD
are a mixture of patched and addended Fortran files
going back over thirty years.  They were gathered to
provide the DIEHARD battery of tests and the makewhat
exec file for generating various kinds of RNG's.  The
resulting mix seems to suit my purposes, but if you
want to examine particular *.c files for details of
testing or generating procedures, you may have difficulty.
The f2c translation is meant for machines, not people,
and even if relation to the original Fortran file were
transparent, you might still have problems because the
Fortran file may itself be a mix of earlier working
versions, or because of steps necessary to overcome
Fortran's quaint insistence on signed integers and the
un-mathematical way it decides I < J? for 32-bit integers.

Whatever, if you compile these *.c files you should get
a working version of DIEHARD.   Be sure to link with
the f2c and math libraries.  A typical command might go
(for the gnu c compiler):

    gcc -o diehard diehard.c -lf2c -lm

which would produce the executable file that does the
DIEHARD battery of tests.   Similar commands would
create the other exec files---`makewhat' for creating
test files, `asc2bin' for changing an ascii (hex) file
of integers to the binary form required by `diehard',
`meld' for merging two binary files of random bits,
`getrnor' for creating a file of standard normal variables.
And `diequick' for a terse version of `diehard'.

These exec files need to read from three text or data files
     tests.txt   makef.txt  operm5d.ata
that must be available to them.

There is also a diehard.doc file in this package, as well
as three postscript files that give theory and descriptions
of tests and random number generators.   There is a fourth
postcript file, cdmake.ps, that describes the CDROM I will
be distributing soon.

George Marsaglia
geo@stat.fsu.edu
 
