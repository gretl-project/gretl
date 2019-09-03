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

/* Debugging routines for matrix_subspec: may be used in geneval.c
   at the initial building on the spec or in usermat.c when
   processing the spec with a matrix in hand.
*/

struct spec_typemap {
    SelType t;
    const char *name;
};

static struct spec_typemap tmap[] = {
    { SEL_NULL,   "SEL_NULL" },
    { SEL_RANGE,  "SEL_RANGE" },
    { SEL_MATRIX, "SEL_MATRIX" },
    { SEL_ALL,    "SEL_ALL" },
    { SEL_DIAG,   "SEL_DIAG" },
    { SEL_UPPER,  "SEL_UPPER" },
    { SEL_LOWER,  "SEL_LOWER" },
    { SEL_REAL,   "SEL_REAL" },
    { SEL_IMAG,   "SEL_IMAG" },
    { SEL_ELEMENT, "SEL_ELEMENT" },
    { SEL_CONTIG, "SEL_CONTIG" },
    { SEL_EXCL,   "SEL_EXCL" },
    { SEL_SINGLE, "SEL_SINGLE" },
    { SEL_STR,    "SEL_STR" }
};

static const char *seltype_name (SelType t)
{
    int n = G_N_ELEMENTS(tmap);
    int i;

    for (i=0; i<n; i++) {
	if (tmap[i].t == t) {
	    return tmap[i].name;
	}
    }

    return "???";
}

#if defined(IN_GENEVAL)

static void print_mspec (matrix_subspec *mspec)
{
    fprintf(stderr, "mspec at %p:\n", (void *) mspec);

    if (mspec != NULL) {
	fprintf(stderr, "  ltype = %s\n", seltype_name(mspec->ltype));
        if (mspec->ltype == SEL_ELEMENT) {
	    fprintf(stderr, "    element = %d\n", mspec->lsel.range[0]);
	} else if (mspec->ltype == SEL_RANGE) {
	    fprintf(stderr, "    lsel.range[0] = %d\n", mspec->lsel.range[0]);
	    fprintf(stderr, "    lsel.range[1] = %d\n", mspec->lsel.range[1]);
	} else if (mspec->ltype == SEL_MATRIX) {
	    gretl_matrix_print(mspec->lsel.m, "sel matrix");
	}
	fprintf(stderr, "  rtype = %s\n", seltype_name(mspec->rtype));
	if (mspec->rtype == SEL_ELEMENT) {
	    fprintf(stderr, "    element = %d\n", mspec->rsel.range[0]);
	} else if (mspec->rtype == SEL_RANGE) {
	    fprintf(stderr, "    rsel.range[0] = %d\n", mspec->rsel.range[0]);
	    fprintf(stderr, "    rsel.range[1] = %d\n", mspec->rsel.range[1]);
	} else if (mspec->rtype == SEL_MATRIX) {
	    gretl_matrix_print(mspec->rsel.m, "sel matrix");
	}
	fputs("end of mspec\n", stderr);
    }
}

#elif defined(IN_USERMAT)

static void subspec_debug_print (matrix_subspec *spec,
				 const gretl_matrix *m)
{
    fprintf(stderr, "matrix_subspec: types = (%s, %s), m is %d x %d\n",
	    seltype_name(spec->ltype), seltype_name(spec->rtype),
	    m->rows, m->cols);
    if (spec->ltype == SEL_MATRIX) {
	fputs(" vector sel,", stderr);
    } else if (spec->ltype != SEL_NULL) {
	fprintf(stderr, " lsel->range = (%d,%d),", spec->lsel.range[0],
		spec->lsel.range[1]);
    }
    if (spec->rtype == SEL_MATRIX) {
	fputs(" vector sel,", stderr);
    } else if (spec->rtype != SEL_NULL) {
	fprintf(stderr, " rsel->range = (%d,%d),", spec->rsel.range[0],
		spec->rsel.range[1]);
    }
    fprintf(stderr, " lh scalar %d, rh scalar %d\n",
	    lhs_is_scalar(spec, m), rhs_is_scalar(spec, m));
}

#endif /* IN_GENEVAL or IN_USERMAT */
