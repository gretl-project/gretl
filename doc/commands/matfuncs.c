/* Grab info from gretl function reference, and format a
   table of matrix-related functions arranged by topic.

   Allin Cottrell, April 2008.
*/

#include "libgretl.h"

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#include <libxslt/xslt.h>
#include <libxslt/xsltInternals.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

#include <glib.h>

#define ROOTNODE "funcref"
#define UTF const xmlChar *

#define VERBOSE 0

#define SECTLEN 64

/* the order in which the topics should appear */
enum {
    MATBUILD,
    MATSHAPE,
    MATH,
    LINALG,
    STATS,
    DATA_UTILS,
    FILTERS,
    NUMERICAL,
    PROBDIST,
    STRINGS,
    TRANSFORMS,
    ACCESS,
    TAB_MAX
};

struct tab_labeler {
    int ID;
    const char *tag;
    const char *title;
};

struct tab_labeler labelers[] = {
    { MATBUILD,     "matbuild",     N_("Creation and I/O") },
    { MATSHAPE,     "matshape",     N_("Shape/size/arrangement") },
    { MATH,         "math",         N_("Element by element") },
    { LINALG,       "linalg",       N_("Matrix algebra") },
    { STATS,        "stats",        N_("Statistics/transformations") },
    { DATA_UTILS,   "data-utils",   N_("Data utilities") },
    { FILTERS,      "filters",      N_("Filters") },
    { NUMERICAL,    "numerical",    N_("Numerical methods") },
    { PROBDIST,     "probdist",     N_("Probability distributions") },
    { STRINGS,      "strings",      N_("Strings") },
    { TRANSFORMS,   "transforms",   N_("Transformations") },
    { ACCESS,       "access",       N_("Accessors") },
    { TAB_MAX,      NULL,           NULL }
};

char reffile[FILENAME_MAX];

typedef struct _sectlist sectlist;
typedef struct _section section;
typedef struct _function function;

struct _sectlist {
    int nsects;
    section **sections;
};

struct _section {
    int ID;
    char name[SECTLEN];
    int nfuns;
    function **funs;
};

struct _function {
    char name[16];
    int matrix_ok;
};

static void missing_attrib (const char *element, const char *attrib)
{
    fprintf(stderr, "Required attribute '%s' missing for element '%s'\n",
	    attrib, element);
}

static int approved_section_title (const char *s)
{
    int i;

    for (i=0; labelers[i].tag != NULL; i++) {
	if (!strcmp(s, labelers[i].tag)) {
	    return 1;
	}
    }

    return 0;
}

static section *section_new (const char *name)
{
    section *sect;
    int i;

    sect = malloc(sizeof *sect);
    if (sect != NULL) {
	*sect->name = 0;
	strncat(sect->name, name, SECTLEN - 1);
	sect->nfuns = 0;
	sect->funs = NULL;

	for (i=0; i<TAB_MAX; i++) {
	    if (!strcmp(sect->name, labelers[i].tag)) {
		sect->ID = labelers[i].ID;
		break;
	    }
	}
    }

    return sect;
}

static function *function_new (const char *name, int matrix_ok)
{
    function *fun = malloc(sizeof *fun);

    if (fun != NULL) {
	*fun->name = 0;
	strncat(fun->name, name, 15);
	fun->matrix_ok = matrix_ok;
    }

    return fun;
}

static int section_already_recorded (sectlist *s, const char *name)
{
    int i;

    for (i=0; i<s->nsects; i++) {
	if (!strcmp(name, s->sections[i]->name)) {
	    return 1;
	}
    }

    return 0;
}

static int function_in_gretl (const char *funword)
{
#if 1
    return 1;
#else
    int i;

    for (i=1; i<NC; i++) {
	if (!strcmp(funword, gretl_command_word(i))) {
	    return 1;
	}
    }

    return 0;
#endif
}

static int 
maybe_add_section (xmlDocPtr doc, xmlNodePtr node, sectlist *slist)
{
    section *sect, **sects = NULL;
    char *tmp;
    int ns, err = 0;

    tmp = (char *) xmlGetProp(node, (UTF) "section");
    if (tmp == NULL) {
	missing_attrib("function", "section");
	return 1;
    } 

    if (!approved_section_title(tmp)) {
	fprintf(stderr, "*** Found unapproved section heading '%s'\n", tmp);
    }

#if VERBOSE
    fprintf(stderr, "looking at section '%s'\n", tmp);
#endif

    if (!section_already_recorded(slist, tmp)) {
#if VERBOSE
	fprintf(stderr, " not recorded: adding new section\n");
#endif

	sect = section_new(tmp);
	if (sect == NULL) {
	    return 1;
	}

	ns = slist->nsects + 1;

	sects = realloc(slist->sections, ns * sizeof *sects);
	if (sects == NULL) {
	    fprintf(stderr, "Out of memory\n");
	    return 1;
	}

	slist->sections = sects;
	slist->sections[ns - 1] = sect;
	slist->nsects = ns;
    } 

    free(tmp);

    return err;
}

static int matrix_return (char *rtype)
{
    if (!strcmp(rtype, "matrix") ||
	!strcmp(rtype, "smatrix") ||
	!strcmp(rtype, "rvec") ||
	!strcmp(rtype, "cvec")) {
	return 1;
    } else {
	return 0;
    }
}

static int ok_matrix_arg (char *atype)
{
    int ret = 0;

    if (atype != NULL) {
	if (!strcmp(atype, "vector") ||
	    !strcmp(atype, "rvec") ||
	    !strcmp(atype, "cvec") ||
	    !strcmp(atype, "matrix") ||
	    !strcmp(atype, "smatrix") ||
	    !strcmp(atype, "symmat") ||
	    !strcmp(atype, "matrixref") ||
	    !strcmp(atype, "anyfloat") ||
	    !strcmp(atype, "series-or-vec") ||
	    !strcmp(atype, "series-or-mat") ||
	    !strcmp(atype, "anyfloat-or-list")) {
	    ret = 1;
	}
    }

    return ret;
}

static int 
place_function (xmlDocPtr doc, xmlNodePtr node, sectlist *s)
{
    xmlNodePtr n1, n2;
    function *fun, **funs = NULL;
    char *sname, *fname, *atype, *retval;
    int matrix_ok = 0;
    int gotargs = 0;
    int nc, n = -1;
    int i, err = 0;

    sname = (char *) xmlGetProp(node, (UTF) "section");
    fname = (char *) xmlGetProp(node, (UTF) "name");
    retval = (char *) xmlGetProp(node, (UTF) "output");

    if (!function_in_gretl(fname)) {
	fprintf(stderr, "*** '%s': obsolete function, skipping\n", fname);
	goto bailout;
    }
    
    if (sname == NULL || fname == NULL || retval == NULL) {
	fprintf(stderr, "Error parsing function\n");
	return 1;
    }

    if (matrix_return(retval)) {
	matrix_ok = 1;
    } else {
	n1 = node->xmlChildrenNode;
	while (n1 != NULL && !err && !gotargs) {
	    if (!xmlStrcmp(n1->name, (UTF) "fnargs")) {
		gotargs = 1;
		n2 = n1->xmlChildrenNode;
		while (n2 != NULL && !err && !matrix_ok) {
		    if (!xmlStrcmp(n2->name, (UTF) "fnarg")) {
			atype = (char *) xmlGetProp(n2, (UTF) "type");
			if (ok_matrix_arg(atype)) {
			    matrix_ok = 1;
			}
			free(atype);
		    }
		    n2 = n2->next;
		}
	    }
	    n1 = n1->next;
	} 
    }  

    if (!matrix_ok) {
	goto bailout;
    }

    for (i=0; i<s->nsects; i++) {
	if (!strcmp(sname, s->sections[i]->name)) {
	    n = i;
	    break;
	}
    }

    if (n < 0) {
	fprintf(stderr, "Couldn't place function '%s'\n", fname);
	return 1;
    }

    fun = function_new(fname, matrix_ok);
    if (fun == NULL) {
	fprintf(stderr, "Out of memory\n");
	return 1;
    }

    nc = s->sections[n]->nfuns + 1;
    funs = realloc(s->sections[n]->funs, nc * sizeof *funs);
    if (funs == NULL) {
	fprintf(stderr, "Out of memory\n");
	return 1;
    }
	
    s->sections[n]->funs = funs;
    s->sections[n]->funs[nc-1] = fun;
    s->sections[n]->nfuns = nc;

#if VERBOSE
    fprintf(stderr, "Added function '%s' to section '%s'\n", fname, sname);
#endif

 bailout:

    free(fname);
    free(sname);

    return err;
}

static int parse_ref_file (sectlist *slist)
{
    xmlDocPtr doc;
    xmlNodePtr cur, flist = NULL;
    char *listname;
    int err = 0;

    LIBXML_TEST_VERSION 
	xmlKeepBlanksDefault(0);

    xmlSubstituteEntitiesDefault(1);
    xmlLoadExtDtdDefaultValue = 1;

    doc = xmlParseFile(reffile);
    if (doc == NULL) {
	err = 1;
	goto bailout;
    }

    cur = xmlDocGetRootElement(doc);
    if (cur == NULL) {
	err = 1;
	goto bailout;
    }

    if (xmlStrcmp(cur->name, (UTF) ROOTNODE)) {
	fprintf(stderr, "File of the wrong type, root node not %s\n", 
		ROOTNODE);
	err = 1;
	goto bailout;
    }

    /* first pass: find the right funclist */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (UTF) "funclist")) {
	    listname = (char *) xmlGetProp(cur, (UTF) "ref");
	    if (!strcmp(listname, "functions")) {
		flist = cur;
		free(listname);
		break;
	    } else {
		free(listname);
	    }
	}
	cur = cur->next;
    } 

    if (flist == NULL) {
	fprintf(stderr, "Couldn't find functions list\n");
	err = 1;
	goto bailout;
    }

    /* first pass: walk the tree, picking up section headings */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (UTF) "function")) {
	    err = maybe_add_section(doc, cur, slist);
	}
	cur = cur->next;
    }   

    /* second pass: assemble functions in sections */
    if (!err) {
	cur = flist;
	cur = cur->xmlChildrenNode;
	while (cur != NULL && !err) {
	    if (!xmlStrcmp(cur->name, (UTF) "function")) {
		err = place_function(doc, cur, slist);
	    }
	    cur = cur->next;
	}
    } 

 bailout:

    xmlFreeDoc(doc);
    xmlCleanupParser();

    return err;
}

static void sectlist_init (sectlist *s)
{
    s->nsects = 0;
    s->sections = NULL;
}

void free_section (section *sect)
{
    int i;

    for (i=0; i<sect->nfuns; i++) {
	free(sect->funs[i]);
    }
    free(sect->funs);
    free(sect);
}

void free_sectlist (sectlist *s)
{
    int i;

    for (i=0; i<s->nsects; i++) {
	free_section(s->sections[i]);
    }
    free(s->sections);
}

static section *get_section_by_id (sectlist *s, int ID)
{
    int i;

    for (i=0; i<s->nsects; i++) {
	if (s->sections[i]->ID == ID) {
	    return s->sections[i];
	}
    }

    return NULL;
}

/* print out the matrix-ok functions as a TeX table, organized
   by "theme" */

const char *tabline =
    "\\begin{tabular}{p{\\cwid}p{\\cwid}p{\\cwid}p{\\cwid}p{\\cwid}p{\\cwid}}";

static int print_topic_lists (sectlist *s)
{
    const section *sect;
    int i, j, k;
    int err = 0;

    for (i=0; i<TAB_MAX; i++) {
	sect = get_section_by_id(s, i);
	if (sect == NULL || sect->nfuns == 0) {
	    fprintf(stderr, "Section '%s': no entries\n", labelers[i].title);
	    continue;
	}
	printf("\\textbf{%s}\n\\hrulefill\n\n", _(labelers[i].title));
	puts(tabline);
	for (j=0; j<sect->nfuns; j++) {
	    printf("\\texttt{%s}", sect->funs[j]->name);
	    if (j == sect->nfuns - 1) {
		for (k=j+1; k%6; k++) {
		    fputs( " &", stdout);
		}
		puts(" \\\\ [4pt]");
	    } else if ((j+1) % 6 == 0) {
		puts(" \\\\");
	    } else {
		puts(" &");
	    }
	}
	puts("\\end{tabular}\n");
    }

    return err;
}

/* so we can localize the section titles */

void nls_init (void)
{
    setlocale(LC_ALL, "");
    bindtextdomain(PACKAGE, LOCALEDIR);
    textdomain(PACKAGE);
    /* note: for TeX output */
    bind_textdomain_codeset(PACKAGE, "ISO-8859-1");
}

int main (int argc, char **argv)
{
    sectlist slist;
    int err = 0;

    if (argc != 2) {
	fprintf(stderr, "Please supply one argument: the name of a "
		"file to process\n");
	exit(EXIT_FAILURE);
    }

    nls_init();

    strcpy(reffile, argv[1]);

    sectlist_init(&slist);

    err = parse_ref_file(&slist);
    if (err) {
	fprintf(stderr, "Error parsing %s\n", reffile);
	exit(EXIT_FAILURE);
    }

    fprintf(stderr, "Found %d sections in '%s'\n", slist.nsects, reffile);

    print_topic_lists(&slist);

    free_sectlist(&slist);

    return err;
}
