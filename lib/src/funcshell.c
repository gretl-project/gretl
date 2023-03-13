struct funcshell_ {
    char name[FN_NAMELEN]; /* identifier */
    UfunAttrs flags;       /* private, etc. */
    char *code;            /* actual code */
    int n_params;          /* number of parameters */
    fn_param *params;      /* parameter info array */
    int rettype;           /* return type (if any) */
};

typedef struct funcshell_ funcshell;

struct pkgshell_ {
    char *fname;           /* filename */
    funcshell **pub;       /* pointers to public interfaces */
    funcshell **priv;      /* pointers to private functions */
    int n_pub;             /* number of public functions */
    int n_priv;            /* number of private functions */
};

typedef struct pkgshell_ pkgshell;

static pkgshell *pkgshell_new (const char *fname)
{
    pkgshell *pks = calloc(1, sizeof *pks);

    if (pks != NULL) {
	pks->fname = gretl_strdup(fname);
    }

    return pks;
}

static int read_funcshell_from_xml (xmlNodePtr node,
				    xmlDocPtr doc,
				    pkgshell *pks)
{
    funcshell *fns = calloc(1, sizeof *fns);
    xmlNodePtr cur;
    char *tmp;
    int err = 0;

    if (fns == NULL) {
	return E_ALLOC;
    }

    if (gretl_xml_get_prop_as_string(node, "name", &tmp)) {
	strncat(fns->name, tmp, FN_NAMELEN - 1);
	free(tmp);
    } else {
	return E_DATA;
    }

    if (gretl_xml_get_prop_as_string(node, "type", &tmp)) {
	fns->rettype = return_type_from_string(tmp, &err);
	free(tmp);
    } else {
	fun->rettype = GRETL_TYPE_VOID;
    }

    if (gretl_xml_get_prop_as_bool(node, "private")) {
	fns->flags |= UFUN_PRIVATE;
    }

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "params")) {
	    err = funcshell_read_params(cur, doc, fns);
	} else if (!xmlStrcmp(cur->name, (XUC) "code")) {
	    err = funcshell_read_code(cur, doc, fns);
	}
	cur = cur->next;
    }

    if (err) {
	fns_free(fns);
    } else {
	; /* add to @pks */
    }

    return err;
}

static pkgshell *real_read_pkgshell (xmlDocPtr doc,
				     xmlNodePtr node,
				     const char *fname,
				     int *err)
{
    xmlNodePtr cur;
    pkgshell *pks;
    char *tmp = NULL;

    pks = pkgshell_new(fname);
    if (pks == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    cur = node->xmlChildrenNode;

    // cur = node->xmlChildrenNode;
    while (cur != NULL && !*err) {
	if (!xmlStrcmp(cur->name, (XUC) "gretl-function")) {
	    *err = read_funcshell_from_xml(cur, doc, pks);
	}
	cur = cur->next;
    }

    return pks;
}

static pkgshell *read_pkg_as_shell (const char *fname, int *err)
{
    pkgshell *pks = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr node = NULL;
    xmlNodePtr cur;

    *err = gretl_xml_open_doc_root(fname, "gretl-functions", &doc, &node);
    if (*err) {
	return NULL;
    }

    cur = node->xmlChildrenNode;
    
    while (cur != NULL && !*err) {
	if (!xmlStrcmp(cur->name, (XUC) "gretl-function-package")) {
	    pks = real_read_pkgshell(doc, cur, fname, err);
	    break;
	}
	cur = cur->next;
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    if (!*err && pks == NULL) {
	*err = E_DATA;
    }

    return pks;
}

int print_all_gfn_code (const char *fname, PRN *prn)
{
    pkgshell *pks = NULL;
    int err = 0;

    pks = read_pkg_as_shell(fname, &err);

    if (!err) {
	print_pkgshell_code(pks, prn);
    }

    free_pkgshell(pks);

    return err;
}
