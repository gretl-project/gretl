/* Lightweight version of user-function struct for use when
   we're just extracting a function from a gfn for printing.
*/

struct funcshell_ {
    char name[FN_NAMELEN]; /* identifier */
    UfunAttrs flags;       /* private, etc. */
    char *code;            /* actual code */
    int n_params;          /* number of parameters */
    fn_param *params;      /* parameter info array */
    int rettype;           /* return type (if any) */
};

typedef struct funcshell_ funcshell;

/* Lightweight version of function-package struct for use when
   we're just printing the code content of a gfn.
*/

struct pkgshell_ {
    funcshell **pub;       /* pointers to public interfaces */
    funcshell **priv;      /* pointers to private functions */
    int n_pub;             /* number of public functions */
    int n_priv;            /* number of private functions */
};

typedef struct pkgshell_ pkgshell;

static void funcshell_destroy (funcshell *fns)
{
    free_params_array(fns->params, fns->n_params);
    free(fns->code);
    free(fns);
}

static void pkgshell_destroy (pkgshell *pks)
{
    int i;

    for (i=0; i<pks->n_priv; i++) {
        funcshell_destroy(pks->priv[i]);
    }
    for (i=0; i<pks->n_pub; i++) {
        funcshell_destroy(pks->pub[i]);
    }
    free(pks->priv);
    free(pks->pub);
    free(pks);
}

/* Attach the function-shell @fns to the package-shell @pfs,
   under its private or public array as appropriate.
*/

static int pkgshell_attach (funcshell *fns,
                            pkgshell *pks)
{
    funcshell **tmp = NULL;
    funcshell ***pfs;
    int n, *pnf;

    pfs = (fns->flags & UFUN_PRIVATE)? &pks->priv : &pks->pub;
    pnf = (fns->flags & UFUN_PRIVATE)? &pks->n_priv : &pks->n_pub;
    n = *pnf + 1;
    tmp = realloc(*pfs, n * sizeof *tmp);

    if (tmp == NULL) {
	return E_ALLOC;
    } else {
	tmp[n-1] = fns;
	*pfs = tmp;
	*pnf = n;
	return 0;
    }
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
	fns->rettype = GRETL_TYPE_VOID;
    }

    if (gretl_xml_get_prop_as_bool(node, "private")) {
	fns->flags |= UFUN_PRIVATE;
    }

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "params")) {
            fns->params = func_read_params(cur, doc, fns->name,
                                           &fns->n_params, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "code")) {
            err = !gretl_xml_node_get_trimmed_string(cur, doc, &fns->code);
	}
	cur = cur->next;
    }

    if (err) {
	funcshell_destroy(fns);
    } else {
	err = pkgshell_attach(fns, pks);
    }

    return err;
}

static pkgshell *real_read_pkgshell (xmlDocPtr doc,
				     xmlNodePtr node,
				     const char *fname,
                                     PRN *prn, int *err)
{
    xmlNodePtr cur;
    pkgshell *pks;
    char *tmp;

    pks = calloc(1, sizeof *pks);
    if (pks == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    if (gretl_xml_get_prop_as_string(node, "name", &tmp)) {
        pprintf(prn, "# hansl code from package %s", tmp);
        free(tmp);
    } else {
        *err = E_DATA;
        free(pks);
        return NULL;
    }

    cur = node->xmlChildrenNode;

    while (cur != NULL && !*err) {
        if (!xmlStrcmp(cur->name, (XUC) "version")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &tmp);
            pprintf(prn, " %s", tmp);
            free(tmp);
        } else if (!xmlStrcmp(cur->name, (XUC) "date")) {
            gretl_xml_node_get_trimmed_string(cur, doc, &tmp);
            pprintf(prn, " (%s)", tmp);
            free(tmp);
        } else if (!xmlStrcmp(cur->name, (XUC) "gretl-function")) {
	    *err = read_funcshell_from_xml(cur, doc, pks);
	}
	cur = cur->next;
    }

    pputs(prn, "\n\n");

    return pks;
}

static pkgshell *read_pkg_as_shell (const char *fname,
                                    PRN *prn, int *err)
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
	    pks = real_read_pkgshell(doc, cur, fname, prn, err);
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

static void print_pkgshell_code (pkgshell *pks, PRN *prn)
{
    funcshell *fns;
    int i;

    if (pks->n_priv > 0) {
        pputs(prn, "# private functions\n\n");
        for (i=0; i<pks->n_priv; i++) {
            fns = pks->priv[i];
            print_function_start(fns->name, fns->rettype,
                                 fns->params, fns->n_params,
                                 prn);
            pputs(prn, fns->code);
            pputs(prn, "\nend function\n\n");
        }
    }

    if (pks->n_pub > 0) {
        pputs(prn, "# public functions\n\n");
        for (i=0; i<pks->n_pub; i++) {
            fns = pks->pub[i];
            print_function_start(fns->name, fns->rettype,
                                 fns->params, fns->n_params,
                                 prn);
            pputs(prn, fns->code);
            pputs(prn, "\nend function\n\n");
        }
    }
}

/* This function is designed for printing the hansl code
   from a function package, preserving any comments and
   vertical space in the gfn "code" nodes. We do a minimal
   read of the package from @fname and output the code to
   @prn.
*/

int print_all_gfn_code (const char *fname, PRN *prn)
{
    pkgshell *pks = NULL;
    int err = 0;

    pks = read_pkg_as_shell(fname, prn, &err);

    if (!err) {
	print_pkgshell_code(pks, prn);
    }

    pkgshell_destroy(pks);

    return err;
}
