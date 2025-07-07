
/* addendum to session.c, for handling the saving of session info to
   an XML file, and the re-building of a session from same.
*/

#include "usermat.h"
#include "boxplots.h"
#include "libset.h"

static int check_graph_file (const char *fname, int type)
{
    char fullname[MAXLEN];
    FILE *fp;
    int err = 0;

    session_file_make_path(fullname, fname, NULL);
    fp = gretl_fopen(fullname, "r");

    if (fp == NULL) {
	file_read_errbox(fname);
	err = 1;
    } else {
	fprintf(stderr, "Opened '%s' OK\n", fullname);
	if (type == GRETL_OBJ_PLOT) {
	    char line[32] = {0};

	    if (fgets(line, sizeof line, fp) != NULL) {
		if (!strncmp(line, "# boxplot generated", 19)) {
		    fprintf(stderr, "Ignoring old boxplot file\n");
		    err = 1;
		}
	    } else {
		err = 1;
	    }
	    fclose(fp);
	} else {
	    fclose(fp);
	}
    }

    return err;
}

/* Arrange things so that when a session file is opened, any graph
   files are numbered consecutively, starting at 1.  This avoids a
   situation where adding a graph to an existing session results in
   over-writing an existing graph file.  (The graph files in a saved
   session may not be numbered consecutively if some graphs were
   initially saved, then deleted; and it's easier to fix this on
   opening a session file than on saving one, since on opening we're
   constructing the graph list de novo.)
*/

static void normalize_graph_filename (char *fname, int gnum)
{
    char s[16];
    int i;

    if (sscanf(fname, "%15[^.].%d", s, &i) == 2) {
	char oldname[MAXLEN];
	char newname[MAXLEN];
	gchar *tmp = NULL;

	if (i != gnum && (!strcmp(s, "graph") || !strcmp(s, "plot"))) {
	    session_file_make_path(oldname, fname, NULL);
	    tmp = g_strdup_printf("graph.%d", gnum);
	    session_file_make_path(newname, tmp, NULL);
	}

	if (tmp != NULL) {
	    gretl_rename(oldname, newname);
	    strcpy(fname, tmp);
	    g_free(tmp);
	}
    }
}

static int restore_session_graphs (xmlNodePtr node)
{
    xmlNodePtr cur;
    int gnum = 0;
    int errs = 0;

    /* reset prior to parsing */
    session.ngraphs = 0;

    cur = node->xmlChildrenNode;

    fprintf(stderr, "Doing restore_session_graphs\n");

    while (cur != NULL) {
	SESSION_GRAPH *sg;
	xmlChar *name = NULL;
	xmlChar *fname = NULL;
	int type, err = 0;

	name = xmlGetProp(cur, (XUC) "name");
	if (name == NULL) {
	    err = 1;
	}

	if (!err && !gretl_xml_get_prop_as_int(cur, "type", &type)) {
	    err = 1;
	}

	if (!err) {
	    fname = xmlGetProp(cur, (XUC) "fname");
	    if (fname == NULL) {
		err = 1;
	    } else {
		fprintf(stderr, "checking '%s'\n", name);
		err = check_graph_file((const char *) fname, type);
		if (!err) {
		    normalize_graph_filename((char *) fname, ++gnum);
		}
	    }
	}

	if (!err) {
	    sg = session_append_graph((const char *) name,
				      (const char *) fname,
				      type);
	    err = (sg == NULL);
	}

	if (!err) {
	    int inpage, has_datafile;

	    if (gretl_xml_get_prop_as_int(cur, "inpage", &inpage)) {
		graph_page_add_file((const char *) fname); /* FIXME path? */
	    }
	    if (gretl_xml_get_prop_as_int(cur, "has_datafile", &has_datafile)) {
		sg->has_datafile = 1;
	    }
	}

	if (err) {
	    errs++;
	}

	free(name);
	free(fname);

	cur = cur->next;
    }

    return errs;
}

static int restore_session_texts (xmlNodePtr node, xmlDocPtr doc)
{
    xmlNodePtr cur;
    int errs = 0;

    session.ntexts = 0;

    cur = node->xmlChildrenNode;

    while (cur != NULL) {
	xmlChar *name = NULL;
	xmlChar *buf = NULL;
	int err = 0;

	name = xmlGetProp(cur, (XUC) "name");
	if (name == NULL) {
	    err = 1;
	} else {
	    buf = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	    if (buf == NULL) {
		err = 1;
	    } else {
		err = session_append_text((const char *) name,
					  (char *) buf);
	    }
	}

	if (err) {
	    errs++;
	}

	free(name);

	cur = cur->next;
    }

    return errs;
}

static int data_submask_from_xml (xmlNodePtr node, xmlDocPtr doc,
				  struct sample_info *sinfo)
{
    char *mask;
    int err;

    err = gretl_xml_get_submask(node, doc, &mask);
    if (!err) {
	sinfo->mask = mask;
    }

    return err;
}

static int data_restrict_from_xml (xmlNodePtr node, xmlDocPtr doc,
				   struct sample_info *sinfo)
{
    char *s = NULL;

    if (gretl_xml_node_get_trimmed_string(node, doc, &s)) {
	sinfo->restriction = s;
    }

    return 0;
}

static int rebuild_session_model (const char *fname,
				  const char *name,
				  GretlObjType type,
				  int tablepos)
{
    gpointer ptr = NULL;
    SavedObjectFlags flags = 0;
    xmlDocPtr doc;
    xmlNodePtr node;
    int iflag, err;

    err = gretl_xml_open_doc_root(fname,
				  (type == GRETL_OBJ_EQN)? "gretl-model" :
				  (type == GRETL_OBJ_VAR)? "gretl-VAR" :
				  "gretl-equation-system", &doc, &node);
    if (err) {
	return err;
    }

    if (gretl_xml_get_prop_as_int(node, "saveflags", &iflag)) {
	flags = iflag;
    } else {
	flags = IN_GUI_SESSION | IN_NAMED_STACK;
    }

    if (type == GRETL_OBJ_EQN) {
	ptr = gretl_model_from_XML(node, doc, dataset, &err);
    } else if (type == GRETL_OBJ_VAR) {
	ptr = gretl_VAR_from_XML(node, doc, dataset, &err);
    } else {
	ptr = equation_system_from_XML(node, doc, &err);
    }

    xmlFreeDoc(doc);

    if (!err && ptr != NULL) {
#if SESSION_DEBUG
	fprintf(stderr, "rebuild_session_model: IN_GUI_SESSION %s, IN_MODEL_TABLE %s\n"
		" IN_NAMED_STACK %s, IS_LAST_MODEL %s\n",
		(flags & IN_GUI_SESSION)? "yes" : "no",
		(flags & IN_MODEL_TABLE)? "yes" : "no",
		(flags & IN_NAMED_STACK)? "yes" : "no",
		(flags & IS_LAST_MODEL)? "yes" : "no");
#endif

	if (flags & IN_GUI_SESSION) {
	    SESSION_MODEL *smod;

	    smod = session_model_new(ptr, name, type);
	    if (smod == NULL) {
		fprintf(stderr, "error from session_model_new\n");
		err = E_ALLOC;
	    } else {
		err = session_append_model(smod);
		fprintf(stderr, "error %d from session_append_model\n", err);
	    }
	}

	if (!err && (flags & IN_MODEL_TABLE)) {
	    add_to_model_table(ptr, MODEL_ADD_BY_CMD, tablepos, NULL);
	}

	if (!err && (flags & IN_NAMED_STACK)) {
	    err = gretl_stack_object(ptr, type);
	}

	if (!err && (flags & IS_LAST_MODEL)) {
	    set_as_last_model(ptr, type);
	}
    }

    /* need to clean up on error here (also: clean up XML parser?) */

    return err;
}

static int restore_session_models (xmlNodePtr node, xmlDocPtr doc)
{
    char fullname[MAXLEN];
    xmlNodePtr cur;
    int errs = 0;

    /* reset prior to parsing */
    session.nmodels = 0;

    cur = node->xmlChildrenNode;

    while (cur != NULL) {
	xmlChar *fname = NULL;
	xmlChar *name = NULL;
	int type = GRETL_OBJ_EQN;
	int tablepos = 0;
	int err = 0;

	fname = xmlGetProp(cur, (XUC) "fname");
	if (fname == NULL) {
	    err = 1;
	} else {
	    name = xmlGetProp(cur, (XUC) "name");
	    if (name == NULL) {
		err = 1;
	    }
	}

	if (!err) {
	    session_file_make_path(fullname, (const char *) fname, NULL);
	    gretl_xml_get_prop_as_int(cur, "type", &type);
	    if (type == GRETL_OBJ_EQN) {
		gretl_xml_get_prop_as_int(cur, "tablepos", &tablepos);
	    }
	    fprintf(stderr, "model file: fname='%s', type=%d\n", fullname, type);
	    err = rebuild_session_model(fullname, (const char *) name,
					type, tablepos);
	}

	if (!err) {
	    model_count_plus();
	} else {
	    fprintf(stderr, "rebuild_session_model: failed on %s (err = %d)\n",
		    fullname, err);
	    errs++;
	}

	free(fname);
	free(name);

	cur = cur->next;
    }

#if SESSION_DEBUG
    fprintf(stderr, "restore_session_models: %d errors\n", errs);
#endif

    return errs;
}

#if SESSION_BUNDLE

static gretl_bundle *session_model_to_bundle (const char *fname,
					      const char *name,
					      GretlObjType otype,
					      DATASET *dset,
					      int *err)
{
    gretl_bundle *mb = NULL;
    void *ptr = NULL;
    xmlDocPtr doc;
    xmlNodePtr node;

    *err = gretl_xml_open_doc_root(fname,
				   (otype == GRETL_OBJ_EQN)? "gretl-model" :
				   (otype == GRETL_OBJ_VAR)? "gretl-VAR" :
				   "gretl-equation-system", &doc, &node);
    if (*err) {
	return NULL;
    }

    if (otype == GRETL_OBJ_EQN) {
	ptr = gretl_model_from_XML(node, doc, dset, err);
    } else if (otype == GRETL_OBJ_VAR) {
	ptr = gretl_VAR_from_XML(node, doc, dset, err);
    } else {
	ptr = equation_system_from_XML(node, doc, err);
    }

    xmlFreeDoc(doc);

    if (!*err && ptr != NULL) {
	if (otype == GRETL_OBJ_EQN) {
	    mb = bundle_from_model(ptr, dset, err);
	} else {
	    mb = bundle_from_system(ptr, otype, dset, err);
	}
	/* add the model's name to its bundle */
	gretl_bundle_set_data(mb, "name", (char *) name, GRETL_TYPE_STRING, 0);
    }

    if (ptr != NULL) {
	/* the structs obtained above are disposable */
	if (otype == GRETL_OBJ_EQN) {
	    gretl_model_free(ptr);
	} else if (otype == GRETL_OBJ_VAR) {
	    gretl_VAR_free(ptr);
	} else {
	    equation_system_destroy(ptr);
	}
    }

    return mb;
}

static int session_models_to_bundles (xmlNodePtr node,
				      xmlDocPtr doc,
				      gretl_array *am,
				      DATASET *dset,
				      const char *sdir)
{
    char fullname[MAXLEN];
    xmlNodePtr cur;
    int i = 0;
    int errs = 0;

    cur = node->xmlChildrenNode;

    while (cur != NULL) {
	gretl_bundle *mb;
	xmlChar *fname = NULL;
	xmlChar *name = NULL;
	int type = GRETL_OBJ_EQN;
	int err = 0;

	fname = xmlGetProp(cur, (XUC) "fname");
	name = xmlGetProp(cur, (XUC) "name");
	if (fname == NULL || name == NULL) {
	    err = 1;
	} else {
	    session_file_make_path(fullname, (const char *) fname, sdir);
	    gretl_xml_get_prop_as_int(cur, "type", &type);
#if SESSION_DEBUG
	    fprintf(stderr, "model file: fname='%s', name='%s', type=%d\n",
		    fname, name, type);
#endif
	    mb = session_model_to_bundle(fullname, (const char *) name,
					 type, dset, &err);
	}
	if (!err) {
	    gretl_array_set_data(am, i++, mb);
	}
	if (err) {
	    fprintf(stderr, "session_model_to_bundle: failed on %s (err = %d)\n",
		    fullname, err);
	    errs++;
	}

	free(fname);
	free(name);

	cur = cur->next;
    }

#if SESSION_DEBUG
    fprintf(stderr, "session_models_to_bundles: %d errors\n", errs);
#endif

    return errs;
}

static gretl_bundle *session_xml_to_bundle (const char *fname,
					    const char *sdir,
					    DATASET *dset,
					    int *err)
{
    gretl_bundle *sb = NULL;
    gretl_array *am = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    xmlChar *tmp;
    int obj_errs = 0;

    *err = gretl_xml_open_doc_root(fname, "gretl-session", &doc, &cur);
    if (*err) {
	gui_errmsg(*err);
	return NULL;
    }

    sb = gretl_bundle_new();
    if (sb == NULL) {
	*err = E_ALLOC;
    }

    /* FIXME add some dataset info? */

    if (!*err) {
	/* Now walk the tree */
	cur = cur->xmlChildrenNode;
	while (cur != NULL && !*err) {
	    if (!xmlStrcmp(cur->name, (XUC) "models")) {
		int n_models = 0;

		tmp = xmlGetProp(cur, (XUC) "count");
		if (tmp != NULL) {
		    n_models = atoi((const char *) tmp);
		    free(tmp);
		}
#if SESSION_DEBUG
		fprintf(stderr, "session_xml_to_bundle: n_models = %d\n", n_models);
#endif
		if (n_models > 0) {
		    am = gretl_array_new(GRETL_TYPE_BUNDLES, n_models, err);
		    if (!*err) {
			obj_errs += session_models_to_bundles(cur, doc, am, dset, sdir);
			*err = gretl_bundle_donate_data(sb, "models", am,
							GRETL_TYPE_ARRAY, 0);
		    }
		}
		break;
	    }
	    cur = cur->next;
	}
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    if (!*err && obj_errs > 0) {
	errbox_printf("%d session object(s) could not be rebuilt", obj_errs);
    }

    return sb;
}

#endif /* SESSION_BUNDLE */

/* peek inside the session file and retrieve the name of the
   data file, only */

static int get_session_datafile_name (const char *fname, struct sample_info *sinfo,
				      int *nodata)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    xmlChar *tmp;
    int err = 0;

    err = gretl_xml_open_doc_root(fname, "gretl-session", &doc, &cur);
    if (err) {
	gui_errmsg(err);
	return 1;
    }

    /* read datafile attribute, if present */
    tmp = xmlGetProp(cur, (XUC) "datafile");
    if (tmp != NULL) {
	if (!strcmp((const char *) tmp, "none")) {
	    *nodata = 1;
	} else {
	    strcpy(sinfo->datafile, (char *) tmp);
	}
	free(tmp);
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    return err;
}

/* (having previously grabbed the data file name) get the rest
   of the info from session.xml */

static int
read_session_xml (const char *fname, struct sample_info *sinfo)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    xmlChar *tmp;
    int object_errs = 0;
    int err = 0;

    err = gretl_xml_open_doc_root(fname, "gretl-session", &doc, &cur);
    if (err) {
	gui_errmsg(err);
	return err;
    }

    /* Now walk the tree */
    cur = cur->xmlChildrenNode;
    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "sample")) {
	    tmp = xmlGetProp(cur, (XUC) "t1");
	    if (tmp != NULL) {
		sinfo->t1 = atoi((const char *) tmp);
		free(tmp);
	    } else {
		err = 1;
	    }
	    tmp = xmlGetProp(cur, (XUC) "t2");
	    if (tmp != NULL) {
		sinfo->t2 = atoi((const char *) tmp);
		free(tmp);
	    } else {
		err = 1;
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "submask")) {
	    err = data_submask_from_xml(cur, doc, sinfo);
	} else if (!xmlStrcmp(cur->name, (XUC) "restriction")) {
	    err = data_restrict_from_xml(cur, doc, sinfo);
	} else if (!xmlStrcmp(cur->name, (XUC) "resample")) {
	    tmp = xmlGetProp(cur, (XUC) "seed");
	    if (tmp != NULL) {
		sinfo->seed = (unsigned) atoi((const char *) tmp);
		free(tmp);
	    }
	    tmp = xmlGetProp(cur, (XUC) "n");
	    if (tmp != NULL) {
		sinfo->resample_n = atoi((const char *) tmp);
		free(tmp);
	    }
	    if (sinfo->resample_n <= 0) {
		sinfo->seed = 0;
		sinfo->resample_n = 0;
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "models")) {
	    tmp = xmlGetProp(cur, (XUC) "count");
	    if (tmp != NULL) {
		session.nmodels = atoi((const char *) tmp);
		free(tmp);
		if (session.nmodels > 0) {
		    object_errs += restore_session_models(cur, doc);
		}
	    }
        } else if (!xmlStrcmp(cur->name, (XUC) "graphs")) {
	    tmp = xmlGetProp(cur, (XUC) "count");
	    if (tmp != NULL) {
		session.ngraphs = atoi((const char *) tmp);
		free(tmp);
		if (session.ngraphs > 0) {
		    object_errs += restore_session_graphs(cur);
		}
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "texts")) {
	    tmp = xmlGetProp(cur, (XUC) "count");
	    if (tmp != NULL) {
		session.ntexts = atoi((const char *) tmp);
		free(tmp);
		if (session.ntexts > 0) {
		    object_errs += restore_session_texts(cur, doc);
		}
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "notes")) {
	    session.notes =
		(char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	    if (session.notes == NULL) {
		object_errs++;
	    } else {
		tmp = xmlGetProp(cur, (XUC) "auto-show");
		if (tmp != NULL) {
		    session.show_notes = 1;
		    free(tmp);
		}
	    }
	}
	if (!err) {
	    cur = cur->next;
	}
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
    }

    if (!err && object_errs > 0) {
	errbox_printf("%d session object(s) could not be rebuilt", object_errs);
    }

    return err;
}

static int maybe_read_functions_file (const char *fname)
{
    FILE *fp;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	/* nothing to be read */
	return 0;
    }

    fclose(fp);

    return read_session_functions_file(fname);
}

static int maybe_read_settings_file (const char *fname)
{
    FILE *fp;

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	/* nothing to be read */
	return 0;
    }

    fclose(fp);

    return libset_read_script(fname);
}

static int model_in_session (const void *ptr)
{
    int i;

    for (i=0; i<session.nmodels; i++) {
	if (session.models[i]->ptr == ptr) {
	    return 1;
	}
    }

    return 0;
}

static SavedObjectFlags model_save_flags (const void *ptr,
					  GretlObjType type)
{
    SavedObjectFlags flags = 0;

    if (model_in_session(ptr)) {
	flags |= IN_GUI_SESSION;
    }

    if (object_is_on_stack(ptr)) {
	flags |= IN_NAMED_STACK;
    }

    if (get_last_model(NULL) == ptr) {
	flags |= IS_LAST_MODEL;
    }

    if (type == GRETL_OBJ_EQN && in_model_table(ptr)) {
	flags |= IN_MODEL_TABLE;
    }

    return flags;
}

static int maybe_write_function_file (char *fullname)
{
    session_file_make_path(fullname, "functions.xml", NULL);
    return write_loaded_functions_file(fullname, OPT_NONE);
}

static int write_settings_file (char *fullname)
{
    session_file_make_path(fullname, "settings.inp", NULL);
    return libset_write_script(fullname);
}

static char *get_xmlname (char *objname, int *err)
{
    char *ret;

    if (gretl_xml_validate(objname)) {
	ret = objname;
    } else {
	ret = gretl_xml_encode(objname);
	if (ret == NULL) {
	    *err = E_ALLOC;
	}
    }

    return ret;
}

static int session_graph_wanted (const char *fname)
{
    int i;

    for (i=0; i<session.ngraphs; i++) {
	if (!strcmp(session.graphs[i]->fname, fname)) {
	    return 1;
	}
    }

    return 0;
}

/* on re-saving a session, avoid keeping model files for
   models that have been dropped from the session, and
   similarly for graph files
*/

static void trash_old_session_files (const char *path)
{
    GDir *sdir = gretl_opendir(path);

    if (sdir != NULL) {
	const gchar *dname;
	char tmp[2*MAXLEN];
	int fnum;

	while ((dname = g_dir_read_name(sdir)) != NULL) {
	    if (!strncmp(dname, "model.", 6)) {
		fnum = atoi(dname + 6);
		sprintf(tmp, "model.%d", fnum);
		if (!strcmp(dname, tmp)) {
		    sprintf(tmp, "%s%cmodel.%d", path, SLASH, fnum);
		    gretl_remove(tmp);
		}
	    } else if (!strncmp(dname, "graph.", 6)) {
		fnum = atoi(dname + 6);
		sprintf(tmp, "graph.%d", fnum);
		if (!strcmp(dname, tmp) &&
		    !session_graph_wanted(tmp)) {
		    sprintf(tmp, "%s%cgraph.%d", path, SLASH, fnum);
		    gretl_remove(tmp);
		}
	    }
	}
	g_dir_close(sdir);
    }
}

static int write_session_xml (const char *datname)
{
    MODEL *pmod;
    char fname[2*MAXLEN];
    char tmpname[2*MAXLEN];
    char *objname, *xmlname;
    PRN *prn;
    int nmodels;
    int tabmodels;
    int i, modnum;
    int err = 0;

    /* we should be in dotdir already when this is called */

    sprintf(fname, "%s%csession.xml", session.dirname, SLASH);
    prn = gretl_print_new_with_filename(fname, &err);

    if (err) {
	fprintf(stderr, " write_session_xml: failed on '%s'\n", fname);
	file_write_errbox(fname);
	return err;
    }

    gretl_xml_header(prn);

    if (*datname != '\0') {
	if (gretl_xml_validate(datname)) {
	    pprintf(prn, "<gretl-session datafile=\"%s\" date=\"%s\">\n",
		    datname, print_today());
	} else {
	    char *xstr = gretl_xml_encode(datname);

	    pprintf(prn, "<gretl-session datafile=\"%s\" date=\"%s\">\n",
		    xstr, print_today());
	    free(xstr);
	}
    } else {
	pprintf(prn, "<gretl-session date=\"%s\">\n", print_today());
    }

    if (data_status) {
	pprintf(prn, " <sample t1=\"%d\" t2=\"%d\"/>\n", dataset->t1, dataset->t2);
	write_dataset_submask(dataset, prn);
    }

    nmodels = session.nmodels;
    tabmodels = model_table_n_models();

    for (i=0; i<tabmodels; i++) {
	pmod = model_table_model_by_index(i);
	if (!model_in_session(pmod)) {
	    nmodels++;
	}
    }

    trash_old_session_files(session.dirname);

    pprintf(prn, " <models count=\"%d\">\n", nmodels);

    modnum = 1;

    for (i=0; i<session.nmodels && !err; i++) {
	int type = session.models[i]->type;
	void *ptr = session.models[i]->ptr;
	SavedObjectFlags sflags;
	int tablepos = 0;
	PRN *pm;

	sprintf(tmpname, "%s%cmodel.%d", session.dirname, SLASH, modnum);
	pm = gretl_print_new_with_filename(tmpname, &err);

	if (err) {
	    file_write_errbox(tmpname);
	} else {
	    sprintf(tmpname, "model.%d", modnum++);
	    objname = session.models[i]->name;
	    xmlname = get_xmlname(objname, &err);
	    if (err) {
		break;
	    }
	    if (type == GRETL_OBJ_EQN) {
		tablepos = model_table_position(ptr);
	    }
	    if (tablepos > 0) {
		pprintf(prn, "  <session-model name=\"%s\" fname=\"%s\" type=\"%d\" tablepos=\"%d\"/>\n",
			xmlname, tmpname, type, tablepos);
	    } else {
		pprintf(prn, "  <session-model name=\"%s\" fname=\"%s\" type=\"%d\"/>\n",
			xmlname, tmpname, type);
	    }
	    if (xmlname != objname) {
		free(xmlname);
	    }
	    gretl_xml_header(pm);
	    sflags = model_save_flags(ptr, type);
	    gretl_push_c_numeric_locale();
	    if (type == GRETL_OBJ_EQN) {
		gretl_model_serialize(ptr, sflags, pm);
	    } else if (type == GRETL_OBJ_VAR) {
		gretl_VAR_serialize(ptr, sflags, pm);
	    } else if (type == GRETL_OBJ_SYS) {
		equation_system_serialize(ptr, sflags, pm);
	    }
	    gretl_pop_c_numeric_locale();
	    gretl_print_destroy(pm);
	}
    }

    for (i=0; i<tabmodels && !err; i++) {
	PRN *pm;

	pmod = model_table_model_by_index(i);
	if (!model_in_session(pmod)) {
	    sprintf(tmpname, "%s%cmodel.%d", session.dirname, SLASH, modnum);
	    pm = gretl_print_new_with_filename(tmpname, &err);

	    if (err) {
		file_write_errbox(tmpname);
	    } else {
		int tablepos = model_table_position(pmod);

		if (pmod->name == NULL) {
		    objname = xmlname = NULL;
		} else {
		    objname = pmod->name;
		    xmlname = get_xmlname(objname, &err);
		}
		if (err) {
		    break;
		}
		sprintf(tmpname, "model.%d", modnum++);
		pprintf(prn, "  <session-model name=\"%s\" fname=\"%s\" type=\"%d\" tablepos=\"%d\"/>\n",
			(xmlname != NULL)? xmlname : "none", tmpname, GRETL_OBJ_EQN, tablepos);
		if (xmlname != NULL && xmlname != objname) {
		    free(xmlname);
		}
		gretl_xml_header(pm);
		gretl_push_c_numeric_locale();
		gretl_model_serialize(pmod, model_save_flags(pmod, GRETL_OBJ_EQN),
				      pm);
		gretl_pop_c_numeric_locale();
		gretl_print_destroy(pm);
	    }
	}
    }

    if (err) {
	gretl_print_destroy(prn);
	gretl_remove(fname);
	return err;
    }

    pputs(prn, " </models>\n");

    pprintf(prn, " <graphs count=\"%d\">\n", session.ngraphs);

    for (i=0; i<session.ngraphs; i++) {
	objname = session.graphs[i]->name;
	xmlname = get_xmlname(objname, &err);
	if (err) {
	    break;
	}
	pprintf(prn, "  <session-graph name=\"%s\" fname=\"%s\" type=\"%d\"",
		xmlname, session.graphs[i]->fname, session.graphs[i]->type);
	if (xmlname != objname) {
	    free(xmlname);
	}
	if (in_graph_page(session.graphs[i]->fname)) {
	    pputs(prn, " inpage=\"1\"");
	}
	if (session.graphs[i]->has_datafile) {
	    pputs(prn, " has_datafile=\"1\"");
	}
	pputs(prn, "/>\n");
    }
    pputs(prn, " </graphs>\n");

    pprintf(prn, " <texts count=\"%d\">\n", session.ntexts);

    for (i=0; i<session.ntexts; i++) {
	objname = session.texts[i]->name;
	xmlname = get_xmlname(objname, &err);
	if (err) {
	    break;
	}
	pprintf(prn, "  <session-text name=\"%s\">", xmlname);
	if (xmlname != objname) {
	    free(xmlname);
	}
	gretl_xml_put_string(session.texts[i]->buf, prn);
	pputs(prn, "</session-text>\n");
    }
    pputs(prn, " </texts>\n");

    if (session.notes != NULL) {
	if (session.show_notes) {
	    pputs(prn, "<notes auto-show=\"true\">");
	} else {
	    pputs(prn, "<notes>");
	}
	gretl_xml_put_string(session.notes, prn);
	pputs(prn, "</notes>\n");
    }

    pputs(prn, "</gretl-session>\n");

    gretl_print_destroy(prn);

    serialize_user_vars(session.dirname);

    maybe_write_function_file(tmpname);
    write_settings_file(tmpname);

    return 0;
}
