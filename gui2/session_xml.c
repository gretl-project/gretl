
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

    session_file_make_path(fullname, fname);
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
	    err = maybe_rewrite_gp_file(fullname);
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
	    session_file_make_path(oldname, fname);
	    tmp = g_strdup_printf("graph.%d", gnum);
	    session_file_make_path(newname, tmp);
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
    int inpage = 0;
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
	    if (gretl_xml_get_prop_as_int(cur, "inpage", &inpage)) {
		graph_page_add_file((const char *) fname); /* FIXME path? */
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
#if 1
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
	    session_file_make_path(fullname, (const char *) fname);
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
	    my_filename_from_utf8(sinfo->datafile);
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
	errbox("%d session object(s) could not be rebuilt", object_errs);
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
    session_file_make_path(fullname, "functions.xml");
    return write_session_functions_file(fullname);
}

static int write_settings_file (char *fullname)
{
    session_file_make_path(fullname, "settings.inp");
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

static int write_session_xml (const char *datname)
{
    MODEL *pmod;
    char fname[MAXLEN];
    char tmpname[MAXLEN];
    char *objname, *xmlname;
    FILE *fp, *fq;
    int nmodels;
    int tabmodels;
    int i, modnum;
    int err = 0;

    /* we should be in dotdir already when this is called */

    sprintf(fname, "%s%csession.xml", session.dirname, SLASH);
    fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	fprintf(stderr, " write_session_xml: failed on '%s'\n", fname);
	file_write_errbox(fname);
	return E_FOPEN;
    }

    gretl_xml_header(fp);

    if (*datname != '\0') {
	/* ensure UTF-8 inside XML file */
	gchar *trname = my_filename_to_utf8(datname);

	fprintf(fp, "<gretl-session datafile=\"%s\" date=\"%s\">\n", 
		trname, print_today());
	g_free(trname);
    } else {
	fprintf(fp, "<gretl-session date=\"%s\">\n", print_today());
    }

    if (data_status) {
	fprintf(fp, " <sample t1=\"%d\" t2=\"%d\"/>\n", dataset->t1, dataset->t2);
	write_datainfo_submask(dataset, fp);
    }

    nmodels = session.nmodels;
    tabmodels = model_table_n_models();

    for (i=0; i<tabmodels; i++) {
	pmod = model_table_model_by_index(i);
	if (!model_in_session(pmod)) {
	    nmodels++;
	}
    }

    fprintf(fp, " <models count=\"%d\">\n", nmodels);

    modnum = 1;

    for (i=0; i<session.nmodels && !err; i++) {
	int type = session.models[i]->type;
	void *ptr = session.models[i]->ptr;
	SavedObjectFlags sflags;
	int tablepos = 0;

	sprintf(tmpname, "%s%cmodel.%d", session.dirname, SLASH, modnum);
	fq = gretl_fopen(tmpname, "w");

	if (fq == NULL) {
	    file_write_errbox(tmpname);
	    err = E_FOPEN;
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
		fprintf(fp, "  <session-model name=\"%s\" fname=\"%s\" type=\"%d\" tablepos=\"%d\"/>\n", 
			xmlname, tmpname, type, tablepos);
	    } else {
		fprintf(fp, "  <session-model name=\"%s\" fname=\"%s\" type=\"%d\"/>\n", 
			xmlname, tmpname, type);
	    }	
	    if (xmlname != objname) {
		free(xmlname);
	    }
	    gretl_xml_header(fq);
	    sflags = model_save_flags(ptr, type);
	    if (type == GRETL_OBJ_EQN) {
		gretl_model_serialize(ptr, sflags, fq);
	    } else if (type == GRETL_OBJ_VAR) {
		gretl_VAR_serialize(ptr, sflags, fq);
	    } else if (type == GRETL_OBJ_SYS) {
		equation_system_serialize(ptr, sflags, fq);
	    }
	    fclose(fq);
	}
    }

    for (i=0; i<tabmodels && !err; i++) {
	pmod = model_table_model_by_index(i);
	if (!model_in_session(pmod)) {
	    sprintf(tmpname, "%s%cmodel.%d", session.dirname, SLASH, modnum);
	    fq = gretl_fopen(tmpname, "w");
	    if (fq == NULL) {
		file_write_errbox(tmpname);
		err = E_FOPEN;
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
		fprintf(fp, "  <session-model name=\"%s\" fname=\"%s\" type=\"%d\" tablepos=\"%d\"/>\n", 
			(xmlname != NULL)? xmlname : "none", tmpname, GRETL_OBJ_EQN, tablepos);
		if (xmlname != NULL && xmlname != objname) {
		    free(xmlname);
		}
		gretl_xml_header(fq);
		gretl_model_serialize(pmod, model_save_flags(pmod, GRETL_OBJ_EQN), 
				      fq);
		fclose(fq);
	    }
	}
    }

    if (err) {
	fclose(fp);
	gretl_remove(fname);
	return err;
    }

    fputs(" </models>\n", fp);

    fprintf(fp, " <graphs count=\"%d\">\n", session.ngraphs);

    for (i=0; i<session.ngraphs; i++) {
	objname = session.graphs[i]->name;
	xmlname = get_xmlname(objname, &err);
	if (err) {
	    break;
	}
	fprintf(fp, "  <session-graph name=\"%s\" fname=\"%s\" type=\"%d\"", 
		xmlname, session.graphs[i]->fname, session.graphs[i]->type);
	if (xmlname != objname) {
	    free(xmlname);
	}
	if (in_graph_page(session.graphs[i]->fname)) {
	    fputs(" inpage=\"1\"/>\n", fp);
	} else {
	    fputs("/>\n", fp);
	}
    } 
    fputs(" </graphs>\n", fp);

    fprintf(fp, " <texts count=\"%d\">\n", session.ntexts);

    for (i=0; i<session.ntexts; i++) {
	objname = session.texts[i]->name;
	xmlname = get_xmlname(objname, &err);
	if (err) {
	    break;
	}	
	fprintf(fp, "  <session-text name=\"%s\">", xmlname);
	if (xmlname != objname) {
	    free(xmlname);
	}	
	gretl_xml_put_raw_string(session.texts[i]->buf, fp);
	fputs("</session-text>\n", fp);
    }    
    fputs(" </texts>\n", fp);

    if (session.notes != NULL) {
	if (session.show_notes) {
	    fputs("<notes auto-show=\"true\">", fp);
	} else {
	    fputs("<notes>", fp);
	}
	gretl_xml_put_raw_string(session.notes, fp);
	fputs("</notes>\n", fp);
    } 

    fputs("</gretl-session>\n", fp);

    fclose(fp);

    serialize_user_vars(session.dirname);

    maybe_write_function_file(tmpname);
    write_settings_file(tmpname);

    return 0;
}
