
/* addendum to session.c, for handling the saving of session info to
   an XML file, and the re-building of a session from same.
*/

static int check_graph_file (const char *fname)
{
    char fullname[MAXLEN];
    FILE *fp;
    int err = 0;

    session_file_make_path(fullname, fname);
    fp = gretl_fopen(fullname, "r");
    if (fp == NULL) {
	errbox(_("Warning: couldn't open graph file %s"), fname);
	err = 1;
    } else {
	fclose(fp);
    }

    return err;
}

static int restore_session_graphs (xmlNodePtr node)
{
    xmlNodePtr cur;
    int i = 0;
    int err = 0;

    session.graphs = mymalloc(session.ngraphs * sizeof *session.graphs);
    if (session.graphs == NULL) {
	return E_ALLOC;
    }

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	xmlChar *name = NULL;
	xmlChar *fname = NULL;
	int type, ID;

	name = xmlGetProp(cur, (XUC) "name");
	if (name == NULL) {
	    err = 1;
	}

	if (!err) {
	    fname = xmlGetProp(cur, (XUC) "fname");
	    if (fname == NULL) {
		err = 1;
	    } else {
		err = check_graph_file((const char *) fname);
	    } 
	}

	if (!err && (!gretl_xml_get_prop_as_int(cur, "ID", &ID) ||
	    !gretl_xml_get_prop_as_int(cur, "type", &type))) {
	    err = 1;
	}

	if (!err) {
	    session.graphs[i] = session_graph_new((const char *) name, 
						  (const char *) fname, 
						  type);
	    if (session.graphs[i] == NULL) {
		err = 1;
	    } else {
		session.graphs[i]->ID = ID;
	    }
	}

	free(name);
	free(fname);

	cur = cur->next;
	i++;
    }

    return err;
}

static int restore_session_texts (xmlDocPtr doc, xmlNodePtr node)
{
    xmlNodePtr cur;
    int i = 0;
    int err = 0;

    session.texts = mymalloc(session.ntexts * sizeof *session.texts);
    if (session.texts == NULL) {
	return E_ALLOC;
    }

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	xmlChar *name = NULL;
	xmlChar *buf = NULL;

	name = xmlGetProp(cur, (XUC) "name");
	if (name == NULL) {
	    err = 1;
	} else {
	    session.texts[i] = session_text_new((const char *) name);
	    if (session.texts[i] == NULL) {
		err = 1;
	    }
	}

	if (!err) {
	    buf = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	    if (buf != NULL) {
		session.texts[i]->buf = (char *) buf;
	    } else {
		err = 1;
	    }
	}

	free(name);

	cur = cur->next;
	i++;
    }

    return err;
}

static MODEL *rebuild_session_model (const char *fname, int *err)
{
    MODEL *pmod;
    xmlDocPtr doc;
    xmlNodePtr node;
    xmlNodePtr cur;
    int got;

#if SESSION_DEBUG
    fprintf(stderr, "rebuild_session_model: trying to open '%s'\n",
	    fname);
#endif

    *err = gretl_xml_open_doc_root(fname, "gretl-model", &doc, &node);
    if (*err) {
	return NULL;
    }

    pmod = gretl_model_new();
    if (pmod == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

#if SESSION_DEBUG
    fprintf(stderr, "rebuild_session_model: allocated model at %p\n",
	    (void *) pmod);
#endif

    got = 0;
    got += gretl_xml_get_prop_as_int(node, "ID", &pmod->ID);
    got += gretl_xml_get_prop_as_int(node, "t1", &pmod->t1);
    got += gretl_xml_get_prop_as_int(node, "t2", &pmod->t2);
    got += gretl_xml_get_prop_as_int(node, "nobs", &pmod->nobs);
    got += gretl_xml_get_prop_as_int(node, "full_n", &pmod->full_n);
    got += gretl_xml_get_prop_as_int(node, "ncoeff", &pmod->ncoeff);
    got += gretl_xml_get_prop_as_int(node, "dfn", &pmod->dfn);
    got += gretl_xml_get_prop_as_int(node, "dfd", &pmod->dfd);
    got += gretl_xml_get_prop_as_int(node, "ifc", &pmod->ifc);
    got += gretl_xml_get_prop_as_int(node, "ci", &pmod->ci);
    got += gretl_xml_get_prop_as_int(node, "nwt", &pmod->nwt);
    got += gretl_xml_get_prop_as_int(node, "order", &pmod->order);
    got += gretl_xml_get_prop_as_int(node, "aux", &pmod->aux);

    if (got < 13) {
	*err = E_DATA;
	goto bailout;
    }

    gretl_push_c_numeric_locale();

    got = 0;
    got += gretl_xml_get_prop_as_double(node, "ess", &pmod->ess);
    got += gretl_xml_get_prop_as_double(node, "tss", &pmod->tss);
    got += gretl_xml_get_prop_as_double(node, "sigma", &pmod->sigma);
    got += gretl_xml_get_prop_as_double(node, "ess_wt", &pmod->ess_wt);
    got += gretl_xml_get_prop_as_double(node, "sigma_wt", &pmod->sigma_wt);
    got += gretl_xml_get_prop_as_double(node, "rsq", &pmod->rsq);
    got += gretl_xml_get_prop_as_double(node, "adjrsq", &pmod->adjrsq);
    got += gretl_xml_get_prop_as_double(node, "fstt", &pmod->fstt);
    got += gretl_xml_get_prop_as_double(node, "lnL", &pmod->lnL);
    got += gretl_xml_get_prop_as_double(node, "ybar", &pmod->ybar);
    got += gretl_xml_get_prop_as_double(node, "sdy", &pmod->sdy);

    got += gretl_xml_get_prop_as_double(node, "crit0", &pmod->criterion[0]);
    got += gretl_xml_get_prop_as_double(node, "crit1", &pmod->criterion[1]);
    got += gretl_xml_get_prop_as_double(node, "crit2", &pmod->criterion[2]);

    got += gretl_xml_get_prop_as_double(node, "dw", &pmod->dw);
    got += gretl_xml_get_prop_as_double(node, "rho", &pmod->rho);

    if (got < 16) {
	*err = E_DATA;
	gretl_pop_c_numeric_locale();
	goto bailout;
    }

    cur = node->xmlChildrenNode;

    while (cur != NULL && !*err) {
	if (!xmlStrcmp(cur->name, (XUC) "coeff")) {
	    pmod->coeff = gretl_xml_get_doubles_array(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "sderr")) {
	    pmod->sderr = gretl_xml_get_doubles_array(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "uhat")) {
	    pmod->uhat = gretl_xml_get_doubles_array(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "yhat")) {
	    pmod->yhat = gretl_xml_get_doubles_array(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "xpx")) {
	    pmod->xpx = gretl_xml_get_doubles_array(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "vcv")) {
	    pmod->vcv = gretl_xml_get_doubles_array(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "list")) {
	    pmod->list = gretl_xml_node_get_list(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "tests")) {
	    *err = attach_model_tests_from_xml(pmod, cur);
	}
	cur = cur->next;
    }

    gretl_pop_c_numeric_locale();

 bailout:

    xmlFreeDoc(doc);

    if (*err) {
	if (pmod != NULL) {
	    gretl_model_free(pmod);
	    pmod = NULL;
	}
    } else {
	gretl_object_ref(pmod, GRETL_OBJ_EQN);
    }

#if SESSION_DEBUG
    fprintf(stderr, "rebuild_session_model: returning with err = %d\n",
	    *err);
#endif
    
    /* need to clean up on error here (also: clean up XML parser?) */

    return pmod;
}

static int restore_session_models (xmlDocPtr doc, xmlNodePtr node)
{
    char fullname[MAXLEN];
    xmlNodePtr cur;
    int i, err = 0;

    session.models = mymalloc(session.nmodels * sizeof *session.models);
    if (session.models == NULL) {
	return E_ALLOC;
    }

#if SESSION_DEBUG
    fprintf(stderr, "rebuilding: session.nmodels = %d\n", session.nmodels);
#endif

    cur = node->xmlChildrenNode;
    i = 0;

    while (cur != NULL && !err) {
	xmlChar *fname = NULL;
	xmlChar *name = NULL;
	MODEL *pmod = NULL;

	fname = xmlGetProp(cur, (XUC) "fname");
	if (fname == NULL) {
	    err = 1;
	} else {
	    session_file_make_path(fullname, (const char *) fname);
	    pmod = rebuild_session_model(fullname, &err);
	}

	if (!err) {
	    name = xmlGetProp(cur, (XUC) "name");
	}

	if (name == NULL) {
	    err = 1;
	} else {
	    session.models[i] = session_model_new(pmod, (const char *) name, 
						  GRETL_OBJ_EQN);
	    if (session.models[i] == NULL) {
		err = E_ALLOC;
	    }
	}

	free(fname);
	free(name);

	cur = cur->next;
	i++;
    }

#if SESSION_DEBUG
    fprintf(stderr, "restore_session_models: returning %d\n", err);
#endif

    /* FIXME clean up on error */

    return err;
}

static int get_data_submask (xmlDocPtr doc, xmlNodePtr cur,
			     struct sample_info *sinfo)
{
    int i, len, err = 0;

    if (!gretl_xml_get_prop_as_int(cur, "length", &len)) {
	return 1;
    }

    if (!gretl_xml_get_prop_as_int(cur, "mode", &sinfo->mode)) {
	return 1;
    }    

    sinfo->mask = calloc(len, 1);
    if (sinfo->mask == NULL) {
	err = 1;
    } else {
	xmlChar *tmp = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	if (tmp == NULL) {
	    err = 1;
	} else {
	    const char *s = (const char *) tmp;
	    int si;

	    for (i=0; i<len; i++) {
		sscanf(s, "%d", &si);
		s += strspn(s, " ");
		s += strcspn(s, " ");
		if (si != 0) {
		    sinfo->mask[i] = si;
		}
	    }
	    free(tmp);
	}
    }

    return err;
}

static int 
read_session_xml (const char *fname, struct sample_info *sinfo) 
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    xmlChar *tmp;
    int err = 0;
    int to_iso_latin = 0;

    LIBXML_TEST_VERSION
	xmlKeepBlanksDefault(0);

    err = gretl_xml_open_doc_root(fname, "gretl-session", &doc, &cur);

    if (err) {
	gui_errmsg(err);
	return 1;
    }

#ifndef USE_GTK2
    if (doc->encoding != NULL && strstr((char *) doc->encoding, "UTF")) {
	to_iso_latin = 1;
    }
#endif

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
	    err = get_data_submask(doc, cur, sinfo);
	} else if (!xmlStrcmp(cur->name, (XUC) "models")) {
	    tmp = xmlGetProp(cur, (XUC) "count");
	    if (tmp != NULL) {
		session.nmodels = atoi((const char *) tmp);
		free(tmp);
		if (session.nmodels > 0) {
		    err = restore_session_models(doc, cur);
		}		
	    }
        } else if (!xmlStrcmp(cur->name, (XUC) "graphs")) {
	    tmp = xmlGetProp(cur, (XUC) "count");
	    if (tmp != NULL) {
		session.ngraphs = atoi((const char *) tmp);
		free(tmp);
		if (session.ngraphs > 0) {
		    err = restore_session_graphs(cur);
		}
	    }	    
	} else if (!xmlStrcmp(cur->name, (XUC) "texts")) {
	    tmp = xmlGetProp(cur, (XUC) "count");
	    if (tmp != NULL) {
		session.ntexts = atoi((const char *) tmp);
		free(tmp);
		if (session.ntexts > 0) {
		    err = restore_session_texts(doc, cur);
		}		
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "notes")) {
	    session.notes = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	    if (session.notes == NULL) {
		err = 1;
	    }
	}
	if (!err) {
	    cur = cur->next;
	}
    }

    if (doc != NULL) {
	xmlFreeDoc(doc);
	xmlCleanupParser();
    }

    return err;
}

static int write_session_xml (void)
{
    char fname[MAXLEN];
    char tmpname[MAXLEN];
    FILE *fp, *fq;
    int i, err = 0;

    chdir(paths.userdir);

    sprintf(fname, "%s%csession.xml", session.dirname, SLASH);
    fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	errbox("Couldn't write session file");
	return E_FOPEN;
    }

    gretl_xml_header(fp);
    fputs("<gretl-session>\n", fp);
    fprintf(fp, " <sample t1=\"%d\" t2=\"%d\"/>\n", datainfo->t1, datainfo->t2);
    write_datainfo_submask(datainfo, fp);

    fprintf(fp, " <models count=\"%d\">\n", session.nmodels);
    for (i=0; i<session.nmodels && !err; i++) {
	if (session.models[i]->type == GRETL_OBJ_EQN) {
	    sprintf(tmpname, "%s%cmodel.%d", session.dirname, SLASH, i+1);
	    fq = fopen(tmpname, "w");
	    if (fq == NULL) {
		errbox("Couldn't write session model file");
		err = E_FOPEN;
	    } else {
		sprintf(tmpname, "model.%d", i+1);
		fprintf(fp, "  <session-model name=\"%s\" fname=\"%s\" addr=\"%p\"/>\n", 
			session.models[i]->name, tmpname, session.models[i]->ptr);
		gretl_model_serialize(session.models[i]->ptr, fq);
		fclose(fq);
	    }
	} else {
	    fprintf(stderr, "FIXME models other than single-equation ones\n");
	}
    }

    if (err) {
	fclose(fp);
	remove(fname);
	return err;
    }

    fputs(" </models>\n", fp);

    fprintf(fp, " <graphs count=\"%d\">\n", session.ngraphs);
    for (i=0; i<session.ngraphs; i++) {
	fprintf(fp, "  <session-graph name=\"%s\" fname=\"%s\" "
		"ID=\"%d\" type=\"%d\"/>\n", 
		session.graphs[i]->name, session.graphs[i]->fname,
		session.graphs[i]->ID, session.graphs[i]->type);
    } 
    fputs(" </graphs>\n", fp);

    fprintf(fp, " <texts count=\"%d\">\n", session.ntexts);
    for (i=0; i<session.ntexts; i++) {
	fprintf(fp, "  <session-text name=\"%s\">", session.texts[i]->name);
	/* XML encoding? */
	fputs(session.texts[i]->buf, fp);
	fputs("</session-text>\n", fp);
    }    
    fputs(" </texts>\n", fp);

    if (session.notes != NULL) {
	fputs(" <notes>", fp);
	fputs(session.notes, fp);
	fputs("</notes>\n", fp);
    }	

    fputs("</gretl-session>\n", fp);

    fclose(fp);

    return 0;
}
