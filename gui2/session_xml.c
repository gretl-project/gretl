

static int check_graph_file (const char *fname)
{
    char fullname[MAXLEN];
    FILE *fp;
    int err = 0;

    sprintf(fullname, "%s%c%s", session.dirname, SLASH, fname);

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
    xmlChar *tmp;
    int i = 0;
    int err = 0;

    session.graphs = mymalloc(session.ngraphs * sizeof *session.graphs);
    if (session.graphs == NULL) {
	return E_ALLOC;
    }

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	
	session.graphs[i] = mymalloc(sizeof **session.graphs);
	if (session.graphs[i] == NULL) {
	    return E_ALLOC;
	}

	session.graphs[i]->name[0] = 0;
	session.graphs[i]->fname[0] = 0;

	tmp = xmlGetProp(cur, (XUC) "name");
	if (tmp != NULL) {
	    strncat(session.graphs[i]->name, (const char *) tmp, 31);
	    free(tmp);
	} else {
	    err = 1;
	}

	tmp = xmlGetProp(cur, (XUC) "fname");
	if (tmp != NULL) {
	    strncat(session.graphs[i]->fname, (const char *) tmp, MAXLEN - 1);
	    free(tmp);
	    err = check_graph_file(session.graphs[i]->fname);
	} else {
	    err = 1;
	}

	cur = cur->next;
	i++;
    }

    return err;
}

static int restore_session_texts (xmlDocPtr doc, xmlNodePtr node)
{
    xmlNodePtr cur;
    xmlChar *tmp;
    char *buf;
    int i = 0;
    int err = 0;

    session.texts = mymalloc(session.ntexts * sizeof *session.texts);
    if (session.texts == NULL) {
	return E_ALLOC;
    }

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {

	session.texts[i] = mymalloc(sizeof **session.texts);
	if (session.texts[i] == NULL) {
	    return E_ALLOC;
	}

	session.texts[i]->name[0] = 0;
	session.texts[i]->buf = NULL;	

	tmp = xmlGetProp(cur, (XUC) "name");
	if (tmp != NULL) {
	    strncat(session.texts[i]->name, (const char *) tmp, OBJNAMLEN - 1);
	    free(tmp);
	    buf = (char *) xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
	    if (buf != NULL) {
		session.texts[i]->buf = buf;
	    } else {
		err = 1;
	    }
	} else {
	    err = 1;
	}

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

    *err = gretl_xml_open_doc_root(fname, "gretl-model", &doc, &node);
    if (*err) {
	return NULL;
    }

    pmod = gretl_model_new();
    if (pmod == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

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

    if (got < 11) {
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
	}
	cur = cur->next;
    }

    /* FIXME some other elements not accounted for */

    gretl_pop_c_numeric_locale();

 bailout:

    xmlFreeDoc(doc);
    
    /* need to clean up on error here */

    return pmod;
}

static int restore_session_models (xmlDocPtr doc, xmlNodePtr node)
{
    char fname[MAXLEN];
    xmlNodePtr cur;
    xmlChar *tmp;
    int i = 0;
    int err = 0;

    session.models = mymalloc(session.nmodels * sizeof *session.models);
    if (session.models == NULL) {
	return E_ALLOC;
    }

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	MODEL *pmod = NULL;

	tmp = xmlGetProp(cur, (XUC) "fname");
	if (tmp == NULL) {
	    err = E_DATA;
	} else {
	    sprintf(fname, "%s%c%s", session.dirname, SLASH, 
		    (const char *) tmp);
	    pmod = rebuild_session_model(fname, &err);
	    free(tmp);
	}

	if (err) {
	    break;
	}

	session.models[i] = mymalloc(sizeof **session.models);
	if (session.models[i] == NULL) {
	    return E_ALLOC;
	}

	session.models[i]->name[0] = 0;
	session.models[i]->type = GRETL_OBJ_EQN;
	session.models[i]->ptr = pmod;	

	tmp = xmlGetProp(cur, (XUC) "name");
	if (tmp != NULL) {
	    strncat(session.models[i]->name, (const char *) tmp, OBJNAMLEN - 1);
	    free(tmp);
	} else {
	    err = 1;
	}

	cur = cur->next;
	i++;
    }

    /* FIXME clean up on error */

    return err;
}

static int read_session_xml (const char *fname, int *t1, int *t2) 
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    xmlChar *tmp;
    int err = 0;
    int to_iso_latin = 0;

    /* COMPAT: Do not generate nodes for formatting spaces */
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
		*t1 = atoi((const char *) tmp);
		free(tmp);
	    } else {
		err = 1;
	    }
	    tmp = xmlGetProp(cur, (XUC) "t2");
	    if (tmp != NULL) {
		*t2 = atoi((const char *) tmp);
		free(tmp);
	    } else {
		err = 1;
	    }
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
