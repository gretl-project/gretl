#ifdef USE_GNOME

static GdkPixbuf *png_mono_pixbuf (const char *fname)
{
    FILE *fsrc, *ftmp;
    char cmd[MAXLEN], temp[MAXLEN], fline[MAXLEN];
    GdkPixbuf *pbuf = NULL;

    sprintf(temp, "%sgpttmp.XXXXXX", paths.userdir);
    if (mktemp(temp) == NULL) {
	return NULL;
    }

    ftmp = gretl_fopen(temp, "w");
    if (ftmp == NULL) {
	return NULL;
    }

    fsrc = gretl_fopen(fname, "r");
    if (fsrc == NULL) {
	fclose(ftmp);
	remove(temp);
	return NULL;
    }

    fprintf(ftmp, "set term pbm mono\n"
	    "set output '%s%s'\n", 
	    paths.userdir, GRETL_PBM_TMP);

    while (fgets(fline, MAXLEN-1, fsrc)) {
	if (strncmp(fline, "set term", 8) && 
	    strncmp(fline, "set output", 10)) {
	    fputs(fline, ftmp);
	}
    }

    fclose(fsrc);
    fclose(ftmp);

    /* run gnuplot on the temp plotfile */
    sprintf(cmd, "\"%s\" \"%s\"", paths.gnuplot, temp);
    if (system(cmd)) {
	remove(temp);
	return NULL;
    }

    remove(temp);

    build_path(paths.userdir, GRETL_PBM_TMP, temp, NULL);
#if GTK_MAJOR_VERSION >= 2
    pbuf = gdk_pixbuf_new_from_file(temp, NULL);
#else
    pbuf = gdk_pixbuf_new_from_file(temp);
#endif
    remove(temp);

    return pbuf;
}

#endif /* USE_GNOME */

/* row format specifications for RTF "tables" */

#define STATS_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx2700\\cellx4000\\cellx6700\\cellx8000\n\\intbl"

/* ............................................................. */

static void printfrtf (double zz, PRN *prn, int endrow)
{
    /* was using "qr", for right alignment */

    if (na(zz)) {
	if (endrow) {
	    pprintf(prn, "\\qc %s\\cell\\intbl \\row\n",
		    I_("undefined"));
	} else {
	    pprintf(prn, "\\qc %s\\cell", I_("undefined"));
	}
	return;
    }

    if (endrow) {
	pprintf(prn, "\\qc %#.*g\\cell\\intbl \\row\n", GRETL_DIGITS, zz);
    } else {
	pprintf(prn, "\\qc %#.*g\\cell", GRETL_DIGITS, zz);
    }
}

#define SUMM_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                   "\\cellx1600\\cellx3200\\cellx4800\\cellx6400" \
                   "\\cellx8000\n"

#define VAR_SUMM_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                      "\\cellx2000\\cellx4000\\cellx6000\\cellx8000\n"

/* ............................................................. */

void rtfprint_summary (Summary *summ,
		       const DATAINFO *pdinfo,
		       PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN], tmp[128];
    int i, vi;

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    sprintf(tmp, I_("Summary Statistics, using the observations %s - %s"),
	    date1, date2);

    pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n", tmp);
    
    if (summ->list[0] == 1) {
	sprintf(tmp, I_("for the variable %s (%d valid observations)"), 
		pdinfo->varname[summ->list[1]], summ->n);
	pprintf(prn, "%s\\par\n\n", tmp);
	pputs(prn, "{" VAR_SUMM_ROW "\\intbl ");
    } else {
	strcpy(tmp, I_("(missing values will be skipped)"));
	pprintf(prn, "%s\\par\n\n", tmp);
	pprintf(prn, "{" SUMM_ROW
		"\\intbl \\qc %s\\cell", I_("Variable"));
    }

    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    I_("Mean"), I_("Median"), I_("Minimum"), I_("Maximum"));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[vi]);
	}
	printfrtf(summ->mean[i], prn, 0);
	printfrtf(summ->median[i], prn, 0);
	printfrtf(summ->low[i], prn, 0);
	printfrtf(summ->high[i], prn, 1);
    }

    if (summ->list[0] > 1) pprintf(prn, "\\intbl \\qc %s\\cell",
				   I_("Variable"));

    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n",
	    I_("Std. Dev."), I_("C.V."), I_("Skewness"), I_("Ex. kurtosis"));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pdinfo->varname[vi]);
	}
	printfrtf(summ->sd[i], prn, 0);
	printfrtf(summ->cv[i], prn, 0);
	printfrtf(summ->skew[i], prn, 0);
	printfrtf(summ->xkurt[i], prn, 1);
    }

    pputs(prn, "}}\n");
}

/* ............................................................. */

static void printftex (double zz, PRN *prn, int endrow)
{
    if (na(zz)) {
	if (endrow) {
	    pprintf(prn, "\\multicolumn{1}{c}{%s}\\\\", I_("undefined"));
	} else {
	    pprintf(prn, "\\multicolumn{1}{c}{%s} & ", I_("undefined"));
	}
    } else {
	char s[32];

	tex_dcolumn_double(zz, s);

	if (endrow) {
	    pprintf(prn, "%s\\\\", s);
	} else {
	    pprintf(prn, "%s & ", s);
	}
    }	
}

/* ............................................................. */

void texprint_summary (Summary *summ,
		       const DATAINFO *pdinfo,
		       PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN], vname[16], tmp[128];
    int i, vi;
    char pt = get_local_decpoint();

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    sprintf(tmp, I_("Summary Statistics, using the observations %s--%s"),
	    date1, date2);

    pprintf(prn, "\\begin{center}\n%s\\\\\n", tmp);
    
    if (summ->list[0] == 1) {
	tex_escape(vname, pdinfo->varname[summ->list[1]]);
	sprintf(tmp, I_("for the variable %s (%d valid observations)"), 
		vname, summ->n);
	pprintf(prn, "%s\\\\[8pt]\n\n", tmp);
	pputs(prn, "\\begin{tabular}{rrrr}\n");
    } else {
	strcpy(tmp, I_("(missing values denoted by $-999$ will be "
		"skipped)"));
	pprintf(prn, "%s\\\\[8pt]\n\n", tmp);
	pprintf(prn, "\\begin{tabular}{lD{%c}{%c}{-1}"
		"D{%c}{%c}{-1}D{%c}{%c}{-1}D{%c}{%c}{-1}}\n", 
		pt, pt, pt, pt, pt, pt, pt, pt);
	pprintf(prn, "%s &", I_("Variable"));
    }

    pprintf(prn, " \\multicolumn{1}{c}{%s}%%\n"
	    " & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{1}{c}{%s}%%\n"
	    "   & \\multicolumn{1}{c}{%s} \\\\[1ex]\n",
	    I_("Mean"), I_("Median"), I_("Minimum"), I_("Maximum"));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    tex_escape(vname, pdinfo->varname[vi]);
	    pprintf(prn, "%s & ", vname);
	}
	printftex(summ->mean[i], prn, 0);
	printftex(summ->median[i], prn, 0);
	printftex(summ->low[i], prn, 0);
	printftex(summ->high[i], prn, 1);
	if (i == summ->list[0] - 1) {
	    pputs(prn, "[10pt]\n\n");
	} else {
	    pputc(prn, '\n');
	}
    }

    if (summ->list[0] > 1) {
	pprintf(prn, "%s & ", I_("Variable"));
    }

    pprintf(prn, " \\multicolumn{1}{c}{%s}%%\n"
	    " & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{1}{c}{%s}%%\n"
	    "   & \\multicolumn{1}{c}{%s} \\\\[1ex]\n",
	    I_("Std.\\ Dev."), I_("C.V."), I_("Skewness"), I_("Ex.\\ kurtosis"));

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (summ->list[0] > 1) {
	    tex_escape(vname, pdinfo->varname[vi]);
	    pprintf(prn, "%s & ", vname);
	}
	printftex(summ->sd[i], prn, 0);
	printftex(summ->cv[i], prn, 0);
	printftex(summ->skew[i], prn, 0);
	printftex(summ->xkurt[i], prn, 1);
	pputc(prn, '\n');
    }

    pputs(prn, "\\end{tabular}\n\\end{center}\n");
    
}

/* ......................................................... */ 

static void tex_outxx (double xx, PRN *prn)
{
    if (na(xx)) {
	pprintf(prn, "%s & ", I_("undefined"));
    } else {
	pprintf(prn, "$%.4f$ & ", xx);
    }
}

/* ......................................................... */ 

static void rtf_outxx (double xx, PRN *prn)
{
    if (na(xx)) {
	pprintf(prn, "\\qc %s\\cell ", I_("undefined"));
    } else {
	pprintf(prn, "\\qc %.4f\\cell ", xx);	
    }
}

/* ......................................................... */ 

static void rtf_corr_row (int lo, PRN *prn)
{
    pputs(prn, "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262");

    if (lo == 2) {
	pputs(prn, "\\cellx1500\\cellx3000\\cellx3500\n");
    }
    else if (lo == 3) {
	pputs(prn, "\\cellx1500\\cellx3000\\cellx4500\\cellx5000\n");
    }
    else if (lo == 4) {
	pputs(prn, "\\cellx1500\\cellx3000\\cellx4500\\cellx6000"
		"\\cellx6500\n");
    }
    else {
	pputs(prn, "\\cellx1500\\cellx3000\\cellx4500\\cellx6000"
		"\\cellx7500\\cellx8000\n");
    }

    pputs(prn, "\\intbl ");
}

/* ......................................................... */ 

static void rtf_table_pad (int pad, PRN *prn)
{
    while (pad--) pputs(prn, "\\cell ");
}

/* ......................................................... */ 

static void
rtfprint_matrix (const double *vec, const int *list,
		 int t1, int t2, int n, int ci,
		 const DATAINFO *pdinfo, PRN *prn)
{
    register int i, j;
    int lo, ljnf, nf, li2, p, k, index, ij2;
    char tmp[128];
    enum { FIELDS = 5 };

    if (ci == CORR) {
	char date1[OBSLEN], date2[OBSLEN];

	ntodate(date1, t1, pdinfo);
	ntodate(date2, t2, pdinfo);

	sprintf(tmp, I_("Correlation coefficients, using the observations "
			"%s - %s"), date1, date2);
	pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n(%s)\\par\n",
		tmp, I_("skipping any missing values"));

	sprintf(tmp, I_("5%% critical value (two-tailed) = %.4f for n = %d"), 
		rhocrit95(n), n);
	pprintf(prn, "%s\\par\n\\par\n{", tmp);
    } 
    else if (ci == COVAR) {
	pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n\\par\n{",
		I_("Coefficient covariance matrix"));
    }
    
    lo = list[0];

    for (i=0; i<=lo/FIELDS; i++) {
	int pad;

	nf = i * FIELDS;
	li2 = lo - nf;
	p = (li2 > FIELDS) ? FIELDS : li2;
	if (p == 0) break;

	pad = (lo > FIELDS)? FIELDS - p : lo - p;

	rtf_corr_row(lo, prn);

	if (pad) rtf_table_pad(pad, prn);

	/* print the varname headings */
	for (j=1; j<=p; ++j)  {
	    ljnf = list[j + nf];
	    pprintf(prn, "%d) %s\\cell %s", ljnf, pdinfo->varname[ljnf],
		    (j == p)? "\\cell \\intbl \\row\n" : "");
	}

	/* print rectangular part, if any, of matrix */
	for (j=0; j<nf; j++) {
	    pputs(prn, "\\intbl "); 
	    if (pad) rtf_table_pad(pad, prn);
	    for (k=0; k<p; k++) {
		index = ijton(j, nf+k, lo);
		if (ci == CORR) {
		    rtf_outxx(vec[index], prn);
		} else {
		    printfrtf(vec[index], prn, 0);
		}
	    }
	    pprintf(prn, "\\ql (%d\\cell \\intbl \\row\n", list[j+1]);
	}

	/* print upper triangular part of matrix */
	for (j=0; j<p; ++j) {
	    pputs(prn, "\\intbl "); 
	    rtf_table_pad(pad + j, prn);
	    ij2 = nf + j;
	    for (k=j; k<p; k++) {
		index = ijton(ij2, nf+k, lo);
		if (ci == CORR) {
		    rtf_outxx(vec[index], prn);
		} else {
		    printfrtf(vec[index], prn, 0);
		}
	    }
	    pprintf(prn, "\\ql (%d\\cell \\intbl \\row\n", list[ij2+1]);
	}
    }
    pputs(prn, "}}\n");
}

/* ........................................................... */

void rtfprint_corrmat (CorrMat *corr,
		       const DATAINFO *pdinfo, 
		       PRN *prn)
{
    rtfprint_matrix(corr->xpx, corr->list, corr->t1, corr->t2,
		    corr->n, CORR, pdinfo, prn);
}

/* ......................................................... */

static void
texprint_matrix (const double *vec, const int *list,
		 int t1, int t2, int n, int ci,
		 const DATAINFO *pdinfo, PRN *prn)
{
    register int i, j;
    int lo, ljnf, nf, li2, p, k, index, ij2;
    char vname[16], tmp[128];
    int fields;

    if (ci == CORR) fields = 5;
    else fields = 4;

    lo = list[0];

    if (ci == CORR) {
	char date1[OBSLEN], date2[OBSLEN];

	ntodate(date1, t1, pdinfo);
	ntodate(date2, t2, pdinfo);

	sprintf(tmp, I_("Correlation coefficients, using the observations "
			"%s--%s"), date1, date2);
	pprintf(prn, "\\begin{center}\n%s\\\\\n(%s)\\\\\n", 
		tmp, I_("skipping any missing values"));

	sprintf(tmp, I_("5\\%% critical value (two-tailed) = %.4f for n = %d"), 
		rhocrit95(n), n);
	pprintf(prn, "%s\\\\\n", tmp);
    }
    else if (ci == COVAR) {
	pprintf(prn, "\\begin{center}\n%s\\\\\n", 
		I_("Coefficient covariance matrix"));
    }

    pputs(prn, "\\vspace{8pt}\n");

    if (ci == CORR) {
	pprintf(prn, "\\begin{tabular}{rrr%s}\n",
		(lo == 3)? "r" : (lo == 4)? "rr" : "rrr");
    } else {
	char pt = get_local_decpoint();

	pputs(prn, "\\begin{tabular}{");
	for (i=0; i<=lo && i<fields; i++) {
	    pprintf(prn, "D{%c}{%c}{-1}", pt, pt);
	}
	pputs(prn, "r}\n");
    }

    for (i=0; i<=lo/fields; i++) {
	nf = i * fields;
	li2 = lo - nf;
	p = (li2 > fields) ? fields : li2;
	if (p == 0) break;

	/* print the varname headings */
	for (j=1; j<=p; ++j)  {
	    ljnf = list[j + nf];
	    tex_escape(vname, pdinfo->varname[ljnf]);
	    if (ci == CORR) {
		pprintf(prn, "%d) %s%s", ljnf, vname,
			(j == p)? " &\\\\" : " & ");
	    } else {
		pprintf(prn, "\\multicolumn{1}{c}{%d) %s}%s", ljnf, vname,
			(j == p)? " &\\\\\n" : " &\n");
	    }
	}
	
	/* insert spacers */
	if (ci == CORR) {
	    for (j=1; j<=p; ++j) {
		pputs(prn, "\\rule{13ex}{0pt} & ");
	    }
	    pputs(prn, "\\\\\[-6pt]\n"); 
	}   

	/* print rectangular part, if any, of matrix */
	for (j=0; j<nf; j++) {
	    for (k=0; k<p; k++) {
		index = ijton(j, nf+k, lo);
		if (ci == CORR) {
		    tex_outxx(vec[index], prn);
		} else {
		    printftex(vec[index], prn, 0);
		}
	    }
	    pprintf(prn, "(%d\\\\\n", list[j+1]);
	}

	/* print upper triangular part of matrix */
	for (j=0; j<p; ++j) {
	    ij2 = nf + j;
	    for (k=0; k<j; k++) pputs(prn, " & ");
	    for (k=j; k<p; k++) {
		index = ijton(ij2, nf+k, lo);
		if (ci == CORR) {
		    tex_outxx(vec[index], prn);
		} else {
		    printftex(vec[index], prn, 0);
		}
	    }
	    pprintf(prn, "(%d\\\\\n", list[ij2+1]);
	}
	pputs(prn, "\\\\\n");
    }
    pputs(prn, "\\end{tabular}\n\\end{center}\n");
}

/* ........................................................... */

void texprint_corrmat (CorrMat *corr,
		       const DATAINFO *pdinfo, 
		       PRN *prn)
{
    texprint_matrix(corr->xpx, corr->list, corr->t1, corr->t2,
		    corr->n, CORR, pdinfo, prn);
}

/* ........................................................... */

static 
void tex_fit_resid_head (const FITRESID *fr, const DATAINFO *pdinfo, 
			 PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN]; 

    ntodate(date1, fr->t1, pdinfo);
    ntodate(date2, fr->t2, pdinfo);

    pputs(prn, "\\begin{raggedright}\n");
    pputs(prn, I_("Model estimation range:"));
    pprintf(prn, " %s--%s", date1, date2);

    pprintf(prn, " ($n$ = %d)\\\\\n", fr->real_nobs); 

    pprintf(prn, I_("Standard error of residuals = %g"), fr->sigma);
    pputs(prn, "\n\\end{raggedright}\n");
}

/* ........................................................... */

static 
void rtf_fit_resid_head (const FITRESID *fr, const DATAINFO *pdinfo, 
			 PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN]; 
    char tmp[128];

    ntodate(date1, fr->t1, pdinfo);
    ntodate(date2, fr->t2, pdinfo);

    pputs(prn, "{\\rtf1\\par\n\\qc ");
    pputs(prn, I_("Model estimation range:")); 
    pprintf(prn, " %s - %s (n = %d)\\par\n", date1, date2, fr->real_nobs);

    sprintf(tmp, I_("Standard error of residuals = %g"), 
	    fr->sigma);
    pprintf(prn, "\\qc %s\\par\n\\par\n", tmp);
}

/* ........................................................... */

void 
texprint_fit_resid (const FITRESID *fr, const DATAINFO *pdinfo, PRN *prn)
{
    int t, anyast = 0;
    double xx;
    char vname[16];

    tex_fit_resid_head(fr, pdinfo, prn); 

    tex_escape(vname, fr->depvar);

    pprintf(prn, "\n\\begin{center}\n"
	    "\\begin{tabular}{rrrrl}\n"
	    " & \n"
	    " \\multicolumn{1}{c}{%s} & \n"
	    "  \\multicolumn{1}{c}{%s} & \n"
	    "   \\multicolumn{1}{c}{%s}\\\\\n",
	    vname, I_("fitted"), I_("residuals"));

    for (t=0; t<fr->nobs; t++) {
	print_obs_marker(t + fr->t1, pdinfo, prn);
	pputs(prn, " & ");

	if (na(fr->actual[t])) {
	    pputs(prn, "\\\\\n");
	} else if (na(fr->fitted[t])) {
	    pprintf(prn, "%13.*f \\\\\n", fr->pmax, fr->actual[t]);
	} else {
	    int ast;

	    xx = fr->actual[t] - fr->fitted[t];
	    ast = (fabs(xx) > 2.5 * fr->sigma);
	    if (ast) anyast = 1;
	    if (fr->pmax != PMAX_NOT_AVAILABLE) {
		pprintf(prn, "%13.*f & %13.*f & %13.*f & %s \\\\\n", 
			fr->pmax, fr->actual[t],
			fr->pmax, fr->fitted[t], fr->pmax, xx,
			(ast)? " *" : "");
	    } else {
		pprintf(prn, "%13g & %13g & %13g & %s \\\\\n", 
			fr->actual[t],
			fr->fitted[t], xx,
			(ast)? " *" : "");
	    }
	}
    }

    pputs(prn, "\\end{tabular}\n"
	  "\\end{center}\n\n");

    if (anyast) pputs(prn, I_("\\textit{Note}: * denotes a residual "
			      "in excess of 2.5 standard errors\n\n"));
}

/* .................................................................. */

#define FR_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                "\\cellx800\\cellx2400\\cellx4000\\cellx5600" \
                "\\cellx6100\n"

void rtfprint_fit_resid (const FITRESID *fr, 
			 const DATAINFO *pdinfo, 
			 PRN *prn)
{
    double xx;
    int anyast = 0;
    int t;

    rtf_fit_resid_head(fr, pdinfo, prn);

    pputs(prn, "{" FR_ROW "\\intbl ");
    pprintf(prn, 
	    " \\qc \\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\ql \\cell"
	    " \\intbl \\row\n",
	    fr->depvar, I_("fitted"), I_("residual"));

    for (t=0; t<fr->nobs; t++) {
	pputs(prn, "\\qr ");
	print_obs_marker(t + fr->t1, pdinfo, prn);
	pputs(prn, "\\cell"); 

	if (na(fr->actual[t])) {
	    pputs(prn, "\\qc \\cell \\qc \\cell \\qc \\cell \\ql \\cell"
		  " \\intbl \\row\n"); 
	} else if (na(fr->fitted[t])) {	 
	    printfrtf(fr->actual[t], prn, 0);
	    pputs(prn, "\\qc \\cell \\qc \\cell \\ql \\cell"
		  " \\intbl \\row\n"); 
	} else {
	    int ast;

	    xx = fr->actual[t] - fr->fitted[t];
	    ast = (fabs(xx) > 2.5 * fr->sigma);
	    if (ast) {
		anyast = 1;
	    }
	    printfrtf(fr->actual[t], prn, 0);
	    printfrtf(fr->fitted[t], prn, 0);
	    printfrtf(xx, prn, 0);
	    pprintf(prn, "\\ql %s\\cell \\intbl \\row\n", 
		    (ast)? "*" : "");
	}
    }

    pputs(prn, "}\n");
    if (anyast) {
	pprintf(prn, "\\par\n\\qc %s \\par\n",
		I_("Note: * denotes a residual in excess of 2.5 standard errors"));
    }
    pputs(prn, "}\n");
}

/* .................................................................. */

static void texprint_fcast_without_errs (const FITRESID *fr, 
					 const DATAINFO *pdinfo, 
					 PRN *prn)
{
    char actual[32], fitted[32];
    char vname[16];
    char pt = get_local_decpoint();
    int t;

    pputs(prn, "%% The table below needs the \"dcolumn\" package\n\n");

    pprintf(prn, "\\begin{center}\n"
	    "\\begin{tabular}{%%\n"
	    "r%% col 1: obs\n"
	    "  l%% col 2: varname\n"
	    "    D{%c}{%c}{-1}%% col 3: fitted\n",
	    pt, pt);

    tex_escape(vname, fr->depvar);

    pprintf(prn, "%s & %s & \\multicolumn{1}{c}{%s} \\\\\n",
	    I_("Obs"), vname, I_("prediction"));

    for (t=0; t<fr->nobs; t++) {
	tex_dcolumn_double(fr->actual[t], actual);
	tex_dcolumn_double(fr->fitted[t], fitted);
	print_obs_marker(t + fr->t1, pdinfo, prn);
	pprintf(prn, " & %s & %s \\\\\n",
		actual, fitted);
    }

    pputs(prn, "\\end{tabular}\n"
	  "\\end{center}\n\n");
}

static void real_texprint_fcast_with_errs (const FITRESID *fr, 
					   const DATAINFO *pdinfo, 
					   PRN *prn)
{
    int t;
    double maxerr;
    char actual[32], fitted[32], sderr[32], lo[32], hi[32];
    char vname[16];
    char pt = get_local_decpoint();

    pprintf(prn, I_("For 95 percent confidence intervals, "
		    "$t(%d, .025) = %.3f$\n\n"), 
	    fr->df, fr->tval);

    pputs(prn, "%% The table below needs the \"dcolumn\" package\n\n");

    pprintf(prn, "\\begin{center}\n"
	    "\\begin{tabular}{%%\n"
	    "r%% col 1: obs\n"
	    "  l%% col 2: varname\n"
	    "    D{%c}{%c}{-1}%% col 3: fitted\n"
	    "      D{%c}{%c}{-1}%% col 4: std error\n"
	    "        D{%c}{%c}{-1}%% col 5: conf int lo\n"
	    "         D{%c}{%c}{-1}}%% col 5: conf int hi\n",
	    pt, pt, pt, pt, pt, pt, pt, pt);

    tex_escape(vname, fr->depvar);

    pprintf(prn, "%s & %s & \\multicolumn{1}{c}{%s}\n"
	    " & \\multicolumn{1}{c}{%s}\n"
	    "  & \\multicolumn{2}{c}{%s} \\\\\n",
	    I_("Obs"), vname,
	    I_("prediction"), I_("std. error"),
	    /* xgettext:no-c-format */
	    I_("95\\% confidence interval"));

    pputs(prn, "& & & & \\multicolumn{1}{c}{low} & "
	  "\\multicolumn{1}{c}{high} \\\\\n");

    for (t=0; t<fr->nobs; t++) {
	maxerr = fr->tval * fr->sderr[t];
	tex_dcolumn_double(fr->actual[t], actual);
	tex_dcolumn_double(fr->fitted[t], fitted);
	tex_dcolumn_double(fr->sderr[t], sderr);
	tex_dcolumn_double(fr->fitted[t] - maxerr, lo);
	tex_dcolumn_double(fr->fitted[t] + maxerr, hi);
	print_obs_marker(t + fr->t1, pdinfo, prn);
	pprintf(prn, " & %s & %s & %s & %s & %s \\\\\n",
		actual, fitted, sderr, lo, hi);
    }

    pputs(prn, "\\end{tabular}\n"
	  "\\end{center}\n\n");
}

void texprint_fcast_with_errs (const FITRESID *fr, 
			       const DATAINFO *pdinfo, 
			       PRN *prn)
{
    if (fr->sderr != NULL) {
	real_texprint_fcast_with_errs(fr, pdinfo, prn);
    } else {
	texprint_fcast_without_errs(fr, pdinfo, prn);
    }
}


/* .................................................................. */

#define FC_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                "\\cellx800\\cellx2200\\cellx3600\n"

#define FCE_ROW  "\\trowd \\trqc \\trgaph60\\trleft-30\\trrh262" \
                 "\\cellx800\\cellx2200\\cellx3600\\cellx5000" \
                 "\\cellx7800\n"

static void rtfprint_fcast_without_errs (const FITRESID *fr, 
					 const DATAINFO *pdinfo, 
					 PRN *prn)
{
    int t;

    pputs(prn, "{\\rtf1\\par\n\n");

    pputs(prn, "{" FC_ROW "\\intbl ");

    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n", 
	    I_("Obs"), fr->depvar, I_("prediction")); 

    for (t=0; t<fr->nobs; t++) {
	pputs(prn, "\\qr ");
	print_obs_marker(t + fr->t1, pdinfo, prn);
	pputs(prn, "\\cell"); 
	printfrtf(fr->actual[t], prn, 0);
	printfrtf(fr->fitted[t], prn, 0);
    }

    pputs(prn, "}}\n");
}

static void real_rtfprint_fcast_with_errs (const FITRESID *fr, 
					   const DATAINFO *pdinfo, 
					   PRN *prn)
{
    int t;
    double maxerr;
    char tmp[128];

    sprintf(tmp, I_("For 95 percent confidence intervals, "
		    "t(%d, .025) = %.3f"), 
	    fr->df, fr->tval);

    pprintf(prn, "{\\rtf1\\par\n\\qc %s\\par\n\\par\n", tmp);

    pputs(prn, "{" FCE_ROW "\\intbl ");
    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n", 
	    I_("Obs"), fr->depvar, I_("prediction"), 
	    I_("std. error"),
	    /* xgettext:no-c-format */
	    I_("95% confidence interval"));

    for (t=0; t<fr->nobs; t++) {
	pputs(prn, "\\qr ");
	print_obs_marker(t + fr->t1, pdinfo, prn);
	pputs(prn, "\\cell"); 
	maxerr = fr->tval * fr->sderr[t];
	printfrtf(fr->actual[t], prn, 0);
	printfrtf(fr->fitted[t], prn, 0);
	printfrtf(fr->sderr[t], prn, 0);
	pprintf(prn, "\\qc (%#.*g, %#.*g)\\cell \\intbl \\row\n", 
		GRETL_DIGITS, fr->fitted[t] - maxerr, 
		GRETL_DIGITS, fr->fitted[t] + maxerr);
    }

    pputs(prn, "}}\n");
}

void rtfprint_fcast_with_errs (const FITRESID *fr, 
			       const DATAINFO *pdinfo, 
			       PRN *prn)
{
    if (fr->sderr != NULL) {
	real_rtfprint_fcast_with_errs(fr, pdinfo, prn);
    } else {
	rtfprint_fcast_without_errs(fr, pdinfo, prn);
    }
}

/* .................................................................. */

static void 
texprint_coeff_interval (const CONFINT *cf, int i, PRN *prn)
{
    char vname[16];

    tex_escape(vname, cf->names[i]);
    pprintf(prn, " %s & ", vname);

    if (isnan(cf->coeff[i])) {
	pprintf(prn, "\\multicolumn{1}{c}{%s} & ", I_("undefined"));
    } else {
	char coeff[32];

	tex_dcolumn_double(cf->coeff[i], coeff);
	pprintf(prn, "%s & ", coeff);
    }

    if (isnan(cf->maxerr[i])) {
	pprintf(prn, "\\multicolumn{2}{c}{%s}", I_("undefined"));
    } else {
	char lo[32], hi[32];

	tex_dcolumn_double(cf->coeff[i] - cf->maxerr[i], lo);
	tex_dcolumn_double(cf->coeff[i] + cf->maxerr[i], hi);
	pprintf(prn, "%s & %s", lo, hi);
    }
    pputs(prn, "\\\\\n");
}

/* .................................................................. */

void texprint_confints (const CONFINT *cf, PRN *prn)
{
    char pt = get_local_decpoint();
    int i;

    pprintf(prn, "$t(%d, .025) = %.3f$\n\n", cf->df, tcrit95(cf->df));

    pputs(prn, "%% The table below needs the \"dcolumn\" package\n\n");

    pprintf(prn, "\\begin{center}\n"
	    "\\begin{tabular}{rD{%c}{%c}{-1}D{%c}{%c}{-1}D{%c}{%c}{-1}}\n",
	    pt, pt, pt, pt, pt, pt);

    pprintf(prn, " %s%%\n"
	    " & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{2}{c}{%s}\\\\\n",
	    I_("Variable"), I_("Coefficient"),
	    /* xgettext:no-c-format */
	    I_("95\\% confidence interval"));

    pprintf(prn, " & & \\multicolumn{1}{c}{%s}%%\n"
	    "  & \\multicolumn{1}{c}{%s}\\\\\n",
	    I_("low"), I_("high"));

    for (i=0; i<cf->ncoeff; i++) {
	texprint_coeff_interval(cf, i, prn);
    }

    pputs(prn, "\\end{tabular}\n"
	  "\\end{center}\n");
}

/* .................................................................. */

static void 
rtfprint_coeff_interval (const CONFINT *cf, int i, PRN *prn)
{
    pprintf(prn, "\\qc %s\\cell", cf->names[i]);

    printfrtf(cf->coeff[i], prn, 0);

    if (isnan(cf->maxerr[i])) {
	pprintf(prn, "\\qc %s\\cell ", I_("undefined"));
    } else {
	pprintf(prn, "\\qc (%#.*g, %#.*g)\\cell ", 
		GRETL_DIGITS, cf->coeff[i] - cf->maxerr[i], 
		GRETL_DIGITS, cf->coeff[i] + cf->maxerr[i]);
    }
    pputs(prn, " \\intbl \\row\n");
}

#define CF_ROW  "\\trowd \\trgaph60\\trleft-30\\trrh262" \
                "\\cellx2400\\cellx4000\\cellx7200\n" 

void rtfprint_confints (const CONFINT *cf, PRN *prn)
{
    int i;

    pprintf(prn, "{\\rtf1\\par\n\\qc t(%d, .025) = %.3f\\par\n\\par\n", 
	    cf->df, tcrit95(cf->df));

    pputs(prn, "{" CF_ROW "\\intbl ");
    pprintf(prn, 
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\qc %s\\cell"
	    " \\intbl \\row\n", 
	    I_("Variable"), I_("Coefficient"), 
	    /* xgettext:no-c-format */
	    I_("95% confidence interval"));

    for (i=0; i<cf->ncoeff; i++) {
	rtfprint_coeff_interval(cf, i, prn);
    }

    pputs(prn, "}}\n");
}

/* .................................................................. */

void texprint_vcv (const VCV *vcv, 
                   const DATAINFO *pdinfo, 
                   PRN *prn)
{
    texprint_matrix(vcv->vec, vcv->list, 0, 0,
                    0, COVAR, pdinfo, prn);
}

/* .................................................................. */

void rtfprint_vcv (const VCV *vcv,
                   const DATAINFO *pdinfo, 
                   PRN *prn)
{
    rtfprint_matrix(vcv->vec, vcv->list, 0, 0,
                    0, COVAR, pdinfo, prn);
}

/* copy data to buffer in CSV format and place on clipboard */

#define SCALAR_DIGITS 12

static int data_to_buf_as_csv (const int *list, PRN *prn)
{
    int i, t, l0 = list[0];
    int *pmax = NULL;
    int tsamp = datainfo->t2 - datainfo->t1 + 1;
    char delim = datainfo->delim;
    double xx;
    char tmp[OBSLEN];

    if (l0 == 0) return 1;

    if (delim == ',' && ',' == datainfo->decpoint) {
	errbox(_("You can't use the same character for "
		 "the column delimiter and the decimal point"));
	return 1;
    }

    pmax = malloc(l0 * sizeof *pmax);
    if (pmax == NULL) return 1;

    for (i=1; i<=l0; i++) {
	if (datainfo->vector[list[i]]) {
	    pmax[i-1] = get_precision(&Z[list[i]][datainfo->t1], 
				      tsamp, 8);
	} else {
	    pmax[i-1] = SCALAR_DIGITS;
	}
    }	

#ifdef ENABLE_NLS
    if (datainfo->decpoint != ',') {
	setlocale(LC_NUMERIC, "C");
    }
#endif

    /* variable names */
    pprintf(prn, "obs%c", delim);
    for (i=1; i<l0; i++) {
	pprintf(prn, "%s%c", datainfo->varname[list[i]], delim);
    }
    pprintf(prn, "%s\n", datainfo->varname[list[l0]]);

    /* actual data values */
    for (t=datainfo->t1; t<=datainfo->t2; t++) {
	if (datainfo->S != NULL) {
	    pprintf(prn, "%s%c", datainfo->S[t], delim);
	} else {
	    ntodate(tmp, t, datainfo);
	    /* Does the "'" work correctly below? */
	    pprintf(prn, "\"'%s\"%c", tmp, delim); 
	}
	for (i=1; i<=l0; i++) { 
	    xx = (datainfo->vector[list[i]])? 
		Z[list[i]][t] : Z[list[i]][0];
	    if (na(xx)) {
		pputs(prn, "NA");
	    } else if (pmax[i-1] == PMAX_NOT_AVAILABLE) {
		pprintf(prn, "%.10g", xx);
	    } else {
		pprintf(prn, "%.*f", pmax[i-1], xx);
	    }
	    pprintf(prn, "%c", (i < l0)? delim : '\n');
	}
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    if (pmax) free(pmax);

    return 0;
}

static int real_csv_to_clipboard (const char *liststr)
{
    char line[MAXLINE];
    int *list = NULL;
    PRN *prn = NULL;
    int err = 0;

    sprintf(line, "store csv %s", liststr);
    list = command_list_from_string(line);

    if (list != NULL) {
	err = bufopen(&prn);
	if (!err) {
	    err = data_to_buf_as_csv(list, prn);
	}
	if (!err) {
	    err = prn_to_clipboard(prn, COPY_CSV);
	}
    }

    free(list);
    gretl_print_destroy(prn);

    return err;
}

int csv_to_clipboard (void)
{
    int err = 0;

    delimiter_dialog();
    data_save_selection_wrapper(COPY_CSV);

    if (storelist != NULL && *storelist != 0) {
	err = real_csv_to_clipboard(storelist);
	free(storelist);
	storelist = NULL;
    }

    return err;
}

int csv_selected_to_clipboard (void)
{
    char *liststr;
    int err = 0;

    liststr = main_window_selection_as_string();

    if (liststr != NULL) {
	delimiter_dialog();
	err = real_csv_to_clipboard(liststr);
	free(liststr);
    }

    return err;
}
