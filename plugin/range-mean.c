#include "libgretl.h"

static void get_range_and_mean (int t1, int t2, double *x,
				double *range, double *mean)
{
    double me = 0.0, mi = x[t1], ma = x[t1];
    int t, n = 0;

    for (t=t1; t<=t2; t++) {
	if (na(x[t])) continue;
	ma = (x[t] > ma)? x[t] : ma;
	mi = (x[t] < mi)? x[t] : mi;
	me += x[t];
	n++;
    }

    if (n > 0) {
	*mean = me / n;
	*range = ma - mi;
    } else {
	*mean = NADBL;
	*range = NADBL;
    }
}


int range_mean_graph (int varnum, double **Z, DATAINFO *pdinfo, 
		      PRN *prn)
{
    double **rmZ;
    DATAINFO *rminfo;
    MODEL rmmod;
    int rmlist[4] = { 3, 1, 2, 0};
    int i, k, t, m, nsamp, err = 0;
    int start, end;
    double mean, range;

    nsamp = pdinfo->t2 - pdinfo->t1 + 1;

    if (nsamp < 16) {
	pprintf(prn, _("Sample is too small for range-mean graph\n"));
	errmsg(err, prn);
	return 1;
    } 	

    if (pdinfo->pd > 1 && nsamp >= 3 * pdinfo->pd) {
	k = pdinfo->pd;
    } else if (nsamp >= 30) {
	k = 10;
    } else {
	k = 5;
    }
	
    m = (nsamp / k) + ((nsamp % k)? 1 : 0);

    fprintf(stderr, "nsamp = %d, k = %d, m = %d\n",
	    nsamp, k, m);
    
    rminfo = create_new_dataset(&rmZ, 3, m, 0);
    if (rminfo == NULL) return E_ALLOC;
    rminfo->extra = 1;

    /* find group means and ranges */
    for (t=0; t<m; t++) {
	start = pdinfo->t1 + t * k;
	end = (start + k > pdinfo->t2)? pdinfo->t2 : start + k;
	get_range_and_mean(start, end, Z[varnum], &range, &mean);
	rmZ[1][t] = range;
	rmZ[2][t] = mean;
    }

    strcpy(rminfo->varname[1], "range");
    strcpy(rminfo->varname[2], "mean");

    rmmod = lsq(rmlist, &rmZ, rminfo, OLS, 0, 0.0);
    if ((err = rmmod.errcode)) {
	pprintf(prn, _("Error estimating range-mean model\n"));
	errmsg(err, prn);
    } else { /* just testing */
	printmodel(&rmmod, rminfo, prn);
    }

    clear_model(&rmmod, NULL);
    free_Z(rmZ, rminfo);
    clear_datainfo(rminfo, CLEAR_FULL);
    free(rminfo);

    return err;
}
