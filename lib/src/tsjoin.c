/* "join" when both datasets are recognized as time series;
   we're looking at annual, quarterly or monthly data; and
   the LHS frequency is greater than the RHS. The cases
   handled are:

   monthly   <- quarterly
   monthly   <- annual
   quarterly <- annual

   Low frequency values are repeated.
*/

typedef struct ts_joiner_ ts_joiner;

struct ts_joiner_ {
    int t1;
    int t2;
    int pdr;
    int rt1;
    int rminor;
};

static int q2m (int q)
{
    return 3 * (q-1) + 1;
}

static int m2q (int m)
{
    return (m + 2) / 3;
}

static int mofq (int m)
{
    m = m % 3;
    return m == 0 ? 3 : m;
}

static guint32 eday (char *s, int t, DATASET *dset)
{
    int y, m = 1;

    ntolabel(s, t, dset);
    y = atoi(s);

    if (dset->pd != 1) {
	s = strchr(s, ':') + 1;
	if (dset->pd == 12) {
	    m = atoi(s);
	} else {
	    m = q2m(atoi(s));
	}
    }

    return epoch_day_from_ymd(y, m, 1);
}

static int ym2obs (char *obs, int y, int m, DATASET *dset)
{
    if (dset->pd == 1) {
	sprintf(obs, "%d", y);
    } else if (dset->pd == 12) {
	sprintf(obs, "%d:%02d", y, m);
    } else {
	sprintf(obs, "%d:%d", y, m2q(m));
    }
    return dateton(obs, dset);
}

static int get_rminor (DATASET *ldset, DATASET *rdset, int m)
{
    if (ldset->pd == 12) {
	return (rdset->pd == 1)? m : mofq(m);
    } else {
	/* quarterly on left */
	return m2q(m);
    }
}

static int fill_ts_joiner (DATASET *ldset, DATASET *rdset,
			   ts_joiner *tjr)
{
    guint32 edl, edr;
    int y, m, d;
    char obsl[OBSLEN];
    char obsr[OBSLEN];

    /* frequency ratio */
    tjr->pdr = ldset->pd / rdset->pd;

    /* start at the later of ldset->t1 and the first
       observation in rdset
    */
    edl = eday(obsl, ldset->t1, ldset);
    edr = eday(obsr, rdset->t1, rdset);
    ymd_bits_from_epoch_day(MAX(edl, edr), &y, &m, &d);
    tjr->t1  = ym2obs(obsl, y, m, ldset);
    tjr->rt1 = ym2obs(obsr, y, m, rdset);
    tjr->rminor = get_rminor(ldset, rdset, m);

    /* stop at the earlier of ldset->t2 and the last
       observation in rdset
    */
    edl = eday(obsl, ldset->t2, ldset);
    edr = eday(obsr, rdset->t2, rdset);
    ymd_bits_from_epoch_day(MIN(edl, edr), &y, &m, &d);
    tjr->t2 = ym2obs(obsl, y, m, ldset);

    /* possible adjustment of the endpoint in case the outer
       dataset is of lower frequency (2021-12-23)
    */
    if (tjr->pdr > 1 && m < 12) {
	int alt_t2 = tjr->t2 + tjr->pdr - 1;

	tjr->t2 = MIN(ldset->t2, alt_t2);
    }

    return 0;
}

static int use_tsjoin (const DATASET *l_dset,
		       const DATASET *r_dset)
{
    if (l_dset->pd == 12) {
	return r_dset->pd == 1 || r_dset->pd == 4;
    } else if (l_dset->pd == 4) {
	return r_dset->pd == 1;
    } else {
	return 0;
    }
}

static int tj_continue (ts_joiner *tjr, int m, int *ps)
{
    if (m == tjr->pdr) {
	m = 1;
	*ps += 1;
    } else {
	m++;
    }

    return m;
}
