void reverse_gradient (double *g, int n);

void BFGS_get_user_values (double *b, int n, int *maxit,
			   double *reltol, gretlopt opt,
			   PRN *prn);

int BFGS_test (double *b, int n, int maxit, double reltol,
	       int *fncount, int *grcount, BFGS_CRIT_FUNC cfunc, 
	       int crittype, BFGS_GRAD_FUNC gradfunc, void *data, 
	       gretlopt opt, PRN *prn);

