int foo (void) {

                                /* 99%    97.5%  95%    90%    10%    5%    2.5%  1% */
    double t_crit_vals[5][8] = {{-3.75, -3.33, -3.00, -2.62, -0.37,  0.00, 0.34, 0.72}, /* T=25 */
				{-3.58, -3.22, -2.93, -2.60, -0.40, -0.03, 0.29, 0.66}, /* T=50 */
				{-3.51, -3.17, -2.89, -2.58, -0.42, -0.05, 0.26, 0.63}, /* T=100 */
				{-3.46, -3.14, -2.88, -2.57, -0.42, -0.06, 0.24, 0.62}, /* T=250 */
				{-3.44, -3.13, -2.87, -2.57, -0.43, -0.07, 0.24, 0.61}, /* T=500 */
				{-3.43, -3.12, -2.86, -2.57, -0.44, -0.07, 0.23, 0.60}} /* T>500 */



    row = (T > 500)? 5 : (T > 450)? 4 : (T > 240)? 3 : (T > 90)? 2 : 
	(T > 40)? 1 : (T > 24)? 0 : -1;
    if (row < 0) {
	sprintf(pval, "significance level unknown");
    } else {
	if (DFt < t_crit_vals[row][0] || DFt > t_crit_vals[row][7])
	    sprintf(pval, "significant at the 1 percent level");
	else if (DFt < t_crit_vals[row][1] || DFt > t_crit_vals[row][6])
	    sprintf(pval, "significant at the 2.5 percent level");
	else if (DFt < t_crit_vals[row][2] || DFt > t_crit_vals[row][5])
	    sprintf(pval, "significant at the 5 percent level");
	else if (DFt < t_crit_vals[row][3] || DFt > t_crit_vals[row][4])
	    sprintf(pval, "significant at the 10 percent level");
	else
	    sprintf(pval, "not significant at the 10 percent level");
    }
}
