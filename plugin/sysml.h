#ifndef SYSML_H_
#define SYSML_H_

double *model_get_Xi (const MODEL *pmod, double **Z, int i);

int fiml_driver (equation_system *sys, double **Z, 
		 DATAINFO *pdinfo, gretlopt opt, PRN *prn);

int liml_driver (equation_system *sys, double **Z, 
		 DATAINFO *pdinfo, PRN *prn);

#endif /* SYSML_H_ */
