#ifndef SYSML_H_
#define SYSML_H_

double *model_get_Xi (const MODEL *pmod, DATASET *dset, int i);

int fiml_driver (equation_system *sys, DATASET *dset, 
		 gretlopt opt, PRN *prn);

int liml_driver (equation_system *sys, DATASET *dset, 
		 PRN *prn);

#endif /* SYSML_H_ */
