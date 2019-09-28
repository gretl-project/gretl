/*
   This is Mike Giles's Poisson inverse code, under GPL. See
   https://people.maths.ox.ac.uk/gilesm/codes/poissinv/

   The method is explicated in Giles's paper "Algorithm 955:
   Approximation of the Inverse Poisson Cumulative Distribution
   Function", ACM Transactions on Mathematical Software, Volume 42
   Issue 1, February 2016.
*/

/*
  This double precision function computes the inverse
  of the Poisson CDF

  u = CDF value in range (0,1)
  lam = Poisson rate

  For lam < 1e15,  max |error| no more than 1
  ave |error| < 1e-16*max(4,lam) for lam < 1e9
              < 1e-6             for lam < 1e15

  For lam > 1e15, the errors will be about 1 ulp.
*/

static double poissinv_core (double U, double V, double Lam)
{
    double Xi, W, T, Del, R, R2, S, S2, Eta, B0, B1;
    double Lami = 1/Lam;
    double X = 0.0;
    int i;

    if (U == 0.0) return 0.0;
    if (V == 0.0) return FP_INFINITE;
    if (!(U > 0.0 && V > 0.0)) return FP_NAN;

    if (Lam > 4.0) {
	W = normal_cdf_inverse(fmin(U, V));
	if (U > V) W = -W;

	/* use polynomial approximations in central region */
	if (fabs(W) < 3.0) {
	    double Lam_root = sqrt(Lam);

	    S = Lam_root*W + (1.0/3.0 + (1.0/6.0)*W*W)*(1.0 - W/(12.0*Lam_root));

	    Del = (1.0 /160.0);
	    Del = (1.0 / 80.0) + Del*(W*W);
	    Del = (1.0 / 40.0) + Del*(W*W);
	    Del = Del * Lami;

	    S = Lam + (S + Del);
	} else {
	    /* otherwise use Newton iteration */
	    S = W / sqrt(Lam);
	    R = 1.0 + S;
	    if (R < 0.1) R = 0.1;

	    do {
		T = log(R);
		R2 = R;
		S2 = sqrt(2.0*((1.0-R) + R*T));
		if (R < 1.0) S2 = -S2;
		R = R2 - (S2-S)*S2/T;
		if (R < 0.1*R2) R = 0.1*R2;
	    } while (fabs(R-R2) > 1e-8);

	    T = log(R);
	    S = Lam*R + log(sqrt(2.0*R*((1.0-R)+R*T))/fabs(R-1.0)) / T;
	    S = S - 0.0218/(S+0.065*Lam);
	    Del = 0.01/S;
	    S = S + Del;
	}

	/* if x > 10, round down to nearest integer, and check accuracy */
	X = floor(S);

	if (S > 10.0 && S < X + 2.0*Del) {
	    /* correction procedure based on Temme approximation */
	    if (X > 0.5*Lam && X < 2.0*Lam) {
		Xi = 1.0 / X;
		Eta = X * Lami;
		Eta = sqrt(2.0*(1.0-Eta+Eta*log(Eta))/Eta);
		if (X > Lam) Eta = -Eta;

		B1 =  8.0995211567045583e-16;              S = B1;
		B0 = -1.9752288294349411e-15;              S = B0 + S*Eta;
		B1 = -5.1391118342426808e-16 + 25.0*B1*Xi; S = B1 + S*Eta;
		B0 =  2.8534893807047458e-14 + 24.0*B0*Xi; S = B0 + S*Eta;
		B1 = -1.3923887224181616e-13 + 23.0*B1*Xi; S = B1 + S*Eta;
		B0 =  3.3717632624009806e-13 + 22.0*B0*Xi; S = B0 + S*Eta;
		B1 =  1.1004392031956284e-13 + 21.0*B1*Xi; S = B1 + S*Eta;
		B0 = -5.0276692801141763e-12 + 20.0*B0*Xi; S = B0 + S*Eta;
		B1 =  2.4361948020667402e-11 + 19.0*B1*Xi; S = B1 + S*Eta;
		B0 = -5.8307721325504166e-11 + 18.0*B0*Xi; S = B0 + S*Eta;
		B1 = -2.5514193994946487e-11 + 17.0*B1*Xi; S = B1 + S*Eta;
		B0 =  9.1476995822367933e-10 + 16.0*B0*Xi; S = B0 + S*Eta;
		B1 = -4.3820360184533521e-09 + 15.0*B1*Xi; S = B1 + S*Eta;
		B0 =  1.0261809784240299e-08 + 14.0*B0*Xi; S = B0 + S*Eta;
		B1 =  6.7078535434015332e-09 + 13.0*B1*Xi; S = B1 + S*Eta;
		B0 = -1.7665952736826086e-07 + 12.0*B0*Xi; S = B0 + S*Eta;
		B1 =  8.2967113409530833e-07 + 11.0*B1*Xi; S = B1 + S*Eta;
		B0 = -1.8540622107151585e-06 + 10.0*B0*Xi; S = B0 + S*Eta;
		B1 = -2.1854485106799979e-06 +  9.0*B1*Xi; S = B1 + S*Eta;
		B0 =  3.9192631785224383e-05 +  8.0*B0*Xi; S = B0 + S*Eta;
		B1 = -0.00017875514403292177 +  7.0*B1*Xi; S = B1 + S*Eta;
		B0 =  0.00035273368606701921 +  6.0*B0*Xi; S = B0 + S*Eta;
		B1 =   0.0011574074074074078 +  5.0*B1*Xi; S = B1 + S*Eta;
		B0 =   -0.014814814814814815 +  4.0*B0*Xi; S = B0 + S*Eta;
		B1 =    0.083333333333333329 +  3.0*B1*Xi; S = B1 + S*Eta;
		B0 =    -0.33333333333333331 +  2.0*B0*Xi; S = B0 + S*Eta;
		S  = S / (1.0 + B1*Xi);

		S = S*exp(-0.5*X*Eta*Eta)/sqrt(2.0*3.141592653589793*X);
		if (X < Lam) {
		    S += 0.5 * erfc(Eta * sqrt(0.5*X));
		    if (S > U) X -= 1.0;
		} else {
		    S -= 0.5 * erfc(-Eta * sqrt(0.5*X));
		    if (S > -V) X -= 1.0;
		}
	    } else {
		/* sum downwards or upwards */
		Xi = 1.0 / X;
		S = -(691.0/360360.0);
		S =  (1.0/1188.0) + S*Xi*Xi;
		S = -(1.0/1680.0) + S*Xi*Xi;
		S =  (1.0/1260.0) + S*Xi*Xi;
		S = -(1.0/360.0)  + S*Xi*Xi;
		S =  (1.0/12.0)   + S*Xi*Xi;
		S = S*Xi;
		S = (X - Lam) - X * log(X*Lami) - S;

		if (X < Lam) {
		    T = exp(-0.5*S);
		    S = 1.0 - T*(U*T) * sqrt(2.0*3.141592653589793*Xi) * Lam;
		    T = 1.0;
		    Xi = X;
		    for (i=1; i<50; i++) {
			Xi -= 1.0;
			T *= Xi*Lami;
			S += T;
		    }
		    if (S > 0.0) X -= 1.0;
		} else {
		    T = exp(-0.5*S);
		    S = 1.0 - T*(V*T) * sqrt(2.0*3.141592653589793*X);
		    Xi = X;
		    for (i=0; i<50; i++) {
			Xi += 1.0;
			S = S*Xi*Lami + 1.0;
		    }
		    if (S < 0.0) {
			X -= 1.0;
		    }
		}
	    }
	}
    }

    /* bottom-up summation */
    if (X < 10.0) {
	X = 0.0;
	T = exp(0.5*Lam);
	Del = 0.0;
	if (U > 0.5) {
	    Del = T*(1e-13*T);
	}
	S = 1.0 - T*(U*T) + Del;

	while (S < 0.0) {
	    X += 1.0;
	    T = X*Lami;
	    Del = T*Del;
	    S = T*S + 1.0;
	}

	/* top-down summation if needed */
	if (S < 2.0*Del) {
	    Del = 1e13*Del;
	    T = 1e17*Del;
	    Del = V*Del;

	    while (Del < T) {
		X += 1.0;
		Del *= X*Lami;
	    }

	    S = Del;
	    T = 1.0;
	    while (S > 0.0) {
		T *= X*Lami;
		S -= T;
		X -= 1.0;
	    }
	}
    }

    return X;
}

static double poissinv (double U, double Lam)
{
    return poissinv_core(U, 1.0 - U, Lam);
}
