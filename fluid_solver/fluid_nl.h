#ifndef FLUID_NL_H_
#define FLUID_NL_H_

#define PI 3.14159265359

	int fluid_nl ( double *a, 
			 double *a_n, 
		   	 double *u, 
		   	 double *u_n, 
		   	 double *p, 
		   	 double *p_n, 
		   	 double *p_old, 
		   	 int	t, 
		   	 int	N, 
		   	 double kappa, 
		   	 double tau, 
		   	 double gamma
		      );

	int linsolve ( int n, 
			 double **A, 
			 double *b, 
			 double *x 
			);

#endif
