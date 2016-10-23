
#include "array.h"
#include "ClickClust.h"

void runClickClust(int (*p1), int (*K1), int (*n1), int *x1, double *alpha, double *Pi1, double *gamma1, int *id, double (*e1), int (*r1), int (*shortem1), double *l, double (*lowPi1), double (*scaleconst1)){

	int i, j, r, shortem, p, K, n, h, ntotal;
	int ***x, **nj;
	
	double tol, lowPi, scaleconst;
	double ***Pi, **gamma;
	
	p = (*p1);
	K = (*K1);
	n = (*n1);
	h = (*r1);
	shortem = (*shortem1);
	tol = (*e1);
	lowPi = (*lowPi1);
	scaleconst = (*scaleconst1);
	
		
	MAKE_3ARRAY(x, n, p, p);
	MAKE_3ARRAY(Pi, p, p, K);
	
	MAKE_MATRIX(gamma, n, K);
	MAKE_MATRIX(nj, n, p);

	array1to3i_(n, p, p, x1, x);
	array1to3(p, p, K, Pi1, Pi);
	array1to2(n, K, gamma1, gamma);	
	
	ntotal = 0.0;
	for (i=0; i<n; i++){
		for (j=0; j<p; j++){
			nj[i][j] = 0;
			for (r=0; r<p; r++){
				ntotal = ntotal + x[i][j][r];
				nj[i][j] = nj[i][j] + x[i][j][r];
			}
		}
	}
	
	EM(p, n, x, alpha, Pi, gamma, id, K, h, tol, l, lowPi, nj, scaleconst, ntotal, shortem);

	array3to1i_(n, p, p, x1, x);
	array3to1(p, p, K, Pi1, Pi);
	array2to1(n, K, gamma1, gamma);
	
	FREE_MATRIX(nj);
	FREE_3ARRAY(Pi);
	FREE_MATRIX(gamma);
	FREE_3ARRAY(x);
		
}


void runClickClust_(int (*p1), int (*K1), int (*n1), int *x1, int *y, double *alpha, double *beta1, double *Pi1, double *gamma1, int *id, double (*e1), int (*r1), int (*shortem1), double *l, double (*lowbeta1), double (*lowPi1), double (*scaleconst1)){

	int i, j, r, shortem, p, K, n, h, ntotal;
	int ***x, **nj;
	
	double tol, lowbeta, lowPi, scaleconst;
	double **beta, ***Pi, **gamma;
	
	p = (*p1);
	K = (*K1);
	n = (*n1);
	h = (*r1);
	shortem = (*shortem1);
	tol = (*e1);
	lowbeta = (*lowbeta1);
	lowPi = (*lowPi1);
	scaleconst = (*scaleconst1);
	
		
	MAKE_3ARRAY(x, n, p, p);
	MAKE_MATRIX(beta, K, p);
	MAKE_3ARRAY(Pi, p, p, K);
	
	MAKE_MATRIX(gamma, n, K);
	MAKE_MATRIX(nj, n, p);

	array1to3i_(n, p, p, x1, x);
	array1to2(K, p, beta1, beta);
	array1to3(p, p, K, Pi1, Pi);
	array1to2(n, K, gamma1, gamma);
	
	ntotal = 0.0;
	for (i=0; i<n; i++){
		for (j=0; j<p; j++){
			nj[i][j] = 0;
			for (r=0; r<p; r++){
				ntotal = ntotal + x[i][j][r];
				nj[i][j] = nj[i][j] + x[i][j][r];
			}
		}
	}
	
	EM_(p, n, x, y, alpha, beta, Pi, gamma, id, K, h, tol, l, lowbeta, lowPi, nj, scaleconst, ntotal, shortem);

	array3to1i_(n, p, p, x1, x);
	array2to1(K, p, beta1, beta);
	array3to1(p, p, K, Pi1, Pi);
	array2to1(n, K, gamma1, gamma);
	
	FREE_MATRIX(nj);
	FREE_MATRIX(beta);
	FREE_3ARRAY(Pi);
	FREE_MATRIX(gamma);
	FREE_3ARRAY(x);
		
}
