


#include<math.h>
#include "array.h"
#define inf 1e+40;
#include "ClickClust.h"




//	Analysis with prior probabilities estimated




double logL_kernel_(int p, int n, int K, int ***x, int *y, double *alpha, double **beta, double ***Pi, int scaleconst, int ntotal){

	int i, k;
	double sum1, sum2;

	sum1 = 0.0;
	for (i=0; i<n; i++){
		sum2 = 0.0;
		for (k=0; k<K; k++){
			sum2 = sum2 + f_kernel(p, k, i, x, Pi, scaleconst) * alpha[k] * beta[k][y[i]];
		}
		sum1 = sum1 + log(sum2);
	}
	
	return sum1 - ntotal * log(scaleconst);
	
}



void Estep_(int p, int n, int ***x, int *y, double *alpha, double **beta, double ***Pi, double **gamma, int K){
	
	int i, k, j, m, r;
	double prod;
								
	for (i=0; i<n; i++){
		for (k=0; k<K; k++){
			gamma[i][k] = 1.0;
			for (m=0; m<K; m++){
				if (m != k){
					prod = log(alpha[m]) - log(alpha[k]) + log(beta[m][y[i]]) - log(beta[k][y[i]]);
					for (j=0; j<p; j++){
						for (r=0; r<p; r++){	
							prod = prod + x[i][j][r] * (log(Pi[j][r][m]) - log(Pi[j][r][k]));
						}
					}
					prod = exp(prod);
					gamma[i][k] = gamma[i][k] + prod;
				}
			}
			gamma[i][k] = 1.0 / gamma[i][k];
		}
	}
			
}




void Mstep_(int p, int n, int ***x, int *y, double *alpha, double **beta, double ***Pi, double **gamma, int K, double lowbeta, double lowPi, int **nj){
	
	int i, j, r, k;
	
	double a, minPi, minbeta;
	double *gammasum, **denom;
	
	MAKE_VECTOR(gammasum, K);
	MAKE_MATRIX(denom, p, K);
	
	/* compute alpha's */
	
	for (k=0; k<K; k++){
		
		gammasum[k] = 0.0;
		
		for (i=0; i<n; i++){
			gammasum[k] = gammasum[k] + gamma[i][k];
		}
		
		alpha[k] = gammasum[k] / n;
	
	}


	/* compute beta's */

	
	for (k=0; k<K; k++){
		
		for (j=0; j<p; j++){
				beta[k][j] = 0.0;
		}
		
		for (i=0; i<n; i++){
			beta[k][y[i]] = beta[k][y[i]] + gamma[i][k];
		}
		
		for (j=0; j<p; j++){
			beta[k][j] = beta[k][j] / gammasum[k];
		}
	
	}
	

	if (lowbeta != 0.0){

	/* changing small beta's to "lowbeta's" */
	
		for (k=0; k<K; k++){
			minbeta = 1.0;
			for (j=0; j<p; j++){
				if (beta[k][j] < minbeta) minbeta = beta[k][j];
			}
				
			if (minbeta < lowbeta){
				a = (lowbeta - minbeta) / (1 - lowbeta * p);
				for (j=0; j<p; j++){
					beta[k][j] = (beta[k][j] + a) / (1 + a * p);
				}
			}
		}
	}	
	
	
	
	/* compute Pi's */
	
	for (j=0; j<p; j++){
		for (k=0; k<K; k++){
			denom[j][k] = 0.0;
			for (i=0; i<n; i++){	
				denom[j][k] = denom[j][k] + gamma[i][k] * nj[i][j];
			}
		}
	}
	
	for (j=0; j<p; j++){
		for (r=0; r<p; r++){
			for (k=0; k<K; k++){
				Pi[j][r][k] = 0.0;
				for (i=0; i<n; i++){	
					Pi[j][r][k] = Pi[j][r][k] + gamma[i][k] * x[i][j][r];
				}
				Pi[j][r][k] = Pi[j][r][k] / denom[j][k];
			}
		}
	}
				
	if (lowPi != 0.0){
	
	/* changing small Pi's to "lowPi" */

		for (k=0; k<K; k++){
			for (j=0; j<p; j++){
				minPi = 1.0;
				for (r=0; r<p; r++){
					if (Pi[j][r][k] < minPi) minPi = Pi[j][r][k];
				}
				
				if (minPi < lowPi){
					a = (lowPi - minPi) / (1 - lowPi * p);
					for (r=0; r<p; r++){
						Pi[j][r][k] = (Pi[j][r][k] + a) / (1 + a * p);
					}					
				}
			}
		}
	}
			
	FREE_VECTOR(gammasum);
	FREE_MATRIX(denom);
		
}




		
void init_(int p, int n, int K, int ***x, int *y, double *alpha, double **beta, double ***Pi, int h, double lowbeta, double lowPi, int **nj, int scaleconst, int ntotal, int shortem){
	
	int k, j, r, s, a, step;
	int *starts;
	
	double ll, bestll;
	double *bestalpha, **bestbeta, ***bestPi, **gamma1;
		
	MAKE_VECTOR(starts, K);
	
	MAKE_VECTOR(bestalpha, K);
	MAKE_MATRIX(bestbeta, K, p);
	MAKE_3ARRAY(bestPi, p, p, K);
	MAKE_MATRIX(gamma1, n, K);
	

	bestll = -inf;

	for (s=0; s<h; s++){
		
		srswor(K, n, starts);
		
		for (k=0; k<K; k++){
			
			alpha[k] = 1.0 / K;
			
			for (j=0; j<p; j++){
				
				beta[k][j] = 1.0 / p;
				
				a = 0;
				for (r=0; r<p; r++){
					if (nj[starts[k]][j] != 0)
						Pi[j][r][k] = 1.0 * x[starts[k]][j][r] / nj[starts[k]][j];
					else 
						Pi[j][r][k] = 1.0 / p;
					
					if (Pi[j][r][k] <= lowPi){
						Pi[j][r][k] = lowPi;
						a = a + 1;
					}
				}
						
				for (r=0; r<p; r++){
					if (Pi[j][r][k] != lowPi) Pi[j][r][k] = (Pi[j][r][k] -  a * lowPi) / (1 - a * lowPi);
				}					
					
			}
		}
		
	
		step = 0;
		while (step < shortem){

			Estep_(p, n, x, y, alpha, beta, Pi, gamma1, K);
			Mstep_(p, n, x, y, alpha, beta, Pi, gamma1, K, lowbeta, lowPi, nj);
			
			step++;
		}	

		ll = logL_kernel_(p, n, K, x, y, alpha, beta, Pi, scaleconst, ntotal);	
				
		if (ll > bestll){
			
			bestll = ll;
			cpy1(alpha, K, bestalpha);
			cpy2(beta, K, p, bestbeta);
			cpy3(Pi, p, p, K, bestPi);
			
		}
		
	}

	
	cpy1(bestalpha, K, alpha);
	cpy2(bestbeta, K, p, beta);
	cpy3(bestPi, p, p, K, Pi);

		
	FREE_VECTOR(starts);
	
	FREE_VECTOR(bestalpha);
	FREE_MATRIX(bestbeta);
	FREE_3ARRAY(bestPi);
	FREE_MATRIX(gamma1);
	
		
}




void EM_(int p, int n, int ***x, int *y, double *alpha, double **beta, double ***Pi, double **gamma, int *id, int K, int h, double tol, double *l, double lowbeta, double lowPi, int **nj, int scaleconst, int ntotal, int shortem){
	
	int step, i, k, M;
	double ll, llold, gmax;
	
	M = K - 1 + K * (p - 1) + K * p * (p - 1);
	
	step = 0;
		
	init_(p, n, K, x, y, alpha, beta, Pi, h, lowbeta, lowPi, nj, scaleconst, ntotal, shortem);
	
//	parprint(p, K, alpha, beta, Pi);
	
	llold = -inf;
	ll = logL_kernel_(p, n, K, x, y, alpha, beta, Pi, scaleconst, ntotal);
	
	while (fabs((ll - llold) / ll) > tol){
		
		step++;
		llold = ll;
		Estep_(p, n, x, y, alpha, beta, Pi, gamma, K);
		Mstep_(p, n, x, y, alpha, beta, Pi, gamma, K, lowbeta, lowPi, nj);
//		parprint(p, K, alpha, beta, Pi);		

		ll = logL_kernel_(p, n, K, x, y, alpha, beta, Pi, scaleconst, ntotal);
							
	}
				
	l[0] = ll;
	l[1] = BIC(M, n, K, l[0]);
		
	for (i=0; i<n; i++){
		gmax = gamma[i][0];
		id[i] = 0;
		for (k=1; k<K; k++){
			if (gamma[i][k] > gmax){
				gmax = gamma[i][k];
				id[i] = k;
			}
		}
	}
			
}


/*
void parprint(int p, int K, double *alpha, double **beta, double ***Pi){
	
	int k, j, r;
		
	printf("Mixing proportions:\n");
	for (k=0; k<K; k++){
		printf("%lf ", alpha[k]);
	}
	printf("\n\n");

	printf("Prior probabilities:\n");
	for (k=0; k<K; k++){
		for (j=0; j<p; j++){
			printf("%lf ", beta[k][j]);
		}
		printf("\n");
	}
	printf("\n\n");
	
	printf("Transition probabilities:\n");
	for (k=0; k<K; k++){
		printf("Cluster %i\n", k);
		for (r=0; r<p; r++){
			for (j=0; j<p; j++){
				printf("%lf ", Pi[j][r][k]);
			}
			printf("\n");
		}
		printf("\n");
	}
}



void gammaprint(int n, int K, double **gamma){
	
	int i, k;
	
	printf("Posterior probabilities:\n");
	for (i=0; i<n; i++){
		for (k=0; k<K; k++){
			printf("%lf ", gamma[i][k]);
		}
		printf("\n");
	}

}
*/
