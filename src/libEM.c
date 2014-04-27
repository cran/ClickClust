

#include<math.h>
#include "array.h"
#define inf 1e+40;
#include "ClickClust.h"



//	Analysis with no prior probabilities


double logL_kernel(int p, int n, int K, int ***x, double *alpha, double ***Pi, int scaleconst, int ntotal){

	int i, k;
	double sum1, sum2;

	sum1 = 0.0;
	for (i=0; i<n; i++){
		sum2 = 0.0;
		for (k=0; k<K; k++){
			sum2 = sum2 + f_kernel(p, k, i, x, Pi, scaleconst) * alpha[k];
		}
		sum1 = sum1 + log(sum2);
	}
	
	sum1 = sum1 - n * log(p);
	
	return sum1 - ntotal * log(scaleconst);
	
}



void Estep(int p, int n, int ***x, double *alpha, double ***Pi, double **gamma, int K){
	
	int i, k, j, m, r;
	double prod;
								
	for (i=0; i<n; i++){
		for (k=0; k<K; k++){
			gamma[i][k] = 1.0;
			for (m=0; m<K; m++){
				if (m != k){
					prod = log(alpha[m]) - log(alpha[k]);
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




void Mstep(int p, int n, int ***x, double *alpha, double ***Pi, double **gamma, int K, double lowPi, int **nj){
	
	int i, j, r, k;
	
	double a, minPi;
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




		
void init(int p, int n, int K, int ***x, double *alpha, double ***Pi, int h, double lowPi, int **nj, int scaleconst, int ntotal, int shortem){
	
	int k, j, r, s, a, step;
	int *starts;
	
	double ll, bestll;
	double *bestalpha, ***bestPi, **gamma1;
		
	MAKE_VECTOR(starts, K);
	
	MAKE_VECTOR(bestalpha, K);
	MAKE_3ARRAY(bestPi, p, p, K);
	MAKE_MATRIX(gamma1, n, K);
	

	bestll = -inf;

	for (s=0; s<h; s++){
		
		for (k=0; k<K; k++){
			alpha[k] = 1.0 / K;
		}
		
		srswor(K, n, starts);
		
		for (k=0; k<K; k++){
			for (j=0; j<p; j++){
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

			Estep(p, n, x, alpha, Pi, gamma1, K);
			Mstep(p, n, x, alpha, Pi, gamma1, K, lowPi, nj);
			
			step++;
		}	

		ll = logL_kernel(p, n, K, x, alpha, Pi, scaleconst, ntotal);	
				
		if (ll > bestll){
			
			bestll = ll;
			cpy1(alpha, K, bestalpha);
			cpy3(Pi, p, p, K, bestPi);
			
		}
		
	}

	
	cpy1(bestalpha, K, alpha);
	cpy3(bestPi, p, p, K, Pi);

		
	FREE_VECTOR(starts);
	
	FREE_VECTOR(bestalpha);
	FREE_3ARRAY(bestPi);
	FREE_MATRIX(gamma1);
	
		
}




void EM(int p, int n, int ***x, double *alpha, double ***Pi, double **gamma, int *id, int K, int h, double tol, double *l, double lowPi, int **nj, int scaleconst, int ntotal, int shortem){
	
	int step, i, k, M;
	double ll, llold, gmax;
	
	M = K - 1 + K * p * (p - 1);
	
	step = 0;
		
	init(p, n, K, x, alpha, Pi, h, lowPi, nj, scaleconst, ntotal, shortem);
	
	llold = -inf;
	ll = logL_kernel(p, n, K, x, alpha, Pi, scaleconst, ntotal);
	

	while ((ll - llold) / fabs(ll) > tol){
		
		step++;
		llold = ll;
			
		Estep(p, n, x, alpha, Pi, gamma, K);
		
		Mstep(p, n, x, alpha, Pi, gamma, K, lowPi, nj);

		ll = logL_kernel(p, n, K, x, alpha, Pi, scaleconst, ntotal);
			
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





