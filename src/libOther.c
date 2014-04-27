

#include<math.h>
#include "array.h"
#define inf 1e+40;
#include "ClickClust.h"

#ifdef __HAVE_R_
	#include <R.h>
	#include <Rmath.h>
#endif


/* multiplies matrix by itself */

void MatrixProd(double **OO, int p, int m, double **Res){
     
     int i, j, k;

     for (i=0; i<p; i++){
         for (j=0; j<p; j++){
             Res[i][j] = 0.0;
             for (k=0; k<m; k++){
                 Res[i][j] = Res[i][j] + OO[i][k] * OO[j][k];
             }
         }     
     }
     
}


/*copies vector A to vector B */

void cpy1(double *a, int n, double *b){
	
	int i;

	for(i=0; i<n; i++) {
		b[i] = a[i];
	}

}


/*copies matrix A to matrix B*/

void cpy2(double **a, int nrows, int ncols, double **b){
	
	int i, j;

	for(i=0; i<nrows; i++) {
		for (j=0; j<ncols; j++) {
			b[i][j] = a[i][j];
		}
	}

}


/*copies cube A to cube B*/

void cpy3(double ***a, int nrows, int ncols, int nslices, double ***b){

	int i, j, k;

	for(i=0; i<nrows; i++) {
		for (j=0; j<ncols; j++) {
			for (k=0; k<nslices; k++) {
				b[i][j][k] = a[i][j][k];
			}
		}
	}

}


/* sampling without replacement */

void srswor(int M, int n, int *y){

	int flag;
	int k, v;
	int *indy;

	MAKE_VECTOR(indy, n);   
  
    for (k=0; k<n; k++) indy[k]=0;

	for (k=0; k<M; k++){
	
		flag=0;
	
		while (flag==0){
			
			#ifdef __HAVE_R_
				v = floor(runif(0.0, n));
			#else
				v = rand() % n;
			#endif			
			
			
			if (indy[v] == 0){
				y[k] = v;
				indy[v] = 1;
                flag = 1;
            }
        }

	}

	FREE_VECTOR(indy);
        
	return;

}




double BIC(int M, int n, int K, double ll){

	return -2 * ll + M * log(n);
		
}




double f_kernel(int p, int k, int i, int ***x, double ***Pi, int scaleconst){
	
	// scaleconst - Constant used for scaling to avoid numeric issues
	// if scaleconst = 1, no scaling is used
	
	int j, r;
	
	double sum;

	sum = 0.0;
	for (j=0; j<p; j++){
		for (r=0; r<p; r++){
			sum = sum + x[i][j][r] * log(Pi[j][r][k] * scaleconst);
		}
	}
	
//	printf("%lf ", exp(sum));
	
	return exp(sum);
	
}

