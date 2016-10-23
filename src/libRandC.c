
#include "array.h"
#include "ClickClust.h"


void array1to2(int a, int b, double *y, double **x){

	int i, j, k;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){

			x[i][j] = y[k];
			k++;

		}
	}
	
	
}


void array1to3(int a, int b, int c, double *y, double ***x){

	int i, j, k, m;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){
			for (m=0; m<c; m++){

				x[i][j][m] = y[k];
				k++;
			
			}
		}
	}
	
	
}


void array2to1(int a, int b, double *y, double **x){

	int i, j, k;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){

			y[k] = x[i][j];
			k++;

		}
	}
	
	
}


void array3to1(int a, int b, int c, double *y, double ***x){

	int i, j, k, m;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){
			for (m=0; m<c; m++){

				y[k] = x[i][j][m];
				k++;
			
			}
		}
	}
	
	
}








void array1to2i(int a, int b, int *y, int **x){

	int i, j, k;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){

			x[i][j] = y[k];
			k++;

		}
	}
	
	
}


void array1to3i(int a, int b, int c, int *y, int ***x){

	int i, j, k, m;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){
			for (m=0; m<c; m++){

				x[i][j][m] = y[k];
				k++;
			
			}
		}
	}
	
	
}

void array1to3i_(int a, int b, int c, int *y, int ***x){

	int i, j, k, m;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){
			for (m=0; m<c; m++){

				x[i][m][j] = y[k];
				k++;
			
			}
		}
	}
	
	
}

void array2to1i(int a, int b, int *y, int **x){

	int i, j, k;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){

			y[k] = x[i][j];
			k++;

		}
	}
	
	
}


void array3to1i(int a, int b, int c, int *y, int ***x){

	int i, j, k, m;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){
			for (m=0; m<c; m++){

				y[k] = x[i][j][m];
				k++;
			
			}
		}
	}
	
	
}

void array3to1i_(int a, int b, int c, int *y, int ***x){

	int i, j, k, m;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){
			for (m=0; m<c; m++){

				y[k] = x[i][m][j];
				k++;
			
			}
		}
	}
	
	
}

