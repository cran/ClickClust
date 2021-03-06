#ifndef CLICKCLUST_H
#define CLICKCLUST_H

void cpy1(double *a, int n, double *b);
void cpy2(double **a, int nrows, int ncols, double **b);
void cpy3(double ***a, int nrows, int ncols, int nslices, double ***b);
void srswor(int M, int n, int *y);
double BIC(int M, int n, int K, double ll);
void Estep(int p, int n, int ***x, double *alpha, double ***Pi, double **gamma, int K);
void Estep_(int p, int n, int ***x, int *y, double *alpha, double **beta, double ***Pi, double **gamma, int K);
void Mstep(int p, int n, int ***x, double *alpha, double ***Pi, double **gamma, int K, double lowPi, int **nj);
void Mstep_(int p, int n, int ***x, int *y, double *alpha, double **beta, double ***Pi, double **gamma, int K, double lowbeta, double lowPi, int **nj);
void init(int p, int n, int K, int ***x, double *alpha, double ***Pi, int h, double lowPi, int **nj, int scaleconst, int ntotal, int shortem);
void init_(int p, int n, int K, int ***x, int *y, double *alpha, double **beta, double ***Pi, int h, double lowbeta, double lowPi, int **nj, int scaleconst, int ntotal, int shortem);
void EM(int p, int n, int ***x, double *alpha, double ***Pi, double **gamma, int *id, int Kmax, int h, double tol, double *l, double lowPi, int **nj, int scaleconst, int ntotal, int shortem);
void EM_(int p, int n, int ***x, int *y, double *alpha, double **beta, double ***Pi, double **gamma, int *id, int K, int h, double tol, double *l, double lowbeta, double lowPi, int **nj, int scaleconst, int ntotal, int shortem);
double logL_kernel(int p, int n, int K, int ***x, double *alpha, double ***Pi, int scaleconst, int ntotal);
double logL_kernel_(int p, int n, int K, int ***x, int *y, double *alpha, double **beta, double ***Pi, int scaleconst, int ntotal);
double f_kernel(int p, int k, int i, int ***x, double ***Pi, int scaleconst);
void array1to2(int a, int b, double *y, double **x);
void array1to3(int a, int b, int c, double *y, double ***x);
void array2to1(int a, int b, double *y, double **x);
void array3to1(int a, int b, int c, double *y, double ***x);
void array1to2i(int a, int b, int *y, int **x);
void array1to3i_(int a, int b, int c, int *y, int ***x);
void array2to1i(int a, int b, int *y, int **x);
void array3to1i_(int a, int b, int c, int *y, int ***x);

#endif /* CLICKCLUST_H */
