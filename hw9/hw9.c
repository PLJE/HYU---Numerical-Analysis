#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include "nrutil.h"
#include<time.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>  

#define NRANSI
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(float** a, int n, float** b, int m)
{
    int* indxc, * indxr, * ipiv;
    int i, icol, irow, j, k, l, ll;
    float big, dum, pivinv, temp;

    indxc = ivector(1, n);
    indxr = ivector(1, n);
    ipiv = ivector(1, n);
    for (j = 1; j <= n; j++) ipiv[j] = 0;
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++)
            if (ipiv[j] != 1)
                for (k = 1; k <= n; k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j][k]) >= big) {
                            big = fabs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
        ++(ipiv[icol]);
        if (irow != icol) {
            for (l = 1; l <= n; l++) SWAP(a[irow][l], a[icol][l])
                for (l = 1; l <= m; l++) SWAP(b[irow][l], b[icol][l])
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix");
        pivinv = 1.0 / a[icol][icol];
        a[icol][icol] = 1.0;
        for (l = 1; l <= n; l++) a[icol][l] *= pivinv;
        for (l = 1; l <= m; l++) b[icol][l] *= pivinv;
        for (ll = 1; ll <= n; ll++)
            if (ll != icol) {
                dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for (l = 1; l <= n; l++) a[ll][l] -= a[icol][l] * dum;
                for (l = 1; l <= m; l++) b[ll][l] -= b[icol][l] * dum;
            }
    }
    for (l = n; l >= 1; l--) {
        if (indxr[l] != indxc[l])
            for (k = 1; k <= n; k++)
                SWAP(a[k][indxr[l]], a[k][indxc[l]]);
    }
    free_ivector(ipiv, 1, n);
    free_ivector(indxr, 1, n);
    free_ivector(indxc, 1, n);
}
#undef SWAP
#undef NRANSI

void execute(const char* fname) {
    int N = 77;

    float** F = matrix(1, N, 1, 3);
    float** Y = matrix(1, N, 1, 2);

    FILE* fp = fopen(fname, "r");
    if (fp == NULL) {
        printf("failed to dopen file\n");
    }
    printf("\n%s parameters\n", fname);
    int i = 1;
    float x, y, xi, yi;
    while (fscanf(fp, "%f %f %f %f", &x, &y, &xi, &yi) == 4) {
        F[i][1] = x;
        F[i][2] = y;
        F[i][3] = 1;
        Y[i][1] = xi;
        Y[i][2] = yi;
        i++;
    }

    float** A = matrix(1, 3, 1, 3);
    float** b = matrix(1, 3, 1, 2);

    for (int i = 1; i <= 3; i++) { //F F^t a = F^t y
        for (int j = 1; j <= 3; j++) {
            A[i][j] = 0; //F F^t
            for (int k = 1; k <= N; k++) {
                A[i][j] += F[k][i] * F[k][j];
            }
            //printf("%f ",A[i][j]);
        }
        //printf("\n");
    }
    //printf("--------------------");
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 2; j++) {
            b[i][j] = 0; // F^t y 
            for (int k = 1; k <= N; k++) {
                b[i][j] += F[k][i] * Y[k][j];
            }
            //printf("%f ", b[i][j]);
        }
        //printf("\n");
    }
    gaussj(A, 3, b, 2);

    int cnt = 1;
    for (int i = 1; i <= 2; i++) {
        for (int j = 1; j <= 3; j++) { // a
            printf("a%d : %f\n", cnt, b[j][i]);
            cnt++;
        }
    }
    /*for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= 2; j++) {
            printf("%f ", y[i][j]);
        }
        printf("\n");
    }*/
    fclose(fp);
}

int main(void) {

    execute("fitdata1.dat");
    execute("fitdata2.dat");
    execute("fitdata3.dat");

    return 0;
}
