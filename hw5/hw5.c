#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include<stdio.h>
#define NRANSI
#include "nrutil.h"
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
		//if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix");
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

#define NRANSI
#define TINY 1.0e-20

void ludcmp(float** a, int n, int* indx, float* d)
{
	int i, imax, j, k;
	float big, dum, sum, temp;
	float* vv;

	vv = vector(1, n);
	*d = 1.0;
	for (i = 1; i <= n; i++) {
		big = 0.0;
		for (j = 1; j <= n; j++)
			if ((temp = fabs(a[i][j])) > big) big = temp;
		//if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i] = 1.0 / big;
	}
	for (j = 1; j <= n; j++) {
		for (i = 1; i < j; i++) {
			sum = a[i][j];
			for (k = 1; k < i; k++) sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i <= n; i++) {
			sum = a[i][j];
			for (k = 1; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ((dum = vv[i] * fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 1; k <= n; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0) a[j][j] = TINY;
		if (j != n) {
			dum = 1.0 / (a[j][j]);
			for (i = j + 1; i <= n; i++) a[i][j] *= dum;
		}
	}
	free_vector(vv, 1, n);
}
#undef TINY
#undef NRANSI

void lubksb(float** a, int n, int* indx, float b[])
{
	int i, ii = 0, ip, j;
	float sum;

	for (i = 1; i <= n; i++) {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii)
			for (j = ii; j <= i - 1; j++) sum -= a[i][j] * b[j];
		else if (sum) ii = i;
		b[i] = sum;
	}
	for (i = n; i >= 1; i--) {
		sum = b[i];
		for (j = i + 1; j <= n; j++) sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
}

float pythag(a, b)
float a, b;
{
	float absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) return absa * sqrt(1.0 + SQR(absb / absa));
	else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}

#define NRANSI
void svdcmp(float** a, int m, int n, float w[], float** v)
{
	float pythag(float a, float b);
	int flag, i, its, j, jj, k, l, nm;
	float anorm, c, f, g, h, s, scale, x, y, z, * rv1;

	rv1 = vector(1, n);
	g = scale = anorm = 0.0;
	for (i = 1; i <= n; i++) {
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m) {
			for (k = i; k <= m; k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k = i; k <= m; k++) {
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][i] = f - g;
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = i; k <= m; k++) s += a[k][i] * a[k][j];
					f = s / h;
					for (k = i; k <= m; k++) a[k][j] += f * a[k][i];
				}
				for (k = i; k <= m; k++) a[k][i] *= scale;
			}
		}
		w[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m && i != n) {
			for (k = l; k <= n; k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k = l; k <= n; k++) {
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][l] = f - g;
				for (k = l; k <= n; k++) rv1[k] = a[i][k] / h;
				for (j = l; j <= m; j++) {
					for (s = 0.0, k = l; k <= n; k++) s += a[j][k] * a[i][k];
					for (k = l; k <= n; k++) a[j][k] += s * rv1[k];
				}
				for (k = l; k <= n; k++) a[i][k] *= scale;
			}
		}
		anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}
	for (i = n; i >= 1; i--) {
		if (i < n) {
			if (g) {
				for (j = l; j <= n; j++)
					v[j][i] = (a[i][j] / a[i][l]) / g;
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = l; k <= n; k++) s += a[i][k] * v[k][j];
					for (k = l; k <= n; k++) v[k][j] += s * v[k][i];
				}
			}
			for (j = l; j <= n; j++) v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for (i = IMIN(m, n); i >= 1; i--) {
		l = i + 1;
		g = w[i];
		for (j = l; j <= n; j++) a[i][j] = 0.0;
		if (g) {
			g = 1.0 / g;
			for (j = l; j <= n; j++) {
				for (s = 0.0, k = l; k <= m; k++) s += a[k][i] * a[k][j];
				f = (s / a[i][i]) * g;
				for (k = i; k <= m; k++) a[k][j] += f * a[k][i];
			}
			for (j = i; j <= m; j++) a[j][i] *= g;
		}
		else for (j = i; j <= m; j++) a[j][i] = 0.0;
		++a[i][i];
	}
	for (k = n; k >= 1; k--) {
		for (its = 1; its <= 30; its++) {
			flag = 1;
			for (l = k; l >= 1; l--) {
				nm = l - 1;
				if ((float)(fabs(rv1[l]) + anorm) == anorm) {
					flag = 0;
					break;
				}
				if ((float)(fabs(w[nm]) + anorm) == anorm) break;
			}
			if (flag) {
				c = 0.0;
				s = 1.0;
				for (i = l; i <= k; i++) {
					f = s * rv1[i];
					rv1[i] = c * rv1[i];
					if ((float)(fabs(f) + anorm) == anorm) break;
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for (j = 1; j <= m; j++) {
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y * c + z * s;
						a[j][i] = z * c - y * s;
					}
				}
			}
			z = w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j = 1; j <= n; j++) v[j][k] = -v[j][k];
				}
				break;
			}
			//if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
			x = w[l];
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = pythag(f, 1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0;
			for (j = l; j <= nm; j++) {
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				for (jj = 1; jj <= n; jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x * c + z * s;
					v[jj][i] = z * c - x * s;
				}
				z = pythag(f, h);
				w[j] = z;
				if (z) {
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for (jj = 1; jj <= m; jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y * c + z * s;
					a[jj][i] = z * c - y * s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	free_vector(rv1, 1, n);
}
#undef NRANSI

float lineq1[4][4] = {
	{4.0,2.0,3.0,-1.0},
	{-2.0,-1.0,-2.0,2.0},
	{5.0,3.0,4.0,-1.0},
	{11.0,4.0,6.0,1.0}
};
float b1[4] = { 4.0 ,-3.0,4.0,11.0 };

float lineq2[5][5] = {
	{2,-4,-5,5,0},
	{-1,1,2,0,4},
	{-1,6,0,3,2},
	{0,1,3,7,5},
	{5,0,8,7,-2}
};
float b2[5] = { -5,2,0,4,-1 };

float lineq3[6][6] = {
	{0.4, 8.2, 6.7, 1.9, 2.2, 5.3},
	{7.8, 8.3 ,7.7, 3.3, 1.9, 4.8},
	{5.5, 8.8 ,3.0, 1.0, 5.1, 6.4},
	{5.1 ,5.1 ,3.6, 5.8, 5.7, 4.9},
	{3.5, 2.7 ,5.7, 8.2, 9.6, 2.9},
	{3.0 ,5.3, 5.6, 3.5, 6.8, 5.7}
};
float b3[6] = { -2.9, -8.2, 7.7 ,-1.0,5.7, 3.0 };

#define NP 20
#define MP 20
#define MAXSTR 80
void by_gaussj() {
	int j, k, l, m, n;
	float** a, ** ai, ** u, ** b, ** x, ** t;
	char dummy[MAXSTR];
	FILE* fp;

	a = matrix(1, NP, 1, NP);
	ai = matrix(1, NP, 1, NP);
	u = matrix(1, NP, 1, NP);
	b = matrix(1, NP, 1, MP);
	x = matrix(1, NP, 1, MP);
	t = matrix(1, NP, 1, MP);
	if ((fp = fopen("input.dat", "r")) == NULL)
		nrerror("Data file matrx1.dat not found\n");
	int prob = 1;
	while (!feof(fp)) {
		fgets(dummy, MAXSTR, fp);
		fgets(dummy, MAXSTR, fp);
		fscanf(fp, "%d %d ", &n, &m);
		fgets(dummy, MAXSTR, fp);
		for (k = 1; k <= n; k++)
			for (l = 1; l <= n; l++) fscanf(fp, "%f ", &a[k][l]);
		fgets(dummy, MAXSTR, fp);
		for (l = 1; l <= m; l++)
			for (k = 1; k <= n; k++) fscanf(fp, "%f ", &b[k][l]);
		/* save matrices for later testing of results */
		for (l = 1; l <= n; l++) {
			for (k = 1; k <= n; k++) ai[k][l] = a[k][l];
			for (k = 1; k <= m; k++) x[l][k] = b[l][k];
		}
		/* invert matrix */
		gaussj(ai, n, x, m);
		printf("\n\t %d번째 A 행렬의 역행렬\n" , prob);
		for (int i = 1; i <= n; i++) {
			for (int j=1; j <= n; j++) {
				printf("\t %12.6f ", ai[i][j]);
			}
			printf("\n");
		}
		printf("\n\t %d번째의 solution vector\n" , prob++);
		for (int i = 1; i <= n; i++) {
			printf("\t %12.6f ", x[i][1]);
		}
	}
	printf("\n");
}

#define NRANSI
void mprove(float** a, float** alud, int n, int indx[], float b[], float x[])
{
	void lubksb(float** a, int n, int* indx, float b[]);
	int j, i;
	double sdp;
	float* r;

	r = vector(1, n);
	for (i = 1; i <= n; i++) {
		sdp = -b[i];
		for (j = 1; j <= n; j++) sdp += a[i][j] * x[j];
		r[i] = sdp;
	}
	lubksb(alud, n, indx, r);
	for (i = 1; i <= n; i++) x[i] -= r[i];
	free_vector(r, 1, n);
}
#undef NRANSI

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
float ran3(idum)
long* idum;
{
	static int inext, inextp;
	static long ma[56];
	static int iff = 0;
	long mj, mk;
	int i, ii, k;

	if (*idum < 0 || iff == 0) {
		iff = 1;
		mj = labs(MSEED - labs(*idum));
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for (i = 1; i <= 54; i++) {
			ii = (21 * i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if (mk < MZ) mk += MBIG;
			mj = ma[ii];
		}
		for (k = 1; k <= 4; k++)
			for (i = 1; i <= 55; i++) {
				ma[i] -= ma[1 + (i + 30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext = 0;
		inextp = 31;
		*idum = 1;
	}
	if (++inext == 56) inext = 1;
	if (++inextp == 56) inextp = 1;
	mj = ma[inext] - ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext] = mj;
	return mj * FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

void add_mprove() {
	int j, k, l, m, n, * indx;
	float d, * x, ** ainit, ** b, ** aa, ** a;
	char dummy[MAXSTR];
	FILE* fp;
	long idum = (-13);
	indx = ivector(1, NP);
	x = vector(1, NP);
	ainit = matrix(1, NP, 1, NP);
	b = matrix(1, NP, 1, NP);
	if ((fp = fopen("input.dat", "r")) == NULL)
		nrerror("Data file matrx1.dat not found\n");
	
	int prob = 1;
	
	while (!feof(fp)) {
		fgets(dummy, MAXSTR, fp);
		fgets(dummy, MAXSTR, fp);
		fscanf(fp, "%d %d ", &n, &m);
		fgets(dummy, MAXSTR, fp);
		for (k = 1; k <= n; k++)
			for (l = 1; l <= n; l++) fscanf(fp, "%f ", &ainit[k][l]);
		fgets(dummy, MAXSTR, fp);
		for (l = 1; l <= m; l++)
			for (k = 1; k <= n; k++) fscanf(fp, "%f ", &b[k][l]);
		//a = convert_matrix(&ainit[0][0], 1, n, 1, n);
		aa = matrix(1, n, 1, n);
		for (int i = 1; i <= n; i++) {
			x[i] = b[i][1];
			for (int j = 1; j <= n; j++) {
				aa[i][j] = ainit[i][j];
			}
		}
		ludcmp(aa, n, indx, &d);
		lubksb(aa, n, indx, x);

		printf("\n%d 번째 문제\n", prob++);
		printf("\nSolution vector for the equations:\n");
		for (int i = 1; i <= n; i++) printf("%12.6f", x[i]);
		printf("\n");
		/* now phoney up x and let mprove fix it */
		for (int i = 1; i <= n; i++) x[i] *= (1.0 + 0.2 * ran3(&idum));
		printf("\nSolution vector with noise added:\n");
		for (int i = 1; i <= n; i++) printf("%12.6f", x[i]);
		printf("\n");

		float * bb = vector(1, n);
		for (int i = 1; i <= n; i++) {
			bb[i] = b[i][1];
		}

		mprove(ainit, aa, n, indx, bb, x); //
		printf("\nSolution vector recovered by mprove:\n");
		for (int i = 1; i <= n; i++) printf("%12.6f", x[i]);
		printf("\n");
	}
	fclose(fp);
}

void by_ludcmp() {
	int j, k, l, m, n, * indx;
	float p, * x, ** a, ** b, ** c;
	char dummy[MAXSTR];
	FILE* fp;
	
	indx = ivector(1, NP);
	x = vector(1, NP);
	a = matrix(1, NP, 1, NP);
	b = matrix(1, NP, 1, NP);
	c = matrix(1, NP, 1, NP);
	if ((fp = fopen("input.dat", "r")) == NULL)
		nrerror("Data file matrx1.dat not found\n");
	int prob = 1;
	while (!feof(fp)) {
		fgets(dummy, MAXSTR, fp);
		fgets(dummy, MAXSTR, fp);
		fscanf(fp, "%d %d ", &n, &m);
		fgets(dummy, MAXSTR, fp);
		for (k = 1; k <= n; k++)
			for (l = 1; l <= n; l++) fscanf(fp, "%f ", &a[k][l]);
		fgets(dummy, MAXSTR, fp);
		for (l = 1; l <= m; l++)
			for (k = 1; k <= n; k++) fscanf(fp, "%f ", &b[k][l]);
		
		for (l = 1; l <= n; l++)
			for (k = 1; k <= n; k++) c[k][l] = a[k][l];
		
		ludcmp(c, n, indx, &p);
		
		for (k = 1; k <= m; k++) {
			for (l = 1; l <= n; l++) x[l] = b[l][k];
			lubksb(c, n, indx, x);

			printf("\n %d 번째의 solution vector\n", prob++);
			for (int i = 1; i <= n; i++) {
				printf("%12.6f ", x[i]);
			}
		}
	}
	fclose(fp);
}
void get_det() {
	int j, k, l, m, n, dum, * indx, * jndx;
	float d, ** a, ** xl, ** xu, ** x;
	char dummy[MAXSTR];
	FILE* fp;
	int prob = 1;

	indx = ivector(1, NP);
	jndx = ivector(1, NP);
	a = matrix(1, NP, 1, NP);
	xl = matrix(1, NP, 1, NP);
	xu = matrix(1, NP, 1, NP);
	x = matrix(1, NP, 1, NP);
	if ((fp = fopen("input.dat", "r")) == NULL)
		nrerror("Data file matrx1.dat not found\n");
	while (!feof(fp)) {
		fgets(dummy, MAXSTR, fp);
		fgets(dummy, MAXSTR, fp);
		fscanf(fp, "%d %d ", &n, &m);
		fgets(dummy, MAXSTR, fp);
		for (k = 1; k <= n; k++)
			for (l = 1; l <= n; l++) fscanf(fp, "%f ", &a[k][l]);
		fgets(dummy, MAXSTR, fp);
		for (l = 1; l <= m; l++)
			for (k = 1; k <= n; k++) fscanf(fp, "%f ", &x[k][l]);

		ludcmp(a, n, indx, &d);
		float det = 1;
		for (k = 1; k <= n; k++) {
			for (l = 1; l <= n; l++) {
				if (k == l)
					det *= a[k][l];
			}
		}
		printf("\n%d 번째 행렬의 determinant : %12.6f\n" , prob++ , det);
	}
	fclose(fp);
}
void main() {
	printf("\n------------------------------1) Gaussj 로 구한 solution 과 inverse------------------------------\n");
	by_gaussj();
	printf("\n------------------------------2) ludcmp 로 구한 solution-------------------------------\n");
	by_ludcmp();
	printf("\n------------------------------4) 각 행렬의 determinant------------------------------\n");
	get_det();

	printf("\n------------------------------5) mprove 적용------------------------------\n");
	add_mprove();
}
