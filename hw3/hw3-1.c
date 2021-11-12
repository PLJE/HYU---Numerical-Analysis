#include<stdio.h>
#include <math.h>

float bessj0(float x)
{
	float ax, z;
	double xx, y, ans, ans1, ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x * x;
		ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
			+ y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
		ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
			+ y * (59272.64853 + y * (267.8532712 + y * 1.0))));
		ans = ans1 / ans2;
	}
	else {
		z = 8.0 / ax;
		y = z * z;
		xx = ax - 0.785398164;
		ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
			+ y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
		ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
			+ y * (-0.6911147651e-5 + y * (0.7621095161e-6
				- y * 0.934945152e-7)));
		ans = sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
	}
	return ans;
}
float bessj1(float x)
{
	float ax, z;
	double xx, y, ans, ans1, ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x * x;
		ans1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
			+ y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
		ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
			+ y * (99447.43394 + y * (376.9991397 + y * 1.0))));
		ans = ans1 / ans2;
	}
	else {
		z = 8.0 / ax;
		y = z * z;
		xx = ax - 2.356194491;
		ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
			+ y * (0.2457520174e-5 + y * (-0.240337019e-6))));
		ans2 = 0.04687499995 + y * (-0.2002690873e-3
			+ y * (0.8449199096e-5 + y * (-0.88228987e-6
				+ y * 0.105787412e-6)));
		ans = sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}
void derivative(float x, float* fn, float* df) {
	//bessel functions can be found in ch.6 of NR in C (bessj0.c , bessj1.c)
	*fn = bessj0(x);
	*df = -bessj1(x);
}

void zbrak(float (*fx)(float), float x1, float x2, int n, float xb1[],
	float xb2[], int* nb)
{
	int nbb, i;
	float x, fp, fc, dx;

	nbb = 0;
	dx = (x2 - x1) / n;
	fp = (*fx)(x = x1);
	for (i = 1; i <= n; i++) {
		fc = (*fx)(x += dx);
		if (fc * fp <= 0.0) {
			xb1[++nbb] = x - dx;
			xb2[nbb] = x;
			if (*nb == nbb) return;

		}
		fp = fc;
	}
	*nb = nbb;
} 

#define JMAX 40
float rtbis(float (*func)(float), float x1, float x2, float xacc)
{ //bisection rtbis
	int j;
	float dx, f, fmid, xmid, rtb;

	f = (*func)(x1);
	fmid = (*func)(x2);
	rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
	for (j = 1; j <= JMAX; j++) {
		fmid = (*func)(xmid = rtb + (dx *= 0.5));
		if (fmid <= 0.0) rtb = xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	return 0.0;
}
#undef JMAX //cal1

#define MAXIT 30
float rtflsp(float (*func)(float), float x1, float x2, float xacc)
{
	int j;
	float fl, fh, xl, xh, swap, dx, del, f, rtf;

	fl = (*func)(x1);
	fh = (*func)(x2);
	if (fl < 0.0) {
		xl = x1;
		xh = x2;
	}
	else {
		xl = x2;
		xh = x1;
		swap = fl;
		fl = fh;
		fh = swap;
	}
	dx = xh - xl;
	for (j = 1; j <= MAXIT; j++) {
		rtf = xl + dx * fl / (fl - fh);
		f = (*func)(rtf);
		if (f < 0.0) {
			del = xl - rtf;
			xl = rtf;
			fl = f;
		}
		else {
			del = xh - rtf;
			xh = rtf;
			fh = f;
		}
		dx = xh - xl;
		if (fabs(del) < xacc || f == 0.0) return rtf;
	}
	return 0.0;
}
#undef MAXIT //cal1

#define MAXIT 30
float rtsec(float (*func)(float), float x1, float x2, float xacc)
{
	int j;
	float fl, f, dx, swap, xl, rts;

	fl = (*func)(x1);
	f = (*func)(x2);
	if (fabs(fl) < fabs(f)) {
		rts = x1;
		xl = x2;
		swap = fl;
		fl = f;
		f = swap;
	}
	else {
		xl = x1;
		rts = x2;
	}
	for (j = 1; j <= MAXIT; j++) {
		dx = (xl - rts) * f / (f - fl);
		xl = rts;
		fl = f;
		rts += dx;
		f = (*func)(rts);
		if (fabs(dx) < xacc || f == 0.0) return rts;
	}
	return 0.0;
}
#undef MAXIT //cal1

#define JMAX 20
float rtnewt(void (*funcd)(float, float*, float*), float x1, float x2,
	float xacc)
{
	int j;
	float df, dx, f, rtn;

	rtn = 0.5 * (x1 + x2);
	for (j = 1; j <= JMAX; j++) {
		(*funcd)(rtn, &f, &df);
		dx = f / df;
		rtn -= dx;
		if (fabs(dx) < xacc) return rtn;
	}
	return 0.0;
}
#undef JMAX  //cal2

#define MAXIT 100
float rtsafe(void (*funcd)(float, float*, float*), float x1, float x2,
	float xacc)
{
	int j;
	float df, dx, dxold, f, fh, fl;
	float temp, xh, xl, rts;

	(*funcd)(x1, &fl, &df);
	(*funcd)(x2, &fh, &df);

	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl = x1;
		xh = x2;
	}
	else {
		xh = x1;
		xl = x2;
	}
	rts = 0.5 * (x1 + x2);
	dxold = fabs(x2 - x1);
	dx = dxold;
	(*funcd)(rts, &f, &df);
	for (j = 1; j <= MAXIT; j++) {
		if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0)
			|| (fabs(2.0 * f) > fabs(dxold * df))) {
			dxold = dx;
			dx = 0.5 * (xh - xl);
			rts = xl + dx;
			if (xl == rts) return rts;
		}
		else {
			dxold = dx;
			dx = f / df;
			temp = rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		(*funcd)(rts, &f, &df);
		if (f < 0.0)
			xl = rts;
		else
			xh = rts;
	}

	return 0.0;
}
#undef MAXIT //cal2

void cal1(float(*func1) (float(*func) (float), float, float, float)  ,float (*bessel)(float),float* x1,float* x2,float xacc, int iteration) {
	//clock_t start = clock();
	for (int i = 1; i <= iteration; i++) {
		float root = (*func1)(bessel, x1[i], x2[i], xacc);
		printf("%0.8f\n", root);
	}
	//clock_t finish = clock();
	//printf("total time : %.8f\n", (float)(finish-start)/CLOCKS_PER_SEC);
}
void cal2(float(*func2) (float(*funcd) (float, float*, float*), float , float, float),void (*bessel)(float, float*, float*),float* x1,float* x2, float xacc,int iteration) {
	for (int i = 1; i <= iteration; i++) {
		float root = (*func2)(bessel, x1[i], x2[i], xacc);
		printf("%0.8f\n" ,root);
	}
}

void main() {
	//ranges
	double x1[200];
	double x2[200];

	//상대오차변화값
	float xxacc = 0.000001;

	//[1.0 , 10.0]
	float from = 1.0;
	float to = 10.0;

	int iteration = 100;
	int n = 1000;

	zbrak(bessj0, from, to, n, x1, x2, &iteration);
	//first, use bracketing routine (zbrak.c)

	printf("a) Bisection (rtbis.c)\n\n");
	cal1(rtbis, bessj0, x1, x2, xxacc, iteration);
	printf("\nb) Linear interpolation (rtflsp.c)\n\n");
	cal1(rtflsp, bessj0, x1, x2, xxacc, iteration);
	printf("\nc) Secant (rtsec.c)\n\n");
	cal1(rtsec, bessj0, x1, x2, xxacc, iteration);
	
	//for newton methods use the following property of the Bessel function:
	//dJ0(x)/dx = -J1(x)
	printf("\nd) Newton-Raphson (rtnewt.c)\n\n");
	cal2(rtnewt, derivative, x1, x2, xxacc, iteration);
	printf("\ne) Newton with bracketing (rtsafe.c)\n\n");
	cal2(rtsafe, derivative, x1, x2, xxacc, iteration);
	//rtsafe.c in NR => Newton(fast) + Bisection(convergence)
}