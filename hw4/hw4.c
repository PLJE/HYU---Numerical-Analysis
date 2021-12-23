#include<stdio.h>
#define _USE_MATH_DEFINES
#include<math.h>

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
		if (fabs(del) < xacc || f == 0.0) { 
			printf("iteration : %d\n", j - 1);
			return rtf; 
		}
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
		if (fabs(dx) < xacc || f == 0.0) {
			printf("iteration : %d\n", j - 1);
			return rts;
		}
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
		if (fabs(dx) < xacc) {
			printf("iteration : %d\n", j - 1);
			return rtn;
		}
	}
	return 0.0;
}
#undef JMAX  //cal2

#define JMAX 40
float rtbis(float (*func)(float), float x1, float x2, float xacc )
{ //bisection rtbis
	int j;
	float dx, f, fmid, xmid, rtb;

	f = (*func)(x1);
	fmid = (*func)(x2);
	rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
	for (j = 1; j <= JMAX; j++) {
		fmid = (*func)(xmid = rtb + (dx *= 0.5));
		if (fmid <= 0.0) rtb = xmid;
		if (fabs(dx) < xacc || fmid == 0.0) { 
			printf("iteration : %d\n", j-1);
			return rtb; 
		}
	}
	return 0.0;
}
#undef JMAX //cal1


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
			if (xl == rts) {
				printf("iteration : %d\n", j - 1);
				return rts;
			}
		}
		else {
			dxold = dx;
			dx = f / df;
			temp = rts;
			rts -= dx;
			if (temp == rts) {
				printf("iteration : %d\n", j - 1);
				return rts;
			}
		}
		if (fabs(dx) < xacc) {
			printf("iteration : %d\n", j - 1);
			return rts;
		}
		(*funcd)(rts, &f, &df);
		if (f < 0.0)
			xl = rts;
		else
			xh = rts;
	}

	return 0.0;
}
#undef MAXIT //cal2

float fr(float r) {
	return exp(-0.005 * r) * cos(sqrt(2000 - 0.01 * r * r) * 0.05) - 0.01;
}
void dfr(float r, float* fn, float* df) {
	*fn = fr(r);
	float cal = sqrt(2000 - r * r * 0.01) / 20;
	*df = r * exp(-r / 200) * sin(cal) / (2000 * cal * 20) - exp(-r / 200) * cos(cal) / 200;
}

float f832(float x) {
	//return pow(x * x + .9 * .9, 1.5) / x - 100 / (8.85 * M_PI);
	
	return (100 /( 8.85 * M_PI)) * (x / ( pow((x * x + 0.9 * 0.9), 1.5)) ) - 1;
}
void df832(float x, float* fn, float* df) {
	*fn = f832(x);
	*df = -(20 * M_PI * (200 * x * x - 81)) / (177 * pow((x * x + 81 / 100), 2.5));

	//*fn = f832(x);
	//float tmp = x * x + .9 * .9;
	//*df = 3 * pow(tmp, 0.5) - pow(tmp, 1.5) / x / x;
}

float f836(float x) {
	return 0.99403 + (1.671 * pow(10, -4) * x) + (9.7215 * pow(10, -8) * x * x)
		-(9.5838 * pow(10, -11) * x * x * x) + 
		(1.9520 * pow(10, -14) * x * x * x * x) - 1.2;
}
void df836(float x, float* fn, float* df) {
	*fn = f836(x);
	*df = (976 * x * x * x - 3593925 * x * x + 2430375000
		* x + 2088750000000) / 12500000000000000;
}

void pr832() {
	float root;
	printf("------------------problem 8.32-------------------------------\n");
	printf("a) bisection (rtbis.c)\n");
	printf("\txacc : 10^-4\n");
	root = (*rtbis)(f832, 0.2, 0.3, 0.0001);
	printf("root1:%.8f\n", root);
	root = (*rtbis)(f832, 1.5, 1.55, 0.0001);
	printf("root2:%.8f\n", root);
	printf("\txacc : 10^-6\n");
	root = (*rtbis)(f832, 0.2, 0.3, 0.000001);
	printf("root1:%.8f\n", root);
	root = (*rtbis)(f832, 1.5, 1.55, 0.000001);
	printf("root2:%.8f\n\n", root);

	printf("b) linear interpolation (rtflsp.c)\n");
	printf("\txacc : 10^-4\n");
	root = (*rtflsp)(f832, 0.2, 0.3, 0.0001);
	printf("root1:%.8f\n", root);
	root = (*rtflsp)(f832, 1.5, 1.55, 0.0001);
	printf("root2:%.8f\n", root);
	printf("\txacc : 10^-6\n");
	root = (*rtflsp)(f832, 0.2, 0.3, 0.000001);
	printf("root1:%.8f\n", root);
	root = (*rtflsp)(f832, 1.5, 1.55, 0.000001);
	printf("root2:%.8f\n\n", root);

	printf("c) Secant (rtsec.c)\n");
	printf("\txacc : 10^-4\n");
	root = (*rtsec)(f832, 0.2, 0.3, 0.0001);
	printf("root1:%.8f\n", root);
	root = (*rtsec)(f832, 1.5, 1.55, 0.0001);
	printf("root2:%.8f\n", root);
	printf("\txacc : 10^-6\n");
	root = (*rtsec)(f832, 0.2, 0.3, 0.000001);
	printf("root1:%.8f\n", root);
	root = (*rtsec)(f832, 1.5, 1.55, 0.000001);
	printf("root2:%.8f\n\n", root);

	printf("d) Newton-Raphson (rtnewt.c)\n");
	printf("\txacc : 10^-4\n");
	root = (*rtnewt)(df832, 0.2, 0.3, 0.0001);
	printf("root1:%.8f\n", root);
	root = (*rtnewt)(df832, 1.5, 1.54, 0.0001);
	printf("root2:%.8f\n", root);
	printf("\txacc : 10^-6\n");
	root = (*rtnewt)(df832, 0.2, 0.25, 0.000001);
	printf("root1:%.8f\n", root);
	root = (*rtnewt)(df832, 1.5, 1.55, 0.000001);
	printf("root2:%.8f\n\n", root);

	printf("e) rtsafe (rtsafe.c)\n");
	printf("\txacc : 10^-4\n");
	root = (*rtsafe)(df832, 0.2, 0.3, 0.0001);
	printf("root1:%.8f\n", root);
	root = (*rtsafe)(df832, 1.5, 1.54, 0.0001);
	printf("root2:%.8f\n", root);
	printf("\txacc : 10^-6\n");
	root = (*rtsafe)(df832, 0.2, 0.3, 0.000001);
	printf("root1:%.8f\n", root);
	root = (*rtsafe)(df832, 1.5, 1.55, 0.000001);
	printf("root2:%.8f\n\n", root);
}

void pr836() {
	float root;
	printf("------------------problem 8.36-------------------------------\n");
	printf("a) bisection (rtbis.c)\n");
	root = (*rtbis)(f836, 1000, 2000, 0.0001);
	printf("xacc=0.0001 root:%.8f\n", root);
	root = (*rtbis)(f836, 1000, 2000, 0.000001);
	printf("xacc=0.000001 root:%.8f\n\n", root);

	printf("b) linear interpolation (rtflsp.c)\n");
	root = (*rtflsp)(f836, 1000, 2000, 0.0001);
	printf("xacc=0.0001 root:%.8f\n", root);
	root = (*rtflsp)(f836, 1000, 2000, 0.000001);
	printf("xacc=0.000001 root:%.8f\n\n", root);

	printf("c) Secant (rtsec.c)\n");
	root = (*rtsec)(f836, 1000, 2000, 0.0001);
	printf("xacc=0.0001 root:%.8f\n", root);
	root = (*rtsec)(f836, 1000, 2000, 0.000001);
	printf("xacc=0.000001 root:%.8f\n\n", root);

	printf("d) Newton-Raphson (rtnewt.c)\n");
	root = (*rtnewt)(df836, 1000, 2000, 0.0001);
	printf("xacc=0.0001 root:%.8f\n", root);
	root = (*rtnewt)(df836, 1000, 1150, 0.000001);
	printf("xacc=0.000001 root:%.8f\n\n", root);

	printf("e) rtsafe (rtsafe.c)\n");
	root = (*rtsafe)(df836, 1000, 2000, 0.0001);
	printf("xacc=0.0001 root:%.8f\n", root);
	root = (*rtsafe)(df836, 1000, 2000, 0.000001);
	printf("xacc=0.000001 root:%.8f\n\n", root);
}

void main() {
	float root;
	int cnt;

	printf("a) bisection (rtbis.c)\n");
	root = (*rtbis)(fr, 0, 400, 0.0001);
	printf("xacc=0.0001 root:%.8f\n", root);
	root = (*rtbis)(fr, 0, 400, 0.000001);
	printf("xacc=0.000001 root:%.8f\n\n", root);

	printf("b) linear interpolation (rtflsp.c)\n");
	root = (*rtflsp)(fr, 0, 400, 0.0001);
	printf("xacc=0.0001 root:%.8f\n", root);
	root = (*rtflsp)(fr, 0, 400, 0.000001);
	printf("xacc=0.000001 root:%.8f\n\n", root);

	printf("c) Secant (rtsec.c)\n");
	root = (*rtsec)(fr, 0, 400, 0.0001);
	printf("xacc=0.0001 root:%.8f\n", root);
	root = (*rtsec)(fr, 0, 400, 0.000001);
	printf("xacc=0.000001 root:%.8f\n\n", root);

	printf("d) Newton-Raphson (rtnewt.c)\n");
	root = (*rtnewt)(dfr, 0, 400, 0.0001);
	printf("xacc=0.0001 root:%.8f\n", root);
	root = (*rtnewt)(dfr, 0, 400, 0.000001);
	printf("xacc=0.000001 root:%.8f\n\n", root);

	printf("e) rtsafe (rtsafe.c)\n");
	root = (*rtsafe)(dfr, 0, 400, 0.0001);
	printf("xacc=0.0001 root:%.8f\n", root);
	root = (*rtsafe)(dfr, 0, 400, 0.000001);
	printf("xacc=0.000001 root:%.8f\n\n", root);

	pr832();
	pr836();
}
