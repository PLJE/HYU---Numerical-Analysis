#include<stdio.h>
#define _USE_MATH_DEFINES
#include<math.h>

//HW3 - 2 
//주어진 세개 문제 + 원하는 방정식1개 rtsafe로 루트 구하기

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

void cal2(float(*func2) (float(*funcd) (float, float*, float*), float, float, float), void (*bessel)(float, float*, float*), float* x1, float* x2, float xacc, int iteration , int prob) {
	for (int i = 0; i < iteration; i++) {
		float root = (*func2)(bessel, x1[i], x2[i], xacc);
		printf("problem %d : %0.8f\n", prob, root);
	}
}

//problem 1
float fx1(float x) {
	return 10 * exp(-x) * sin(2 * M_PI * x) - 2;
}
void derivative1(float x, float* fn, float* df) {
	*fn = fx1(x);
	*df = -10 * exp(-x) * sin(2 * M_PI * x) + 20 * M_PI * exp(-x) * cos(2 * M_PI * x);
}

//problem 2
float fx2(float x) {
	return pow(x,2) - 2*x * exp(-x) + exp(-2*x);
}
void derivative2(float x, float* fn, float* df) {
	*fn = fx2(x);
	*df = 2 * x - 2 * exp(-x) + 2*x * exp(-x) - 2 * exp((-2)*x);
}

//problem 3
float fx3(float x) {
	return cos(x + sqrt(2)) + x * (x / 2 + sqrt(2));
}
void derivative3(float x, float* fn, float* df) {
	*fn = fx3(x);
	*df = -sin(x + sqrt(2)) + x + sqrt(2);
}

//problem 4
float fx4(float x) {
	return pow(x, 2) + sin(-3 * x);
}
void derivative4(float x, float* fn, float* df) {
	*fn = fx4(x);
	*df = 2 * x - 3 * cos(-3 * x);
}

void main() {
	
	float xacc = 0.000001;

	float p1_x1[1] = {  0.1};
	float p1_x2[1] = {  1.0};

	float p2_x1[1] = { 0 };
	float p2_x2[1] = { 1.0};

	float p3_x1[1] = { -2.0 };
	float p3_x2[1] = { -1.0 };

	float p4_x1[1] = { 1.0};
	float p4_x2[1] = { 2.0 };

	cal2(rtsafe, derivative1, p1_x1, p1_x2, xacc, 1, 1);
	cal2(rtsafe, derivative2, p2_x1, p2_x2, xacc, 1, 2);
	cal2(rtsafe, derivative3, p3_x1, p3_x2, xacc, 1, 3);
	cal2(rtsafe, derivative4, p4_x1, p4_x2, xacc, 1, 4);
}