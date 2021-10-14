/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : ±è½ÂÈ¯
Created          : 26-03-2018
Modified         : 14-10-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/


#include "myNM.h"


// define factorial function
double factorial(double _x)
{
	if (_x <= 1)
		return 1;
	else
		return _x * factorial(_x - 1);
}

// define Taylor sine function[rad]
double sinTaylor(double x)
{
	int N_max = 20;
	double epsilon = 1e-5;

	double S_N = 0, S_N_prev = 0, rel_chg = 0;
	int N = 0;

	do
	{
		N++;
		S_N_prev = S_N;

		S_N = 0;
		for (int k = 0; k <= N - 1; k++)
		{
			S_N += pow(-1, k) * pow(x, 2 * (double)k + 1) / factorial(2 * (double)k + 1);
		}

		rel_chg = fabs((S_N - S_N_prev) / S_N_prev);

	} while (N < N_max && rel_chg >= epsilon);

	CompareSinAns(x, N, S_N); // call compare function[rad]

	return 0;
}

// define Taylor sine function[deg]
double sindTaylor(double d)
{
	double d_to_r = d * PI / 180; // convert sine input from deg to rad

	int N_max = 20;
	double epsilon = 1e-5;

	double S_N_d = 0, S_N_prev_d = 0, rel_chg_d = 0;
	int N_d = 0;

	do
	{
		N_d++;
		S_N_prev_d = S_N_d;

		S_N_d = 0;
		for (int k = 0; k <= N_d - 1; k++)
		{
			S_N_d += pow(-1, k) * pow(d_to_r, 2 * k + 1) / factorial(2 * k + 1);
		}

		rel_chg_d = fabs((S_N_d - S_N_prev_d) / S_N_prev_d);

	} while (N_d < N_max && rel_chg_d >= epsilon);

	CompareSindAns(d, N_d, S_N_d); // call compare function[deg]

	return 0;
}

// define compare function[rad]
void CompareSinAns(double sin_in, int iter, double approx)
{
	printf("\n\n");
	printf("=======================================\n");
	printf("    sin( %f[rad] ) Calculation   \n", sin_in);
	printf("=======================================\n");
	printf("   -  iteration No. = %d        \n", iter);
	printf("   -  My     result = %3.12f    \n", approx);
	printf("   -  Math.h result = %3.12f    \n", sin(sin_in));
	printf("   -  absolute err. = %3.12f    \n", approx - sin(sin_in));
	printf("=======================================\n");
}

// define compare function[deg]
void CompareSindAns(double sin_in_d, int iter_d, double approx_d)
{
	printf("\n\n");
	printf("=======================================\n");
	printf("    sin( %f[deg] ) Calculation   \n", sin_in_d);
	printf("=======================================\n");
	printf("   -  iteration No. = %d        \n", iter_d);
	printf("   -  My     result = %3.12f    \n", approx_d);
	printf("   -  Math.h result = %3.12f    \n", sin(sin_in_d * PI / 180));
	printf("   -  absolute err. = %3.12f    \n", approx_d - sin(sin_in_d * PI / 180));
	printf("=======================================\n");
}

// Define a function that defines the target equation (1 variable)
double myFunc(const double x)
{
	return sqrt(1 - (x * x));
}

// Define a function that defines the target equation (2 variable)
double myFunc2(const double t, const double y)
{
	double tau = 1;
	double T = 1 / tau;
	double f = 10;
	double Vm = 1;
	double w = 2 * PI * f;

	return -1 * T * y + T * Vm * cos(w * t);
}

// Return the dy/dx results for the input data (truncation error: O(h^2))
void gradient1D(double x[], double y[], double dydx[], int m)
{
	double h = x[1] - x[0];

	// Check for length of x and y. They must be equal.
	if (sizeof(x) != sizeof(y))
	{
		printf("WARNING!: Length of X and Y are not equal\n");
		return;
	}

	// Check length of dataset. It must be larger than 2.
	if (m < 3)
	{
		printf("WARNING!: Length of dataset is short\n");
		return;
	}

	// Three-Point Forward Difference
	dydx[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);

	// Two-Point Central Difference
	for (int k = 1; k < m - 1; k++)
		dydx[k] = (y[k + 1] - y[k - 1]) / (2 * h);

	// Three-Point Forward Difference
	dydx[m - 1] = (y[m - 3] - 4 * y[m - 2] + 3 * y[m - 1]) / (2 * h);
}

// Function callback
void func_call(double func(const double x), double xin)
{
	double yout = func(xin);
	printf("Y_out by my_func=%f\n", yout);
}

// Return the dy/dx results for the target equation (truncation error: O(h^2))
void gradientFunc(double func(const double x), double x[], double dydx[], int m)
{
	double* y;

	y = (double*)malloc(sizeof(double) * m);

	for (int i = 0; i < m; i++)
	{
		y[i] = func(x[i]);
	}

	// get gradient 1D
	gradient1D(x, y, dydx, m);

	free(y);
}

// Prints 1D array
void printVec(double* _vec, int _row)
{
	for (int i = 0; i < _row; i++)
		printf("Vector[%d] = %.1f \n", i, _vec[i]);

	printf("\n");
}

// Integration using rectangular method for discrete data inputs
double IntegrateRect(double x[], double y[], int m)
{
	int N = m - 1;
	double I = 0;

	for (int i = 0; i < N; i++)
		I += y[i] * (x[i + 1] - x[i]);

	return I;
}

// Integration using trapezoidal method for discrete data inputs
double trapz(double x[], double y[], int m)
{
	int N = m - 1;
	double I = 0;

	for (int i = 0; i < N; i++)
		I += (y[i] + y[i + 1]) * (x[i + 1] - x[i]);

	return I / 2;
}

// Integration using Simpson1/3 method for function inputs
double integral(double func(const double x), double a, double b, int n)
{
	double h = (b - a) / n;
	double I = func(a) + 4 * func(b - h) + func(b);

	for (int i = 1; i < n - 2; i += 2)
	{
		double x = a + i * h;
		I += 4 * func(x) + 2 * func(x + h);
	}

	return I * h / 3;
}

// Integration using Simpson3/8 method for function inputs
double simpson38(double func(const double x), double a, double b, int n)
{
	double h = (b - a) / n;
	double I = func(a) + 3 * (func(b - 2 * h) + func(b - h)) + func(b);

	for (int i = 1; i < n - 4; i += 3)
	{
		double x = a + i * h;
		I += 3 * (func(x) + func(x + h)) + 2 * func(x + 2 * h);
	}

	return I * h * 3 / 8;
}

// Return the Euler method results for the target equation
void odeEU(double myfunc(const double t, const double y), double y[], double t0, double tf, double h)
{
	int m = (tf - t0) / h;

	for (int i = 0; i < m; i++)
	{
		y[i + 1] = y[i] + myfunc(t0 + i * h, y[i]) * h;
	}
}

// Return the Euler's Modified method results for the target equation
void odeEM(double myfunc(const double t, const double y), double y[], double t0, double tf, double h)
{
	int m = (tf - t0) / h;
	double s1 = 0;
	double s2 = 0;

	for (int i = 0; i < m; i++)
	{
		s1 = myfunc(t0 + i * h, y[i]);
		y[i + 1] = y[i] + s1 * h;
		s2 = myfunc(t0 + (i + 1) * h, y[i + 1]);
		y[i + 1] = y[i] + (s1 + s2) * h / 2;
	}
}

// Return the Runge-Kutta 2nd order method results for the target equation
void odeRK2(double myfunc(const double t, const double y), double y[], double t0, double tf, double h)
{
	int m = (tf - t0) / h;
	double alpha = 1;
	double C2 = 1 / (2 * alpha);
	double C1 = 1 - C2;
	double K1 = 0;
	double K2 = 0;

	for (int i = 0; i < m; i++)
	{
		K1 = myfunc(t0 + i * h, y[i]);
		K2 = myfunc(t0 + (i + alpha) * h, y[i] + alpha * K1 * h);
		y[i + 1] = y[i] + (C1 * K1 + C2 * K2) * h;
	}
}

// Return the Runge-Kutta 3rd order method results for the target equation
void odeRK3(double myfunc(const double t, const double y), double y[], double t0, double tf, double h)
{
	int m = (tf - t0) / h;
	double K1 = 0;
	double K2 = 0;
	double K3 = 0;

	for (int i = 0; i < m; i++)
	{
		K1 = myfunc(t0 + i * h, y[i]);
		K2 = myfunc(t0 + (i + 0.5) * h, y[i] + 0.5 * K1 * h);
		K3 = myfunc(t0 + (i + 1) * h, y[i] + (-1 * K1 + 2 * K2) * h);
		y[i + 1] = y[i] + (K1 + 4 * K2 + K3) * h / 6;
	}
}

// Return the called method results for the target equation
void ode(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, int method)
{
	switch (method)
	{
	case EU:
		odeEU(myfunc, y, t0, tf, h);
		break;
	case EM:
		odeEM(myfunc, y, t0, tf, h);
		break;
	case RK2:
		odeRK2(myfunc, y, t0, tf, h);
		break;
	case RK3:
		odeRK3(myfunc, y, t0, tf, h);
		break;
	default:
		printf("\nWarning! That method does not exist.\n");
		system("pause");
		exit(0);
	}
}

// Governing equation and return vector differentiation results for the target vector
void mckfunc(const double t, const double Y[], double dYdt[], int in_type)
{
	double m = 1;
	double c = 7;
	double k = 6.9;
	double A = 2;
	double f = 5;

	double F_HR = A * cos(2 * PI * f * t);
	double F_FV = 0;
	double F_SR = A;

	dYdt[0] = Y[1];

	switch (in_type)
	{
	case HR:
		dYdt[1] = (-1 * k * Y[0] - 1 * c * Y[1] + F_HR) / m;
		break;
	case FV:
		dYdt[1] = (-1 * k * Y[0] - 1 * c * Y[1] + F_FV) / m;
		break;
	case SR:
		dYdt[1] = (-1 * k * Y[0] - 1 * c * Y[1] + F_SR) / m;
		break;
	default:
		printf("\nWarning! That input type does not exist.\n");
		system("pause");
		exit(0);
	}
}

// Return the Runge-Kutta 2nd order method results for the 2nd order target equation
void sys2RK2(void func(const double t, const double Y[], double dYdt[], int in_type), double y[], double z[], double t0, double tf, double h, double y_init, double z_init, int input_type)
{
	int n = (tf - t0) / h;
	y[0] = y_init;
	z[0] = z_init;

	double Y1[2] = { 0, };
	double dYdt1[2] = { 0, };
	double Y2[2] = { 0, };
	double dYdt2[2] = { 0, };

	for (int i = 0; i < n; i++)
	{
		Y1[0] = y[i];
		Y1[1] = z[i];
		mckfunc(t0 + i * h, Y1, dYdt1, input_type);
		Y2[0] = y[i] + dYdt1[0] * h;
		Y2[1] = z[i] + dYdt1[1] * h;
		mckfunc(t0 + (i + 1) * h, Y2, dYdt2, input_type);
		y[i + 1] = y[i] + (dYdt1[0] + dYdt2[0]) * h / 2;
		z[i + 1] = z[i] + (dYdt1[1] + dYdt2[1]) * h / 2;
	}
}

// Matrix addition
Matrix	addMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Out;
}

// Apply back-substitution
Matrix	backSub(Matrix _U, Matrix _b)
{
	Matrix Out = createMat(_b.rows, 1);

	for (int i = _U.rows; i > 0; i--) {
		// add your code here
		for (int j = i + 1; j <= _U.cols; j++)
			;// add your code here
		// add your code here
	}

	return Out;
}