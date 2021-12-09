/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : ±è½ÂÈ¯
Created          : 26-03-2018
Modified         : 09-12-2021
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
	return sqrt(1 + pow((exp(x / 4491) - exp(-x / 4491)) / 2, 2));
}

// Define a function that defines the target equation (2 variable)
double myFunc2(const double x, const double t)
{
	double A = 1;
	double lamda = 2;
	double k = 2 * PI / lamda;
	double f = 3;
	double w = 2 * PI * f;

	return A * sin(k * x - w * t);
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

// Return the 1st partial differentiation results for the target equation
void partialGradient(double func(double x, double t), double x, double t, double dydx[], double dydt[])
{
	double h = 0.1;

	dydx[0] = (func(x + h, t) - func(x - h, t)) / (2 * h);
	dydt[0] = (func(x, t + h) - func(x, t - h)) / (2 * h);
}

// Return the 2nd partial differentiation results for the target equation
void partialGradient2(double func(double x, double t), double x, double t, double dy2dx2[], double dy2dt2[])
{
	double h = 0.1;

	dy2dx2[0] = (func(x - h, t) - 2 * func(x, t) + func(x + h, t)) / (h * h);
	dy2dt2[0] = (func(x, t - h) - 2 * func(x, t) + func(x, t + h)) / (h * h);
}

// Prints 1D array
void printVec(double* _vec, int _row)
{
	for (int i = 0; i < _row; i++)
		printf("Vector[%d] = %.1f \n", i, _vec[i]);

	printf("\n");
}

// Prints 1D array as range
void printVec_range(double* _vec, int r1, int r2)
{
	for (int i = r1; i <= r2; i++)
		printf("Vector[%d] = %f \n", i, _vec[i]);

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
	double m = 0.2;
	double c = 0.5;
	double k = 2 * 9.81;
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

// Return the Runge-Kutta 3rd order method results for the 2nd order target equation
void sys2RK3(void func(const double t, const double Y[], double dYdt[], int in_type), double y[], double z[], double t0, double tf, double h, double y_init, double z_init, int input_type)
{
	int n = (tf - t0) / h;
	y[0] = y_init;
	z[0] = z_init;

	double Y1[2] = { 0, };
	double dYdt1[2] = { 0, };
	double Y2[2] = { 0, };
	double dYdt2[2] = { 0, };
	double Y3[2] = { 0, };
	double dYdt3[2] = { 0, };

	for (int i = 0; i < n; i++)
	{
		Y1[0] = y[i];
		Y1[1] = z[i];
		mckfunc(t0 + i * h, Y1, dYdt1, input_type);
		Y2[0] = y[i] + 0.5 * dYdt1[0] * h;
		Y2[1] = z[i] + 0.5 * dYdt1[1] * h;
		mckfunc(t0 + (i + 0.5) * h, Y2, dYdt2, input_type);
		Y3[0] = y[i] + (-1 * dYdt1[0] + 2 * dYdt2[0]) * h;
		Y3[1] = z[i] + (-1 * dYdt1[1] + 2 * dYdt2[1]) * h;
		mckfunc(t0 + (i + 1) * h, Y3, dYdt3, input_type);
		y[i + 1] = y[i] + (dYdt1[0] + 4 * dYdt2[0] + dYdt3[0]) * h / 6;
		z[i + 1] = z[i] + (dYdt1[1] + 4 * dYdt2[1] + dYdt3[1]) * h / 6;
	}
}

// Matrix addition
Matrix	addMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols)
	{
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

// Matrix multiply
Matrix	multMat(Matrix _A, Matrix _B)
{
	if (_A.cols != _B.rows)
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'multMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = zeros(_A.rows, _B.cols);

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			for (int k = 0; k < _A.cols; k++)
				Out.at[i][j] += _A.at[i][k] * _B.at[k][j];

	return Out;
}

// Apply Gauss elimination
void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d)
{
	if (_A.rows != _A.cols)
	{
		printf("Warning! Square matrix is not used.\n");
		exit(0);
	}

	double mult = 0;
	copyVal(_A, _U);
	copyVal(_b, _d);

	for (int k = 0; k < _U.rows - 1; k++)
	{
		for (int i = k + 1; i < _U.rows; i++)
		{
			if (_U.at[k][k] == 0)
			{
				printf("Warning! Divided by zero.\n");
				exit(0);
			}

			mult = _U.at[i][k] / _U.at[k][k];
			_U.at[i][k] = 0;
			_d.at[i][0] = _d.at[i][0] - mult * _d.at[k][0];

			for (int j = k + 1; j < _U.cols; j++)
				_U.at[i][j] = _U.at[i][j] - mult * _U.at[k][j];
		}
	}
}

// Apply back-substitution after Gauss elimination
Matrix	backSub(Matrix _U, Matrix _d)
{
	Matrix Out = createMat(_d.rows, 1);
	double s = 0;

	for (int i = _U.rows - 1; i >= 0; i--)
	{
		s = 0;

		for (int j = i + 1; j < _U.cols; j++)
			s += _U.at[i][j] * Out.at[j][0];

		if (_U.at[i][i] == 0)
		{
			printf("Warning! Divided by zero.\n");
			exit(0);
		}

		Out.at[i][0] = (_d.at[i][0] - s) / _U.at[i][i];
	}

	return Out;
}

// Apply LU decomposition
void LUdecomp(Matrix A, Matrix L, Matrix U)
{
	if (A.rows != A.cols)
	{
		printf("Warning! Square matrix is not used.\n");
		exit(0);
	}

	double mult = 0;
	copyVal(A, U);

	for (int k = 0; k < U.rows - 1; k++)
	{
		for (int i = k + 1; i < U.rows; i++)
		{
			if (U.at[k][k] == 0)
			{
				printf("Warning! Divided by zero.\n");
				exit(0);
			}

			mult = U.at[i][k] / U.at[k][k];
			U.at[i][k] = 0;
			L.at[i][k] = mult;

			for (int j = k + 1; j < U.cols; j++)
				U.at[i][j] = U.at[i][j] - mult * U.at[k][j];
		}
	}
}

// Apply solve equation after LU decomposition
void solveLU(Matrix L, Matrix U, Matrix b, Matrix x)
{
	Matrix y = zeros(L.cols, b.cols);

	fwdsub(L, b, y);
	backsub(U, y, x);

	freeMat(y);
}

// Apply forward-substitution after LU decomposition
void fwdsub(Matrix L, Matrix b, Matrix y)
{
	double s = 0;

	for (int i = 0; i < L.rows; i++)
	{
		s = 0;

		for (int j = 0; j < i; j++)
			s += L.at[i][j] * y.at[j][0];

		if (L.at[i][i] == 0)
		{
			printf("Warning! Divided by zero.\n");
			exit(0);
		}

		y.at[i][0] = (b.at[i][0] - s) / L.at[i][i];
	}
}

// Apply back-substitution after LU decomposition
void backsub(Matrix U, Matrix y, Matrix x)
{
	double s = 0;

	for (int i = U.rows - 1; i >= 0; i--)
	{
		s = 0;

		for (int j = i + 1; j < U.cols; j++)
			s += U.at[i][j] * x.at[j][0];

		if (U.at[i][i] == 0)
		{
			printf("Warning! Divided by zero.\n");
			exit(0);
		}

		x.at[i][0] = (y.at[i][0] - s) / U.at[i][i];
	}
}

// Find inverse matrix
void inv(Matrix A, Matrix Ainv)
{
	Matrix U = zeros(A.rows, A.cols);
	Matrix L = eye(A.rows, A.cols);
	Matrix e = zeros(A.rows, 1);
	Matrix colA = zeros(Ainv.cols, 1);

	LUdecomp(A, L, U);

	if (U.at[U.rows - 1][U.cols - 1] == 0)
	{
		printf("Warning! Rank A does not equal to number of A columns.\n");
		exit(0);
	}

	for (int i = 0; i < Ainv.cols; i++)
	{
		initMat(e, 0);
		e.at[i][0] = 1;
		solveLU(L, U, e, colA);

		for (int j = 0; j < Ainv.rows; j++)
		{
			Ainv.at[j][i] = colA.at[j][0];
		}
	}

	freeMat(U);		freeMat(L);		freeMat(e);		freeMat(colA);
}

// Matrix and constant multiply
Matrix	multConst(Matrix _A, double c)
{
	if (c == 0)
	{
		printf("Warning! Divided by zero.\n");
		exit(0);
	}

	Matrix Out = createMat(_A.rows, _A.cols);

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			Out.at[i][j] = _A.at[i][j] * c;

	return Out;
}

// Return norm of vector
double	normVec(Matrix _A)
{
	if (_A.cols != 1)
	{
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'normVec' function");
		printf("\n*************************************************\n");
		exit(0);
	}

	double norm = 0;

	for (int i = 0; i < _A.rows; i++)
		norm += _A.at[i][0] * _A.at[i][0];

	return sqrt(norm);
}

// Apply QR decomposition
void QRdecomp(Matrix _A, Matrix _Q, Matrix _R)
{
	if (_A.rows != _A.cols)
	{
		printf("Warning! Square matrix is not used.\n");
		exit(0);
	}

	int n = _A.rows;
	Matrix I = eye(n, n);
	Matrix H = createMat(n, n);
	Matrix Q2 = copyMat(I);
	Matrix R2 = copyMat(_A);
	Matrix c = createMat(n, 1);
	Matrix e = zeros(n, 1);
	Matrix v = createMat(n, 1);
	Matrix vt = createMat(1, n);
	Matrix vvt = createMat(n, n);
	double vtv = 0;

	for (int k = 0; k < n - 1; k++)
	{
		for (int i = 0; i < n; i++)
			c.at[i][0] = R2.at[i][k];

		for (int i = 0; i < k; i++)
			c.at[i][0] = 0;

		initMat(e, 0);
		e.at[k][0] = 1;

		if (c.at[k][0] < 0)
			e.at[k][0] = -1;

		v = addMat(c, multConst(e, normVec(c)));

		vt = transpose(v);
		vvt = multMat(v, vt);
		vtv = normVec(v) * normVec(v);

		if (vtv == 0)
		{
			printf("Warning! Divided by zero.\n");
			exit(0);
		}

		H = addMat(I, multConst(vvt, -2 / vtv));

		Q2 = multMat(Q2, H);
		R2 = multMat(H, R2);
	}

	copyVal(R2, _R);
	copyVal(Q2, _Q);

	freeMat(I);		freeMat(H);		freeMat(Q2);	freeMat(R2);
	freeMat(c);		freeMat(e);		freeMat(v);		freeMat(vt);	freeMat(vvt);
}

// Return eigenvalues
Matrix eig(Matrix _A)
{
	if (_A.rows != _A.cols)
	{
		printf("Warning! Square matrix is not used.\n");
		exit(0);
	}

	Matrix Q = createMat(_A.rows, _A.cols);
	Matrix R = createMat(_A.rows, _A.cols);
	Matrix U = copyMat(_A);
	Matrix Out = createMat(U.rows, 1);

	for (int i = 0; i < 200; i++)
	{
		QRdecomp(U, Q, R);
		U = multMat(R, Q);
	}

	for (int i = 0; i < U.rows; i++)
		Out.at[i][0] = U.at[i][i];

	freeMat(Q);		freeMat(R);		freeMat(U);

	return Out;
}

// Return eigenvectors
Matrix eigvec(Matrix _A)
{
	if (_A.rows != _A.cols)
	{
		printf("Warning! Square matrix is not used.\n");
		exit(0);
	}

	int n = _A.rows;
	Matrix I = eye(n, n);
	Matrix B = createMat(n, n);
	Matrix b1 = createMat(n - 1, n);
	Matrix b2 = createMat(n - 1, n - 1);
	Matrix b2_inv = createMat(n - 1, n - 1);
	Matrix c1 = createMat(n - 1, 1);
	Matrix c2 = createMat(n - 1, 1);
	double normOut = 0;
	Matrix Out = createMat(n, n);

	Matrix lamda = eig(_A);

	for (int i = 0; i < n; i++)
	{
		B = addMat(_A, multConst(I, -lamda.at[i][0]));

		for (int j = 0; j < n; j++)
		{
			if (i < j)
			{
				for (int k = 0; k < n; k++)
					b1.at[j - 1][k] = B.at[j][k];
			}
			else if (i > j)
			{
				for (int k = 0; k < n; k++)
					b1.at[j][k] = B.at[j][k];
			}
		}

		for (int j = 0; j < n; j++)
		{
			if (i < j)
			{
				for (int k = 0; k < n - 1; k++)
					b2.at[k][j - 1] = b1.at[k][j];
			}
			else if (i > j)
			{
				for (int k = 0; k < n - 1; k++)
					b2.at[k][j] = b1.at[k][j];
			}
			else
			{
				for (int k = 0; k < n - 1; k++)
					c1.at[k][0] = -b1.at[k][j];
			}
		}

		inv(b2, b2_inv);
		c2 = multMat(b2_inv, c1);

		for (int j = 0; j < n; j++)
		{
			if (i < j)
			{
				Out.at[j][i] = c2.at[j - 1][0];
			}
			else if (i > j)
			{
				Out.at[j][i] = c2.at[j][0];
			}
			else
			{
				Out.at[j][i] = 1;
			}
		}

		normOut = 0;
		for (int j = 0; j < n; j++)
			normOut += Out.at[j][i] * Out.at[j][i];

		if (normOut == 0)
		{
			printf("Warning! Divided by zero.\n");
			exit(0);
		}

		normOut = sqrt(normOut);
		for (int j = 0; j < n; j++)
			Out.at[j][i] /= normOut;
	}

	freeMat(I);		freeMat(B);		freeMat(b1);	freeMat(b2);
	freeMat(b2_inv);	freeMat(c1);	freeMat(c2);

	return Out;
}

// Apply linear regression
void linearRegression(Matrix x, Matrix y, Matrix z)
{
	if (x.rows != y.rows)
	{
		printf("\n***********************************************************");
		printf("\n  ERROR: The number of elements of x and y are not same.");
		printf("\n***********************************************************\n");
		exit(0);
	}

	int m = x.rows;
	double Sx = 0;
	double Sy = 0;
	double Sxx = 0;
	double Sxy = 0;

	for (int i = 0; i < m; i++)
	{
		Sx += x.at[i][0];
		Sy += y.at[i][0];
		Sxx += x.at[i][0] * x.at[i][0];
		Sxy += x.at[i][0] * y.at[i][0];
	}

	z.at[0][0] = (m * Sxy - Sx * Sy) / (m * Sxx - Sx * Sx);
	z.at[1][0] = (Sxx * Sy - Sx * Sxy) / (m * Sxx - Sx * Sx);
}

// Apply ploynomial curve fit
void polyfit(Matrix x, Matrix y, Matrix z, int n)
{
	if (x.rows != y.rows)
	{
		printf("\n***********************************************************");
		printf("\n  ERROR: The number of elements of x and y are not same.");
		printf("\n***********************************************************\n");
		exit(0);
	}

	int m = x.rows;
	Matrix Sx = zeros(2 * n + 1, 1);
	Matrix S = createMat(n + 1, n + 1);
	Matrix Sinv = createMat(n + 1, n + 1);
	Matrix b = zeros(n + 1, 1);

	for (int i = 0; i < 2 * n + 1; i++)
		for (int k = 0; k < m; k++)
			Sx.at[i][0] += pow(x.at[k][0], i);

	for (int i = 0; i < n + 1; i++)
		for (int j = 0; j < n + 1; j++)
			S.at[i][j] = Sx.at[2 * n - (i + j)][0];

	for (int i = 0; i < n + 1; i++)
		for (int k = 0; k < m; k++)
			b.at[i][0] += y.at[k][0] * pow(x.at[k][0], n - i);

	inv(S, Sinv);
	Matrix ztemp = multMat(Sinv, b);
	copyVal(ztemp, z);

	freeMat(Sx);		freeMat(S);			freeMat(Sinv);		freeMat(b);
	freeMat(ztemp);
}

// Define a original function
double myfunc(const double x)
{
	double w0 = 20000;
	double L = 4;
	double E = 7 * pow(10, 10);
	double I = 52.9 * pow(10, -6);

	return w0 * (7 * pow(L, 4) * x - 10 * L * L * pow(x, 3) + 3 * pow(x, 5)) / (360 * L * E * I);
}

// Define a function F for Newton-Raphson method
double func(const double x)
{
	double w0 = 20000;
	double L = 4;
	double E = 7 * pow(10, 10);
	double I = 52.9 * pow(10, -6);

	return w0 * (7 * pow(L, 4) - 30 * L * L * x * x + 15 * pow(x, 4)) / (360 * L * E * I);
}

// Define a function dF for Newton-Raphson method
double dfunc(const double x)
{
	double w0 = 20000;
	double L = 4;
	double E = 7 * pow(10, 10);
	double I = 52.9 * pow(10, -6);

	return w0 * (-60 * L * L * x + 60 * pow(x, 3)) / (360 * L * E * I);
}

// Apply Newton-Raphson method
double newtonRaphson(double func(const double x), double dfunc(const double x), double x0, double tol)
{
	int Nmax = 1000;
	double xn = x0;
	double ep = fabs(func(xn));;

	for (int i = 0; i < Nmax; i++)
	{
		printf("Iteration:%d \t", i);
		printf("X(n): %f \t", xn);
		printf("Tolerance: %.10f\n", ep);

		if (ep <= tol)
			break;

		if (dfunc(xn) == 0)
		{
			printf("\nWarning! Divided by zero.\n");
			exit(0);
		}

		xn -= func(xn) / dfunc(xn);
		ep = fabs(func(xn));
	}

	return xn;
}