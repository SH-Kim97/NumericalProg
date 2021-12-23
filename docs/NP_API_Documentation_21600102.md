---
Editor Name: 김승환
Modified Date: 23-12-2021
Description: Numerical Method
---


# Reference of Numerical Programming API

`#include "myNM.h"`

`#include "myMatrix.h"`



## Taylor of Sine Function

### sinTaylor()

Print taylor sine value[rad]

```text
double sinTaylor(double x);
```

**Parameters**

* **x**: angle in radian

**Example code**

```text
double x = PI / 3;

sinTaylor(x);
```



### sindTaylor()

Print taylor sine value[deg]

```text
double sindTaylor(double d);
```

**Parameters**

* **d**: angle in degree

**Example code**

```text
double d = 60;

sindTaylor(d);
```



## Differentiation

### gradient1D()

Returns the dy/dx results for the input data of y(x) with truncation error: O(h^2)

```text
void gradient1D(double x[], double y[], double dydx[], int m);
```

**Parameters**

* **x[]**: data set of x
* **y[]**: data set of y
* **dydx[]**: result of gradient
* **m:** length of **x[]** and **y[]**

**Example code**

```text
int m = 21;
double t[21] = { 0 };
for (int i = 0; i < m; i++)
	t[i] = 0.2 * i;

double x[] = { -5.87, -4.23, -2.55, -0.89, 0.67, 2.09, 3.31, 4.31, 5.06, 5.55, 5.78, 5.77, 5.52, 5.08, 4.46, 3.72, 2.88, 2.00, 1.10, 0.23, -0.59 };
double dxdt[21] = { 0 };

gradient1D(t, x, dxdt, m);
```



### gradientFunc()

Return the dy/dx results for the target equation (truncation error: O(h^2))

```text
void gradientFunc(double func(const double x), double x[], double dydx[], int m);
```

**Parameters**

* **func**: Function **func** is defined.
* **x[]**: data set of x
* **dydx[]**: result of gradient
* **m:** length of **x[]**

**Example code**

```text
int m = 21;
double dydx[21] = { 0 };
double t[21] = { 0 };
for (int i = 0; i < m; i++)
	t[i] = 0.2 * i;

gradientFunc(myFunc, t, dydx, m);
```



### partialGradient()

Return the 1st partial differentiation results for the target equation

```text
void partialGradient(double func(double x, double t), double x, double t, double dydx[], double dydt[]);
```

**Parameters**

* **func**: Function **func** is defined.
* **x**: x value
* **t**: t value
* **dydx[]**: result of partial gradient about x
* **dydt[]**: result of partial gradient about t

**Example code**

```text
double dydx[1] = { 0 };
double dydt[1] = { 0 };

double x3 = 0.2;
double t3 = 1;

partialGradient(myFunc2, x3, t3, dydx, dydt);
```



### partialGradient2()

Return the 2nd partial differentiation results for the target equation

```text
void partialGradient2(double func(double x, double t), double x, double t, double dy2dx2[], double dy2dt2[]);
```

**Parameters**

* **func**: Function **func** is defined.
* **x**: x value
* **t**: t value
* **dy2dx2[]**: result of 2nd partial gradient about x
* **dy2dt2[]**: result of 2nd partial gradient about t

**Example code**

```text
double dy2dx2[1] = { 0 };
double dy2dt2[1] = { 0 };

double x3 = 0.2;
double t3 = 1;

partialGradient2(myFunc2, x3, t3, dy2dx2, dy2dt2);
```



## Integration

### IntegrateRect\(\)

Integration using rectangular method for discrete data inputs

```text
double IntegrateRect(double x[], double y[], int m);
```

**Parameters**

* **x[]**: data set of x
* **y[]**: data set of y
* **m**: length of **x[]** and **y[]**

**Example code**

```text
double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
int M = sizeof(x) / sizeof(x[0]);

double I_rect = IntegrateRect(x, y, M);
```



### trapz\(\)

Integration using trapezoidal method for discrete data inputs

```text
double trapz(double x[], double y[], int m);
```

**Parameters**

* **x[]**: data set of x
* **y[]**: data set of y
* **m**: length of **x[]** and **y[]**

**Example code**

```text
double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
int M = sizeof(x) / sizeof(x[0]);

double I_trapz = trapz(x, y, M);
```



### integral\(\)

Integration using Simpson1/3 method for function inputs

```text
double integral(double func(const double x), double a, double b, int n);
```

**Parameters**

* **func**: Function **func** is defined.
* **a**: start point of x
* **b**: end point of x
* **n**: number of interval

**Example code**

```text
int a = -1, b = 1, N = 12;

double I_simpson13 = integral(myFunc, a, b, N);
```



### simpson38\(\)

Integration using Simpson3/8 method for function inputs

```text
double simpson38(double func(const double x), double a, double b, int n);
```

**Parameters**

* **func**: Function **func** is defined.
* **a**: starting point
* **b**: ending point
* **n**: number of interval

**Example code**

```text
int a = -1, b = 1, N = 12;

double I_simpson38 = simpson38(myFunc, a, b, N);
```



## ODE-IVP

### odeEU\(\)

Return the Euler method results for the target equation

```text
void odeEU(double myfunc(const double t, const double y), double y[], double t0, double tf, double h);
```

**Parameters**

* **myfunc**: Function **myfunc** is defined.
* **y\[\]**: solution of ODE
* **t0**: starting point
* **tf**: ending point
* **h**: length of step

**Example code**

```text
double a = 0;
double b = 0.1;
double h = 0.001;

double yEU[101] = { 0, };
	
odeEU(myFunc2, yEU, a, b, h);
```



### odeEM\(\)

Return the Euler's Modified method results for the target equation

```text
void odeEM(double myfunc(const double t, const double y), double y[], double t0, double tf, double h);
```

**Parameters**

* **myfunc**: Function **myfunc** is defined.
* **y\[\]**: solution of ODE
* **t0**: starting point
* **tf**: ending point
* **h**: length of step

**Example code**

```text
double a = 0;
double b = 0.1;
double h = 0.001;

double yEM[101] = { 0, };

odeEM(myFunc2, yEM, a, b, h);
```



### odeRK2\(\)

Return the Runge-Kutta 2nd order method results for the target equation

```text
void odeRK2(double myfunc(const double t, const double y), double y[], double t0, double tf, double h);
```

**Parameters**

* **myfunc**: Function **myfunc** is defined.
* **y\[\]**: solution of ODE
* **t0**: starting point
* **tf**: ending point
* **h**: length of step

**Example code**

```text
double a = 0;
double b = 0.1;
double h = 0.001;

double yRK2[101] = { 0, };

odeRK2(myFunc2, yRK2, a, b, h);
```



### odeRK3\(\)

Return the Runge-Kutta 3rd order method results for the target equation

```text
void odeRK3(double myfunc(const double t, const double y), double y[], double t0, double tf, double h);
```

**Parameters**

* **myfunc**: Function **myfunc** is defined.
* **y\[\]**: solution of ODE
* **t0**: starting point
* **tf**: ending point
* **h**: length of step

**Example code**

```text
double a = 0;
double b = 0.1;
double h = 0.001;

double yRK3[101] = { 0, };

odeRK3(myFunc2, yRK3, a, b, h);
```



### ode\(\)

Return the called method results for the target equation

```text
void ode(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, int method);
```

**Parameters**

* **myfunc**: Function **myfunc** is defined.
* **y\[\]**: solution of ODE
* **t0**: starting point
* **tf**: ending point
* **h**: length of step
* **method**: selection of method

**Example code**

```text
double a = 0;
double b = 0.1;
double h = 0.001;

double yRK3[101] = { 0, };

ode(myFunc2, yRK3, a, b, h, RK3);
```



### sys2RK2\(\)

Return the Runge-Kutta 2nd order method results for the 2nd order target equation

```text
void sys2RK2(void func(const double t, const double Y[], double dYdt[], int in_type), double y[], double z[], double t0, double tf, double h, double y_init, double z_init, int input_type);
```

**Parameters**

* **func**: Function **func** is defined.
* **y\[\]**: y and gradient of y
* **z[]**: gradient of **y[]**
* **t0**: starting point
* **tf**: ending point
* **h**: length of step
* **y_init**: initial value of **y\[\]**
* **z_init**: initial value of **z\[\]**
* **input_type**: selection of input type

**Example code**

```text
double a = 0;
double b = 1;
double h = 0.01;

double yHR[101] = { 0, };
double zHR[101] = { 0, };
double yFV[101] = { 0, };
double zFV[101] = { 0, };
double ySR[101] = { 0, };
double zSR[101] = { 0, };

sys2RK2(mckfunc, yHR, zHR, a, b, h, 0, 0, HR);
sys2RK2(mckfunc, yFV, zFV, a, b, h, 0, 0.2, FV);
sys2RK2(mckfunc, ySR, zSR, a, b, h, 0, 0, SR);
```



### sys2RK3\(\)

Return the Runge-Kutta 3rd order method results for the 2nd order target equation

```text
void sys2RK3(void func(const double t, const double Y[], double dYdt[], int in_type), double y[], double z[], double t0, double tf, double h, double y_init, double z_init, int input_type);
```

**Parameters**

* **func**: Function **func** is defined.
* **y\[\]**: y and gradient of y
* **z[]**: gradient of **y[]**
* **t0**: starting point
* **tf**: ending point
* **h**: length of step
* **y_init**: initial value of **y\[\]**
* **z_init**: initial value of **z\[\]**
* **input_type**: selection of input type

**Example code**

```text
double a2 = 0;
double b2 = 2;
double h2 = 0.02;

double y2[101] = { 0, };
double z2[101] = { 0, };

sys2RK3(mckfunc, y2, z2, a2, b2, h2, 0.05, 0, FV);
```



## Linear Solver

### gaussElim\(\)

Apply Gauss elimination

```cpp
void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d);
```

**Parameters**

* **A**: matrix **A** in structure Matrix form, should be \(nxn\) square
* **b**: vector **b** in structure Matrix form, should be \(nx1\)
* **U**: matrix **U** in structure Matrix form, should be \(nxn\) square
* **d**: vector **d** in structure Matrix form, should be \(nx1\)

**Example code**

```text
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix vecb = txt2Mat(path, "prob1_vecb");
Matrix matU = zeros(matA.rows, matA.cols);
Matrix vecd = zeros(vecb.rows, vecb.cols);

gaussElim(matA, vecb, matU, vecd);
```



### LUdecomp\(\)

Apply LU decomposition

```cpp
void LUdecomp(Matrix A, Matrix L, Matrix U);
```

**Parameters**

* **A**: matrix **A** in structure Matrix form, should be \(nxn\) square
* **L**: matrix **L** in structure Matrix form, should be \(nxn\) square
* **U**: matrix **U** in structure Matrix form, should be \(nxn\) square

**Example code**

```text
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix matU = zeros(matA.rows, matA.cols);
Matrix matL = eye(matA.rows, matA.cols);

LUdecomp(matA, matL, matU);
```



### solveLU\(\)

Solve equation after LU decomposition

```cpp
void solveLU(Matrix L, Matrix U, Matrix b, Matrix x);
```

**Parameters**

* **L**: matrix **L** in structure Matrix form, should be \(nxn\) square
* **U**: matrix **U** in structure Matrix form. should be \(nxn\) square
* **b**: vector **b** in structure Matrix form, should be \(nx1\)
* **x**: vector **x** in structure Matrix form, should be \(nx1\)

**Example code**

```text
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix vecb = txt2Mat(path, "prob1_vecb");
Matrix matU = zeros(matA.rows, matA.cols);
Matrix matL = eye(matA.rows, matA.cols);
Matrix vecx = zeros(matA.cols, vecb.cols);

LUdecomp(matA, matL, matU);
solveLU(matL, matU, vecb, vecx);
```



### fwdsub\(\)

Apply forward-substitution after LU decomposition

```cpp
void fwdsub(Matrix L, Matrix b, Matrix y);
```

**Parameters**

* **L**: matrix **L** in structure Matrix form, should be \(nxn\) square
* **b**: vector **b** in structure Matrix form, should be \(nx1\)
* **y**: vector **y** in structure Matrix form, should be \(nx1\)

**Example code**

```text
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix vecb = txt2Mat(path, "prob1_vecb");
Matrix matU = zeros(matA.rows, matA.cols);
Matrix matL = eye(matA.rows, matA.cols);

LUdecomp(matA, matL, matU);

Matrix y = zeros(L.cols, b.cols);

fwdsub(L, b, y);
```



### backsub\(\)

Apply back-substitution after LU decomposition

```cpp
void backsub(Matrix U, Matrix y, Matrix x);
```

**Parameters**

* **U**: matrix **U** in structure Matrix form, should be \(nxn\) square
* **y**: vector **y** in structure Matrix form, should be \(nx1\)
* **x**: vector **x** in structure Matrix form, should be \(nx1\)

**Example code**

```text
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix vecb = txt2Mat(path, "prob1_vecb");
Matrix matU = zeros(matA.rows, matA.cols);
Matrix matL = eye(matA.rows, matA.cols);
Matrix vecx = zeros(matA.cols, vecb.cols);

LUdecomp(matA, matL, matU);

Matrix y = zeros(L.cols, b.cols);

fwdsub(L, b, y);
backsub(U, y, x);
```



### inv\(\)

Find inverse matrix

```text
void inv(Matrix A, Matrix Ainv);
```

**Parameters**

* **A**: matrix **A** in structure Matrix form, should be \(nxn\) square
* **Ainv**: matrix **Ainv** in structure Matrix form, should be \(nxn\) square

**Example code**

```text
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix matAinv = zeros(matA.rows, matA.cols);

inv(matA, matAinv);
```



## Eigen Value Problem

### QRdecomp()

Apply QR decomposition

```text
void QRdecomp(Matrix _A, Matrix _Q, Matrix _R);
```

**Parameters**

* **A**: matrix **A** in structure Matrix form, should be \(nxn\) square
* **Q**: matrix **Q** in structure Matrix form, should be \(nxn\) square
* **R**: matrix **R** in structure Matrix form, should be \(nxn\) square

**Example code**

```text
Matrix matA = txt2Mat(path, "matA");
Matrix matQ = createMat(matA.rows, matA.cols);
Matrix matR = createMat(matA.rows, matA.cols);

QRdecomp(matA, matQ, matR);
```



### eig()

Return eigenvalues

```text
Matrix eig(Matrix _A);
```

**Parameters**

* **A**: matrix **A** in structure Matrix form, should be \(nxn\) square

**Example code**

```text
Matrix matA = txt2Mat(path, "matA");
Matrix valEig;

valEig = eig(matA);
```



### eigvec()

Return eigenvectors

```text
Matrix eigvec(Matrix _A);
```

**Parameters**

* **A**: matrix **A** in structure Matrix form, should be \(nxn\) square

**Example code**

```text
Matrix matA = txt2Mat(path, "matA");
Matrix vecEig;

vecEig = eigvec(matA);
```



## Curve Fitting

### linearRegression()

Apply linear regression

```text
void linearRegression(Matrix x, Matrix y, Matrix z);
```

**Parameters**

* **x**: vector **x** in structure Matrix form, should be \(mx1\)
* **y**: vector **y** in structure Matrix form, should be \(mx1\)
* **z**: vector **z** in structure Matrix form, should be \(2x1\)

**Example code**

```text
Matrix vecx1 = txt2Mat(path, "vecx1");
Matrix vecy1 = txt2Mat(path, "vecy1");
Matrix z1 = zeros(2, 1);

linearRegression(vecx1, vecy1, z1);
```



### polyfit()

Apply ploynomial curve fit

```text
void polyfit(Matrix x, Matrix y, Matrix z, int n);
```

**Parameters**

* **x**: vector **x** in structure Matrix form, should be \(mx1\)
* **y**: vector **y** in structure Matrix form, should be \(mx1\)
* **z**: vector **z** in structure Matrix form, should be \((n+1)x1\)
* **n**: order of polynomial equation

**Example code**

```text
Matrix vecx2 = txt2Mat(path, "vecx2");
Matrix vecy2 = txt2Mat(path, "vecy2");
int n1 = 3;
Matrix z2 = zeros(n1 + 1, 1);

polyfit(vecx2, vecy2, z2, n1);
```



### logfit()

Apply log curve fit

```text
void logfit(Matrix x, Matrix y, double z[]);
```

**Parameters**

* **x**: vector **x** in structure Matrix form, should be \(mx1\)
* **y**: vector **y** in structure Matrix form, should be \(mx1\)
* **z[]**: parameters of log equation

**Example code**

```text
double X[] = { 20, 22,    24,    26,    28,    30,    32,    34,    36,    38 };
double S[] = { 16.6314,   19.7724,   23.3386,   27.3063,   28.3004,   30.5335,   34.7548,   39.9840,   43.5914,   47.3012 };
Matrix X1 = createMat(sizeof(X) / sizeof(double), 1);
Matrix S1 = createMat(sizeof(S) / sizeof(double), 1);
double zhat[2] = { 0,0 };

for (int i = 0; i < sizeof(X) / sizeof(double); i++)
{
	X1.at[i][0] = X[i];
	S1.at[i][0] = S[i];
}

logfit(X1, S1, zhat);
```



## Non-Linear Solver

### newtonRaphson\(\)

Apply Newton-Raphson method

```text
double newtonRaphson(double func(const double x), double dfunc(const double x), double x0, double tol);
```

**Parameters**

* **func**: Function **func** is defined.
* **dfunc**: Function **dfunc** is defined.
* **x0**: initial value
* **tol**: tolerance error

**Example code**

```text
double x0 = 2.1;
double tol = pow(10, -9);
double x_sol = 0;

x_sol = newtonRaphson(func, dfunc, x0, tol);
```



### newtonRaphson2\(\)

Apply Newton-Raphson method for system

```text
void newtonRaphson2(double x0, double y0, double tol, double ans[]);
```

**Parameters**

* **x0**: initial value of x
* **y0**: initial value of y
* **tol**: tolerance error
* **ans[]**: solution of non-linear equation

**Example code**

```text
double x = -1;
double y = 5;
double tol = pow(10, -9);
double sys_ans[2] = { 0,0 };

newtonRaphson2(x, y, tol, sys_ans);
```

