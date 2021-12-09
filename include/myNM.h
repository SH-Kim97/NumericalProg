/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : ±è½ÂÈ¯
Created          : 26-03-2018
Modified         : 09-12-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/


#define	PI 3.14159265358979323846264338327950288419716939937510582
#define EU 0
#define EM 1
#define RK2 2
#define RK3 3
#define HR 0
#define FV 1
#define SR 2
#define _CRT_SECURE_NO_WARNINGS

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include "myMatrix.h"


// declare factorial function
extern double factorial(double _x);

// declare Taylor sine function[rad]
extern double sinTaylor(double x);

// declare Taylor sine function[deg]
extern double sindTaylor(double d);

// declare compare function[rad]
extern void CompareSinAns(double sin_in, int iter, double approx);

// declare compare function[deg]
extern void CompareSindAns(double sin_in_d, int iter_d, double approx_d);

// Define a function that defines the target equation (1 variable)
extern double myFunc(const double x);

// Define a function that defines the target equation (2 variable)
extern double myFunc2(const double x, const double t);

// Return the dy/dx results for the input data (truncation error: O(h^2))
extern void gradient1D(double x[], double y[], double dydx[], int m);

// Function callback
extern void func_call(double func(const double x), double xin);

// Return the dy/dx results for the target equation (truncation error: O(h^2))
extern void gradientFunc(double func(const double x), double x[], double dydx[], int m);

// Return the 1st partial differentiation results for the target equation
extern void partialGradient(double func(double x, double t), double x, double t, double dydx[], double dydt[]);

// Return the 2nd partial differentiation results for the target equation
extern void partialGradient2(double func(double x, double t), double x, double t, double dy2dx2[], double dy2dt2[]);

// Prints 1D array
extern void printVec(double* _vec, int _row);

// Prints 1D array as range
extern void printVec_range(double* _vec, int r1, int r2);

// Integration using rectangular method for discrete data inputs
extern double IntegrateRect(double x[], double y[], int m);

// Integration using trapezoidal method for discrete data inputs
extern double trapz(double x[], double y[], int m);

// Integration using Simpson1/3 method for function inputs
extern double integral(double func(const double x), double a, double b, int n);

// Integration using Simpson3/8 method for function inputs
extern double simpson38(double func(const double x), double a, double b, int n);

// Return the Euler method results for the target equation
extern void odeEU(double myfunc(const double t, const double y), double y[], double t0, double tf, double h);

// Return the Euler's Modified method results for the target equation
extern void odeEM(double myfunc(const double t, const double y), double y[], double t0, double tf, double h);

// Return the Runge-Kutta 2nd order method results for the target equation
extern void odeRK2(double myfunc(const double t, const double y), double y[], double t0, double tf, double h);

// Return the Runge-Kutta 3rd order method results for the target equation
extern void odeRK3(double myfunc(const double t, const double y), double y[], double t0, double tf, double h);

// Return the called method results for the target equation
extern void ode(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, int method);

// Governing equation and return vector differentiation results for the target vector
extern void mckfunc(const double t, const double Y[], double dYdt[], int in_type);

// Return the Runge-Kutta 2nd order method results for the 2nd order target equation
extern void sys2RK2(void func(const double t, const double Y[], double dYdt[], int in_type), double y[], double z[], double t0, double tf, double h, double y_init, double z_init, int input_type);

// Return the Runge-Kutta 3rd order method results for the 2nd order target equation
extern void sys2RK3(void func(const double t, const double Y[], double dYdt[], int in_type), double y[], double z[], double t0, double tf, double h, double y_init, double z_init, int input_type);

// Matrix addition
extern	Matrix	addMat(Matrix _A, Matrix _B);

// Matrix multiply
extern	Matrix	multMat(Matrix _A, Matrix _B);

// Apply Gauss elimination
extern void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d);

// Apply back-substitution after Gauss elimination
extern	Matrix	backSub(Matrix _U, Matrix _d);

// Apply LU decomposition
extern void LUdecomp(Matrix A, Matrix L, Matrix U);

// Solve equation after LU decomposition
extern void solveLU(Matrix L, Matrix U, Matrix b, Matrix x);

// Apply forward-substitution after LU decomposition
extern void fwdsub(Matrix L, Matrix b, Matrix y);

// Apply back-substitution after LU decomposition
extern void backsub(Matrix U, Matrix y, Matrix x);

// Find inverse matrix
extern void inv(Matrix A, Matrix Ainv);

// Matrix and constant multiply
extern Matrix	multConst(Matrix _A, double c);

// Return norm of vector
extern double	normVec(Matrix _A);

// Apply QR decomposition
extern void QRdecomp(Matrix _A, Matrix _Q, Matrix _R);

// Return eigenvalues
extern Matrix eig(Matrix _A);

// Return eigenvectors
extern Matrix eigvec(Matrix _A);

// Apply linear regression
extern void linearRegression(Matrix x, Matrix y, Matrix z);

// Apply ploynomial curve fit
extern void polyfit(Matrix x, Matrix y, Matrix z, int n);

// Define a original function
extern double myfunc(const double x);

// Define a function F for Newton-Raphson method
extern double func(const double x);

// Define a function dF for Newton-Raphson method
extern double dfunc(const double x);

// Apply Newton-Raphson method
extern double newtonRaphson(double func(const double x), double dfunc(const double x), double x0, double tol);

#endif