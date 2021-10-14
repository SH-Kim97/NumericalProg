/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : ���ȯ
Created          : 26-03-2018
Modified         : 14-10-2021
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
double myFunc(const double x);

// Define a function that defines the target equation (2 variable)
double myFunc2(const double t, const double y);

// Return the dy/dx results for the input data (truncation error: O(h^2))
void gradient1D(double x[], double y[], double dydx[], int m);

// Return the dy/dx results for the target equation (truncation error: O(h^2))
void gradientFunc(double func(const double x), double x[], double dydx[], int m);

// Function callback
void func_call(double func(const double x), double xin);

// Prints 1D array
void printVec(double* _vec, int _row);

// Integration using rectangular method for discrete data inputs
double IntegrateRect(double x[], double y[], int m);

// Integration using trapezoidal method for discrete data inputs
double trapz(double x[], double y[], int m);

// Integration using Simpson1/3 method for function inputs
double integral(double func(const double x), double a, double b, int n);

// Integration using Simpson3/8 method for function inputs
double simpson38(double func(const double x), double a, double b, int n);

// Return the Euler method results for the target equation
void odeEU(double myfunc(const double t, const double y), double y[], double t0, double tf, double h);

// Return the Euler's Modified method results for the target equation
void odeEM(double myfunc(const double t, const double y), double y[], double t0, double tf, double h);

// Return the Runge-Kutta 2nd order method results for the target equation
void odeRK2(double myfunc(const double t, const double y), double y[], double t0, double tf, double h);

// Return the Runge-Kutta 3rd order method results for the target equation
void odeRK3(double myfunc(const double t, const double y), double y[], double t0, double tf, double h);

// Return the called method results for the target equation
void ode(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, int method);

// Governing equation and return vector differentiation results for the target vector
void mckfunc(const double t, const double Y[], double dYdt[], int in_type);

// Return the Runge-Kutta 2nd order method results for the 2nd order target equation
void sys2RK2(void func(const double t, const double Y[], double dYdt[], int in_type), double y[], double z[], double t0, double tf, double h, double y_init, double z_init, int input_type);

// Matrix addition
extern	Matrix	addMat(Matrix _A, Matrix _B);

// Apply back-substitution
extern	Matrix	backSub(Matrix _U, Matrix _b);

#endif