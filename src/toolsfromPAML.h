/*************************************
 *  toolsfromPAML.h
 *  functions from PAML
 *  Originally provided by **EvoRadical (Wong, Sainudiin, & Nielsen 2006)**
 *  Add annotations by Hayate Takeuchi (Final Edit: March 31 2025)
 *************************************/

#ifndef TOOLSFROMPAML_H
#define TOOLSFROMPAML_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* Macro definitions */
#define FOR(i,n) for(i=0; i<n; i++)       /* A loop cycles i from 0 to n-1 */
#define FPN(file) fputc('\n', file)       /* Output newline to file */
#define min2(a,b) ((a)<(b)?(a):(b))       /* min() */
#define max2(a,b) ((a)>(b)?(a):(b))       /* Max() */

/* For error handling */
void error(char *message);

/* Basic operations on vectors and matrices */
int zero(double x[], int n);                        /* Initialize the array x with zeros */
double sum(double x[], int n);                      /* Calculate the sum of the array x */
int fillxc(double x[], double c, int n);            /* Assign the value c to the array x */
int xtoy(double x[], double y[], int n);            /* Copy x to y */
int axtoy(double a, double x[], double y[], int n); /* Multiply x by the scalar a and assign the result to y */
int identity(double x[], int n);                    /* Generates a unit matrix (1-dimensional representation) */

/* Matrix arithmetic functions */
int PMatUVRoot(double P[], double t, int n, double U[], double V[], double Root[]); /* P(t) = U * exp(D*t) * V */
int matinv(double x[], int n, int m, double space[]); /* Find the inverse of the matrix x */

/* Eigenvalue/Eigenvector calculation */
int eigen(int job, double A[], int n, double rr[], double ri[], double vr[], double vi[], double w[]);

/* Complex types and operation definitions */
typedef struct { double re, im; } complex; /* Complex type definition */
#define csize(a) (fabs(a.re) + fabs(a.im)) /* Size of a complex */

complex compl(double re, double im);       /* Generate a complex */
//complex conj(complex a);                   /* Returns the complex conjugate */
complex cplus(complex a, complex b);       /* Complex addition */
complex cminus(complex a, complex b);      /* Complex subtraction */
complex cby(complex a, complex b);         /* Complex multiplication */
complex cdiv(complex a, complex b);        /* Complex division */
//complex cexp(complex a);                   /* Complex exponential function */
complex cfactor(complex x, double a);      /* Multiply complex number x by scalar a */

int cxtoy(complex x[], complex y[], int n);                             /* Complex vector copy */
int cmatby(complex a[], complex b[], complex c[], int n,int m,int k);   /* Complex matrix multiplication */
int cmatout(FILE *fout, complex x[], int n, int m);                     /* Output complex matrix to file */
int cmatinv(complex x[], int n, int m, double space[]);                 /* Inverse of a complex matrix */

#endif