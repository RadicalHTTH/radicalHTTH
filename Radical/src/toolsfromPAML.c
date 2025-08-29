/*************************************
 *  toolsfromPAML.c
 *  functions from PAML
 *  Originally provided by **EvoRadical (Wong, Sainudiin, & Nielsen 2006)**
 *  Add annotations by Hayate Takeuchi (Final Edit: March 31 2025)
 *************************************/

#include "toolsfromPAML.h"

/* Genetic code table: GeneticCode[icode][codon] GenetCode[icode][#codon] */
int GeneticCode[][64] = {
    {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,-1,17,
     10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
      9, 9, 9,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
     19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},  // 0:universal

    {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,17,17,
     10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
      9, 9,12,12,16,16,16,16, 2, 2,11,11,15,15,-1,-1,
     19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},  // 1:vertebrate mt.

    {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,-1,17,
     16,16,16,16,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
      9, 9,12,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
     19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},  // 2:yeast mt.

    {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,17,17,
     10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
      9, 9, 9,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
     19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},  // 3:mold mt.

    {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,17,17,
     10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
      9, 9,12,12,16,16,16,16, 2, 2,11,11,15,15,15,15,
     19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},  // 4:invertebrate mt.

    {13,13,10,10,15,15,15,15,18,18, 5, 5, 4, 4,-1,17,
     10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
      9, 9, 9,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
     19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},  // 5:ciliate nucl.

    {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,17,17,
     10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
      9, 9, 9,12,16,16,16,16, 2, 2, 2,11,15,15,15,15,
     19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},  // 6:echinoderm mt.

    {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4, 4,17,
     10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
      9, 9, 9,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
     19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},  // 7:euplotid mt.

    {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,-1,17,
     10,10,10,15,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
      9, 9, 9,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
     19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},  // 8:alternative yeast nucl.

    {13,13,10,10,15,15,15,15,18,18,-1,-1, 4, 4,17,17,
     10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
      9, 9,12,12,16,16,16,16, 2, 2,11,11,15,15, 7, 7,
     19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7},  // 9:ascidian mt.

    {13,13,10,10,15,15,15,15,18,18,-1, 5, 4, 4,-1,17,
     10,10,10,10,14,14,14,14, 8, 8, 5, 5, 1, 1, 1, 1,
      9, 9, 9,12,16,16,16,16, 2, 2,11,11,15,15, 1, 1,
     19,19,19,19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7}   // 10:blepharisma nucl.
    } ;

/***************************************
 * P(t) = U * exp(Root * t) * V
 ***************************************/
int PMatUVRoot (double P[], double t, int n, double U[], double V[], double Root[])
{
    int i,j,k;
    double expt, uexpt, *pP;

    if (t<1e-100) { identity (P, n); return(0); }
    for (k=0,zero(P,n*n); k<n; k++)
        for (i=0,pP=P,expt=exp(t*Root[k]); i<n; i++)
            for (j=0,uexpt=U[i*n+k]*expt; j<n; j++)
                *pP++ += uexpt*V[k*n+j];

    return (0);
}


/**************************
 * Basic Numerical Tools
 **************************/

/* Function initializes all the first n elements 
 * of the array x given as an argument to 0 */
int zero (double x[], int n)
{ int i; FOR (i,n) x[i]=0; return (0);}

/* Function calculates the sum of the first n elements 
 * of array x and returns the sum (double type) */
double sum (double x[], int n)
{ int i; double t=0;  for(i=0; i<n; i++) t += x[i];    return(t); }

/* Function assigns the value c to all the first n elements of an array x */
int fillxc (double x[], double c, int n)
{ int i; FOR (i,n) x[i]=c; return (0);}

/* Function copies the first n elements of array x to array y */
int xtoy (double x[], double y[], int n)
{ int i; for (i=0; i<n; y[i]=x[i],i++) ;  return(0); }

/* Function multiplies each element of array x 
 * by scalar a and writes the results to array y */
int axtoy(double a, double x[], double y[], int n)
{ int i; for (i=0; i<n; y[i] = a*x[i],i++) ;  return(0);}

/* Function stores an n × n identity matrix 
 * in the one-dimensional array x */
int identity (double x[], int n)
{ int i,j;  FOR (i,n)  { FOR(j,n)  x[i*n+j]=0;  x[i*n+i]=1; }  return (0); }

/* Function that sorts x in descending or ascending order 
 * and stores the original indices in rank */
int sort1 (double x[], int n, int rank[], int descending, int space[])
{
    int i,j, it=0, *mark=space;
    double t=0;
    FOR (i, n) mark[i]=1;
    FOR (i, n) {
        for (j=0; j<n; j++)  if (mark[j]) { t=x[j]; it=j++; break; }
        if (descending) {
            for ( ; j<n; j++)
                if (mark[j] && x[j]>t) { t=x[j]; it=j; }
        }
        else {
            for ( ; j<n; j++)
                if (mark[j] && x[j]<t) { t=x[j]; it=j; }
        }
        mark[it]=0;   rank[i]=it;
    }
    return (0);
}

/* Error handler function prints an error message and terminates the program */
void error (char * message)
{ printf("\nError: %s.\n", message); exit(-1); }

int matinv (double x[], int n, int m, double space[])
{

    register int i,j,k;
    int *irow=(int*) space;
    double ee=1.0e-30, t,t1,xmax;
    double det=1.0;

    FOR (i,n)  {
        xmax = 0.;
        for (j=i; j<n; j++)
            if (xmax < fabs(x[j*m+i])) { xmax = fabs(x[j*m+i]); irow[i]=j; }
        det *= xmax;
        if (xmax < ee)   {
            printf("\nDet becomes zero at %3d!\t\n", i+1);
            return(-1);
        }
        if (irow[i] != i) {
            FOR (j,m) {
                t = x[i*m+j];
                x[i*m+j] = x[irow[i]*m+j];
                x[irow[i]*m+j] = t;
            }
        }
        t = 1./x[i*m+i];
        FOR (j,n) {
            if (j == i) continue;
            t1 = t*x[j*m+i];
            FOR(k,m)  x[j*m+k] -= t1*x[i*m+k];
            x[j*m+i] = -t1;
        }
        FOR(j,m)   x[i*m+j] *= t;
            x[i*m+i] = t;
    }
    for (i=n-1; i>=0; i--) {
        if (irow[i] == i) continue;
        FOR(j,n)  {
            t = x[j*m+i];
            x[j*m+i] = x[j*m + irow[i]];
            x[j*m + irow[i]] = t;
        }
    }
    return (0);
}