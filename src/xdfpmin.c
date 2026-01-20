/***********************************************************
 *  This is driver for the DFP optimization routine (Davidon–Fletcher–Powell method)
 *  (C) Copr. 1986-92 Numerical Recipes Software
 *  Originally provide by Tianlin Wang **EvoRadical (Wong, Sainudiin, & Nielsen 2006)**
 *  Add annotations and minor revised by Hayate Takeuchi (Final Edit: December 14 2025)
 ***********************************************************/

#include <stddef.h>

#define ITMAX 300             // Maximum number of iterations
#define EPS 3.0e-6            // Convergence criterion for numerical precision
#define ALF 1.0e-4            // Sufficient decrease parameter for line search
#define TOLX (4*EPS)          // Step size tolerance
#define STPMX 100.0           // Maximum allowed step length

#define NR_END 1              // For 1-based indexing used in dvector/dmatrix
#define FREE_ARG char*

// Macro to free all dynamically allocated memory used in dfpmin()
#define FREEALL \
    free_dvector(xi, 1, n); \
    free_dvector(pnew, 1, n); \
    free_dmatrix(hessin, 1, n, 1, n); \
    free_dvector(hdg, 1, n); \
    free_dvector(g, 1, n); \
    free_dvector(dg, 1, n); \
    return;

#define FMAX(a, b) ((a) > (b) ? (a) : (b))   // Max of two values
#define SQR(a) ((a) * (a))                   // Square of a value

// Whether to use central difference (1) or forward difference (0) in gradient
int AlwaysCenter = 1;

// Step size for finite difference derivative approximation
double Small_Diff = .5e-6;

/****** External Functions ******/

// User-defined function to minimize (must be defined elsewhere)
extern double func(double x[], int n);

/* Computes gradient numerically using finite differences */
int gradient(int n, double x[], double f0, double g[],
             double (*fun)(double x[], int n), double space[],
             int Central);

/* Main DFP optimization function */
void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
            double (*func)(double[], int n));

/* Line search algorithm used within DFP */
void lnsrch(int n, double xold[], double fold, double g[], double p[],
            double x[], double *f, double stpmax, int *check,
            double (*func)(double[], int));

/****** Memory Allocation Utilities ******/

// Allocate/free double vectors and matrices with 1-based indexing
double *dvector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dvector(double *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

/****** Compute the gradient of the objective function at point x ******/
int gradient(int n, double x[], double f0, double g[],
             double (*fun)(double x[], int n), double space[],
             int Central)
{
    int i, j;
    double *x0 = space;           // Temporary copy of x (for backward)
    double *x1 = space + n;       // Temporary copy of x (for forward)
    double eh0 = Small_Diff;      // Base step size
    double eh;

    if (Central) {
        // Central difference: more accurate but requires 2 function evaluations per variable
        for (i = 1; i <= n; i++) {
            // Copy x into x0 and x1
            for (j = 1; j <= n; j++)
                x0[j] = x1[j] = x[j];

            // Adaptive step size depending on magnitude of x[i]
            eh = pow(eh0 * (fabs(x[i]) + 1), 0.67);
            x0[i] -= eh;
            x1[i] += eh;

            // Central difference approximation: [f(x+h) - f(x-h)] / 2h
            g[i] = ((*fun)(x1, n) - (*fun)(x0, n)) / (2.0 * eh);
        }
    } else {
        // Forward difference: cheaper but less accurate
        for (i = 1; i <= n; i++) {
            for (j = 1; j <= n; j++)
                x1[j] = x[j];

            // Step size proportional to the magnitude of x[i]
            eh = eh0 * (fabs(x[i]) + 1);
            x1[i] += eh;

            // Forward difference: [f(x+h) - f(x)] / h
            g[i] = ((*fun)(x1, n) - f0) / eh;
        }
    }

    // Debug print: show computed gradients
    printf("\nGradients: ");
    for (i = 1; i <= n; i++)
        printf("\t%f", g[i]);

    return 0;
}

/****** Minimize a multivariable function using the DFP method ******/
void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
            double (*func)(double[], int n))
{
    int check, i, its, j;
    double den, fac, fad, fae, fp, stpmax, sum = 0.0, sumdg, sumxi, temp, test;
    double *dg, *g, *hdg, **hessin, *pnew, *xi;

    // Gradient and update vectors
    double *tv;
    
    dg = dvector(1, n);       // g_new - g_old
    g  = dvector(1, n);       // current gradient
    hdg = dvector(1, n);      // Hessian * dg
    hessin = dmatrix(1, n, 1, n);  // Approximate inverse Hessian
    pnew = dvector(1, n);     // new position
    xi = dvector(1, n);       // current search direction
    
    tv = dvector(1, 2 * n);   // workspace for gradient computation

    // Initial function and gradient evaluation
    fp = (*func)(p, n);
    gradient(n, p, fp, g, func, tv, AlwaysCenter);

    // Initialize Hessian to identity matrix and direction to steepest descent
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++)
            hessin[i][j] = 0.0;
        hessin[i][i] = 1.0;
        xi[i] = -g[i];        // initial direction: negative gradient
        sum += p[i] * p[i];   // for step size scaling
    }

    stpmax = STPMX * FMAX(sqrt(sum), (double) n);  // maximum allowed step

    for (its = 1; its <= ITMAX; its++) {
        printf("\nDfpmin Step %03d: Value of the function = %10.5f\n", its, fp);
        for (i = 1; i <= n; i++)
            printf("\t%f", p[i]);

        *iter = its;

        // Perform line search to update p → pnew
        lnsrch(n, p, fp, g, xi, pnew, fret, stpmax, &check, func);
        fp = *fret;

        // Compute step taken
        for (i = 1; i <= n; i++) {
            xi[i] = pnew[i] - p[i];
            p[i] = pnew[i];
        }

        // Check for small change in position → convergence
        test = 0.0;
        for (i = 1; i <= n; i++) {
            temp = fabs(xi[i]) / FMAX(fabs(p[i]), 1.0);
            if (temp > test)
                test = temp;
        }

        if (test < TOLX) {
            printf("\nGradients in the last step (TOLX exit): ");
            for (i = 1; i <= n; i++)
                printf("\t%f", g[i]);
            FREEALL
        }

        // Save old gradient and compute new gradient
        for (i = 1; i <= n; i++)
            dg[i] = g[i];
        gradient(n, p, fp, g, func, tv, AlwaysCenter);

        // Check for small gradient → convergence
        test = 0.0;
        den = FMAX(*fret, 1.0);
        for (i = 1; i <= n; i++) {
            temp = fabs(g[i]) * FMAX(fabs(p[i]), 1.0) / den;
            if (temp > test)
                test = temp;
        }

        if (test < gtol) {
            printf("\ngraidents in the last step: ");
			printf("\nexiting from gtol");
            for (i = 1; i <= n; i++)
                printf("\t%f", g[i]);
            FREEALL
        }

        // Compute difference of gradients
        for (i = 1; i <= n; i++)
            dg[i] = g[i] - dg[i];

        // Compute Hessian × dg
        for (i = 1; i <= n; i++) {
            hdg[i] = 0.0;
            for (j = 1; j <= n; j++)
                hdg[i] += hessin[i][j] * dg[j];
        }

        // Update inverse Hessian matrix using DFP formula
        fac = fae = sumdg = sumxi = 0.0;
        for (i = 1; i <= n; i++) {
            fac   += dg[i] * xi[i];
            fae   += dg[i] * hdg[i];
            sumdg += SQR(dg[i]);
            sumxi += SQR(xi[i]);
        }

        if (fac > sqrt(EPS * sumdg * sumxi)) {
            fac = 1.0 / fac;
            fad = 1.0 / fae;
            for (i = 1; i <= n; i++)
                dg[i] = fac * xi[i] - fad * hdg[i];

            for (i = 1; i <= n; i++) {
                for (j = i; j <= n; j++) {
                    hessin[i][j] += fac * xi[i] * xi[j]
                                  - fad * hdg[i] * hdg[j]
                                  + fae * dg[i] * dg[j];
                    hessin[j][i] = hessin[i][j];  // enforce symmetry
                }
            }
        } else {
            printf("\nSkipped Hessian update due to small curvature condition.\n");
        }

        // Compute next direction: -H * g
        for (i = 1; i <= n; i++) {
            xi[i] = 0.0;
            for (j = 1; j <= n; j++)
                xi[i] -= hessin[i][j] * g[j];
        }
    }

    /* Reached ITMAX iterations without satisfying gtol.
       Treat the current point as the final solution instead of aborting. */
    fprintf(stderr,
            "\n[WARN] dfpmin reached ITMAX=%d without meeting gtol.\n"
            "       Proceeding with the last parameter set as final.\n",
            ITMAX);

    *iter = ITMAX;
    *fret = fp;

    FREEALL
}

/****** Line search algorithm used within DFP ******/
void lnsrch(int n, double xold[], double fold, double g[], double p[],
            double x[], double *f, double stpmax, int *check,
            double (*func)(double[], int)) {

    int i;
    double a, b, disc, f2;
    double alam, alam2 = 0, alamin, rhs1, rhs2;
    double slope, sum = 0, temp, test, tmplam;

    *check = 0;

    // Calculate max step size based on p[] vector norm
    for (i = 1; i <= n; i++)
        sum += p[i] * p[i];
    sum = sqrt(sum);
    if (sum > stpmax)
        for (i = 1; i <= n; i++)
            p[i] *= stpmax / sum;

    // Compute slope = ∇f ⋅ p (directional derivative)
    slope = 0.0;
    for (i = 1; i <= n; i++)
        slope += g[i] * p[i];

    // If slope is non-negative, direction is not a descent direction
    if (slope >= 0.0)
        error("Roundoff problem in lnsrch: not a descent direction.");

    // Determine minimum step size (based on machine tolerance and relative size)
    test = 0.0;
    for (i = 1; i <= n; i++) {
        temp = fabs(p[i]) / FMAX(fabs(xold[i]), 1.0);
        if (temp > test) test = temp;
    }
    alamin = TOLX / test;

    // Start with full step length
    alam = 1.0;

    // Main backtracking loop
    for (;;) {
        for (i = 1; i <= n; i++)
            x[i] = xold[i] + alam * p[i];  // new trial point

        *f = (*func)(x, n);  // evaluate function at new point

        // If step is too small, restore original point and report failure
        if (alam < alamin) {
            for (i = 1; i <= n; i++)
                x[i] = xold[i];
            *check = 1;
            return;
        }

        // Check Armijo (sufficient decrease) condition
        if (*f <= fold + ALF * alam * slope)
            return;

        // Use quadratic or cubic backtracking to update alam
        if (alam == 1.0) {
            tmplam = -slope / (2.0 * (*f - fold - slope));
        } else {
            rhs1 = *f - fold - alam * slope;
            rhs2 = f2 - fold - alam2 * slope;
            a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
            b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);

            if (a == 0.0) {
                tmplam = -slope / (2.0 * b);
            } else {
                disc = b * b - 3.0 * a * slope;
                if (disc < 0.0) {
                    tmplam = 0.5 * alam;
                } else if (b <= 0.0) {
                    tmplam = (-b + sqrt(disc)) / (3.0 * a);
                } else {
                    tmplam = -slope / (b + sqrt(disc));
                }
            }

            if (tmplam > 0.5 * alam)
                tmplam = 0.5 * alam;
        }

        alam2 = alam;
        f2 = *f;
        alam = FMAX(tmplam, 0.1 * alam);  // Avoid step too small
    }
}

/****** Memory allocation utilities for vectors and matrices ******/

double *dvector(long nl, long nh)
/* Allocate a double vector with subscript range v[nl..nh] */
{
    double *v;
    v = (double *)malloc((unsigned int)((nh - nl + 1 + NR_END) * sizeof(double)));
    if (!v) error("allocation failure in dvector()");
    return v - nl + NR_END;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* Allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double **m;

    // Allocate row pointers
    m = (double **)malloc((unsigned int)((nrow + NR_END) * sizeof(double *)));
    if (!m) error("allocation failure 1 in dmatrix()");
    m += NR_END;
    m -= nrl;

    // Allocate actual data block and point rows to slices
    m[nrl] = (double *)malloc((unsigned int)((nrow * ncol + NR_END) * sizeof(double)));
    if (!m[nrl]) error("allocation failure 2 in dmatrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    return m;
}

int *ivector(long nl, long nh)
/* Allocate an int vector with subscript range v[nl..nh] */
{
    int *v;
    v = (int *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int)));
    if (!v) error("allocation failure in ivector()");
    return v - nl + NR_END;
}

void free_dvector(double *v, long nl, long nh)
/* Free a double vector allocated by dvector() */
{
    free((FREE_ARG)(v + nl - NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* Free a double matrix allocated by dmatrix() */
{
    free((FREE_ARG)(m[nrl] + ncl - NR_END));
    free((FREE_ARG)(m + nrl - NR_END));
}

#undef ALF
#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX
#undef FREEALL