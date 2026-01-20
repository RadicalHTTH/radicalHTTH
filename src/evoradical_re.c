/*************************************
 *  evoradical_re.h
 *  Copyright, Wendy Wong, 2004 
 *  Originally provided by **EvoRadical (Wong, Sainudiin, & Nielsen 2006)**
 *  source code for finding eigen values of a real nonsymmetric matrix
 *  was written by Tianlin Wang at University of Illinois. I 
 *  Source code for local optimization was adopted from Numerical Recipes
 *  in C (Press et al 1992).
 *  The organization of the program code was greatly inspired by codeml.c
 *  from the PAML package (Yang, 1997-2002). 
 *  Add branch model and remove some models by Hayate Takeuchi (Final Edit: December 14 2025)
 *************************************/

#include "evoradical_re.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// --- Global variables for the Branch model ---
int is_foreground_branch[MAXBRANCHES];  // foreground branch
int branch_gamma_class[MAXBRANCHES];    // gamma-class per branch: 0/1/2
int num_fg1_branches = 0; // counts of labeled branches (#1) for output & df
int num_fg2_branches = 0; // counts of labeled branches (#2) for output & df
double gamma0, gamma1;                  // background, foreground gamma
int eigen_error_flag = 0;
int timer_initialized = 0;
int use_foreground_branches = 1;
static void make_null_out_name(const char *orig, char *buf, size_t bufsize);
void DebugPrintAsciiTree(FILE *fp);

/**
  * Creating a new node
 **/
SEQUENCE new_node(void)
{
    return (malloc(sizeof(NODE)));
}

int Qmatrix(double Root[], double U[], double V[],
            double kappa, double omega, double gamma,
            int AAPartition[], double Q[], int n, int branch_i)
{
    int i, j, k, codon1, codon2, aa1, aa2;
    int ndiff, pos = 0, from[3], to[3];
    double space[64 * 3], *ri = space;

    /* Use normalized 61-codon frequencies as stationary distribution */
    double *pi_cod = SeqData.normCodon_pi;   /* sum_i pi_cod[i] = 1 */

    /* 1. Initialize Q to 0 */
    for (i = 0; i < n * n; i++)
        Q[i] = 0.0;

    /* 2. Fill off-diagonal entries */
    for (i = 0; i < n; i++) {
        codon1 = FROM61[i];
        from[0] = codon1 / 16;
        from[1] = (codon1 / 4) % 4;
        from[2] = codon1 % 4;

        for (j = 0; j < n; j++) {
            codon2 = FROM61[j];
            to[0] = codon2 / 16;
            to[1] = (codon2 / 4) % 4;
            to[2] = codon2 % 4;

            /* only single nucleotide differences */
            for (k = 0, ndiff = 0; k < 3; k++) {
                if (from[k] != to[k]) {
                    ndiff++;
                    pos = k;
                }
            }
            if (ndiff != 1)
                continue;

            /* base mutation rate */
            Q[i * n + j] = 1.0;

            /* transition → multiply by kappa */
            if ((from[pos] + to[pos] == 1) || (from[pos] + to[pos] == 5))
                Q[i * n + j] = kappa;

            /* scale by target codon stationary frequency */
            Q[i * n + j] *= pi_cod[j];

            /* amino acid change handling */
            aa1 = GeneticCode[0][codon1];
            aa2 = GeneticCode[0][codon2];

            if (aa1 != aa2) {
                int is_radical =
                    (AAPartition[aa1 * (aa1 - 1) / 2 + aa2] == 0 ||
                     AAPartition[aa2 * (aa2 - 1) / 2 + aa1] == 0);

                switch (SeqData.Model) {
                    case 1: /* site model alt */
                        if (is_radical)
                            Q[i * n + j] *= gamma;
                        else
                            Q[i * n + j] *= omega;
                        break;

                    case 2: /* site model null */
                        if (is_radical)
                            Q[i * n + j] *= gamma;
                        else
                            Q[i * n + j] *= omega;
                        break;

                    case 3: { /* branch model alt */
                        if (is_radical)
                            Q[i * n + j] *= gamma;
                        Q[i * n + j] *= omega;
                    }   break;

                    case 4: /* branch model null */
                        if (is_radical)
                            Q[i * n + j] *= gamma;
                        else
                            Q[i * n + j] *= omega;
                        break;
                }
            }

            /* enforce reversibility: Q[j,i] = Q[i,j] */
            Q[j * n + i] = Q[i * n + j];
        }
    }

    /* 3. Put diagonal so that each row sums to 0 */
    for (i = 0; i < n; i++)
        Q[i * n + i] = -sum(Q + i * n, n);

    /* 4. Normalize mean substitution rate to 1:
          rate = - sum_i pi_i * Q_ii
     */
    double rate = 0.0;
    for (i = 0; i < n; i++) {
        rate += pi_cod[i] * (-Q[i * n + i]);
    }

    if (rate <= 0.0 || !isfinite(rate)) {
        fprintf(stderr,
                "\n[Qmatrix] invalid mean rate (rate=%g, model=%d, kappa=%g, omega=%g, gamma=%g, branch=%d)\n",
                rate, SeqData.Model, kappa, omega, gamma, branch_i);
        eigen_error_flag = 1;
        return -1;
    }

    double inv_rate = 1.0 / rate;
    for (i = 0; i < n * n; i++)
        Q[i] *= inv_rate;

    /* 5. Eigen decomposition of normalized Q */
    {
        int estatus = eigen(1, Q, n, Root, ri, U, V, space + n);
        if (estatus != 0) {
            fprintf(stderr,
                    "\n[Qmatrix/eigen] estatus=%d (branch=%d, model=%d, kappa=%g, omega=%g, gamma=%g, rate=%g)\n",
                    estatus, branch_i, SeqData.Model, kappa, omega, gamma, rate);
            eigen_error_flag = 1;
            return -1;
        }
    }

    xtoy(U, V, n * n);
    matinv(V, n, n, space);

    return 0;
}

/**
  * Function that constructs a transition probability matrix based on the site model 
  * and initializes the likelihood vector for each node
**/
int InitializeLkl_site(double omega, double gamma, int AAPartition[])
{
    int pattern_i, node_i, temp, n, no_pattern;

    n = SeqData.ncode;
    no_pattern = SeqData.CNo_pattern;

    // Construct the Q matrix and its eigendecomposition
    if (Qmatrix(Root, U, V, SeqData.kappa, omega, gamma, AAPartition, PMat, n, -1) != 0) {
        return -1;
    }        

    // Initialize the likelihood vector for each node
    for (pattern_i = 0; pattern_i < no_pattern; pattern_i++) {
        for (node_i = 0; node_i < tree.no_node; node_i++) {
            if (pattern_i == 0)
                fillxc(nodes[node_i].likelihood,
                       (double)(node_i >= SeqData.No_seq),
                       no_pattern * n);

            if (node_i < SeqData.No_seq) {
                temp = SeqData.Codonseq[node_i][pattern_i];
                if (temp != -2)
                    nodes[node_i].likelihood[temp] = 1.0;
                else
                    printf("\n nodes[%d].likelihood[%d]=-2", node_i, temp);
            }
        }
    }

    return 0;
}

/**
 * Function that constructs a transition probability matrix based on the branch model 
 * and initializes the likelihood vector for each node
 */

int InitializeLkl_branch(double omega, double gamma0_in, double gamma1_in, int AAPartition[])
{
    int pattern_i, node_i, temp, n, no_pattern;

    n = SeqData.ncode;
    no_pattern = SeqData.CNo_pattern;

    // Assign gamma0 and gamma1 to global variables
    gamma0 = gamma0_in;
    gamma1 = gamma1_in;

    // Initializing the likelihood vector
    for (pattern_i = 0; pattern_i < no_pattern; pattern_i++) {
        for (node_i = 0; node_i < tree.no_node; node_i++) {
            if (pattern_i == 0) {
                fillxc(nodes[node_i].likelihood,
                       (double)(node_i >= SeqData.No_seq),
                       no_pattern * n);
            }

            if (node_i < SeqData.No_seq) {
                temp = SeqData.Codonseq[node_i][pattern_i];
                if (temp != -2)
                    nodes[node_i].likelihood[temp] = 1.0;
                else
                    printf("\n nodes[%d].likelihood[%d]=-2", node_i, temp);
            }
        }
    }

    return 0;
}

/**
  * Calculates the likelihood of the nodes up to the root
 **/

int return_p(int node_i)
{
    int i, j, k;
    int child_i, pattern_i, temp;
    int n = 61;
    int no_pattern = SeqData.CNo_pattern;
    double t, lkl;

    for (i = 0; i < nodes[node_i].no_child; i++) {
        child_i = nodes[node_i].child[i];
        if (nodes[child_i].no_child > 0) {
            if (return_p(child_i) != 0) {
                return -1;
            }
        }
    }

    for (i = 0; i < nodes[node_i].no_child; i++) {
        child_i = nodes[node_i].child[i];
        t = nodes[child_i].time;

        // model3
        if (SeqData.Model == 3) {
            int branch_i = nodes[child_i].ibranch;

            /* branch gamma class (0 = background, 1 = #1, 2 = #2) */
            int class_id = 0;

            if (use_foreground_branches) {
                if (branch_i >= 0 && branch_i < MAXBRANCHES) {
                    class_id = branch_gamma_class[branch_i];
                }
            }

            /* default is background gamma */
            double gamma_branch = SeqData.gamma[0][0];
            if (class_id == 1) {
                gamma_branch = SeqData.gamma[0][1];   /* for #1 */
            } else if (class_id == 2) {
                gamma_branch = SeqData.gamma[0][2];   /* for #2 */
            }

            if (Qmatrix(Root, U, V, SeqData.kappa, SeqData.omega[0],
                        gamma_branch, SeqData.AAPartition[0],
                        PMat, n, branch_i) != 0) {
                return -1; // eigen fail
            }
        }

        // model 4
        else if (SeqData.Model == 4) {
            if (child_i == nodes[node_i].child[0]) {
                if (Qmatrix(Root, U, V, SeqData.kappa, SeqData.omega[0],
                            SeqData.gamma[0][0], SeqData.AAPartition[0],
                            PMat, n, -1) != 0) {
                    return -1;
                }
            }
        }

        PMatUVRoot(PMat, t, n, U, V, Root);

#ifdef _OPENMP
#pragma omp parallel for private(pattern_i,j,k,t,lkl,temp) schedule(static)
#endif
        for (pattern_i = 0; pattern_i < no_pattern; pattern_i++) {
            if (nodes[child_i].no_child < 1) {
                temp = SeqData.Codonseq[child_i][pattern_i];
                if (temp != -2) {
                    for (j = 0; j < n; j++) {
                        nodes[node_i].likelihood[n * pattern_i + j]
                            *= PMat[j * n + temp];
                    }
                }
            } else {
                for (j = 0; j < n; j++) {
                    t = 0.0;
                    for (k = 0; k < n; k++) {
                        lkl = nodes[child_i].likelihood[pattern_i * n + k];
                        t += PMat[j * n + k] * lkl;
                    }
                    nodes[node_i].likelihood[pattern_i * n + j] *= t;
                }
            }
        }
    }

    return 0;
}

/**
  * Calculate the log-likelihood of the gamma partition model (site model: model==1,2)
  * Returns the likelihood for each category (gamma, omega combination)
 **/
double return_gammaPartsP(int node_i)
{
    int i, j, k, l;
    int n = 61;
    int no_pattern = SeqData.CNo_pattern;
    int *patterncounts = SeqData.CPatternCounts;
    double *pi = SeqData.normCodon_pi;
    double prob, lkl = 0.0;
    double eps = 1e-300;

    /* only for site models */
    if (!(SeqData.Model == 1 || SeqData.Model == 2)) {
        error("return_gammaPartsP() is only available for model==1,2 (site model).");
    }

    eigen_error_flag = 0;

    /* ----- 1. compute pattern-wise numerators for each (omega, gamma) category ----- */
    for (i = 0; i < SeqData.omegaCats; i++) {
        for (j = 0; j < SeqData.gammaCats[0]; j++) {

            if (InitializeLkl_site(SeqData.omega[i],
                                   SeqData.gamma[0][j],
                                   SeqData.AAPartition[0]) != 0) {
                eigen_error_flag = 1;
                return 0.0;
            }

            if (return_p(tree.root) != 0) {
                eigen_error_flag = 1;
                return 0.0;
            }

            for (l = 0; l < no_pattern; l++) {
                prob = 0.0;
                for (k = 0; k < n; k++)
                    prob += nodes[tree.root].likelihood[l * n + k] * pi[k];

                /* store numerator * prior for each category */
                SeqData.postdistGamma[
                    l * SeqData.omegaCats * SeqData.gammaCats[0]
                    + i * SeqData.gammaCats[0] + j
                ] = prob * SeqData.probGamma[i * SeqData.gammaCats[0] + j];
            }
        }
    }

    /* ----- 2. sum categories for each pattern and accumulate log-likelihood ----- */

    int zero_patterns = 0;            /* how many patterns needed epsilon */
    int warn_printed  = 0;            /* limit the number of warnings printed */

    for (j = 0; j < no_pattern; j++) {
        prob = 0.0;
        for (i = 0; i < SeqData.omegaCats * SeqData.gammaCats[0]; i++)
            prob += SeqData.postdistGamma[j * SeqData.omegaCats * SeqData.gammaCats[0] + i];

        if (prob < eps) {
            if (warn_printed < 10) {
                printf("⚠️ pattern %d: probability is zero, applying epsilon.\n", j);
                warn_printed++;
            }
            prob = eps;
            zero_patterns++;
        }

        lkl += log(prob) * patterncounts[j];
    }

    /* ----- 3. if too many patterns collapsed to epsilon, treat as invalid point ----- */

    if (zero_patterns > 0) {
        double frac_zero = (double)zero_patterns / (double)no_pattern;

        /* threshold: if >= 1% of patterns are numerically zero, we regard this
           parameter set as numerically broken and ask the optimizer to move away. */
        if (frac_zero >= 0.01) {
            fprintf(stderr,
                    "\n[WARN] %d / %d patterns (%.3f) collapsed to epsilon in return_gammaPartsP().\n"
                    "       Treating this parameter point as invalid (huge penalty).\n",
                    zero_patterns, no_pattern, frac_zero);
            eigen_error_flag = 1;
            return 0.0;
        }
    }

    return lkl;
}

/**
 * Function to calculate the log-likelihood for the branch model (model==3,4)
 */
double return_fixGammaP(int node_i)
{
    int i, k, l;
    int n = 61;
    int no_pattern = SeqData.CNo_pattern;
    int *patterncounts = SeqData.CPatternCounts;
    double *pi = SeqData.normCodon_pi;
    double prob, lkl = 0.0;
    double eps = 1e-300;

    if (!(SeqData.Model == 3 || SeqData.Model == 4)) {
        error("return_fixGammaP() is only available for model==3 or 4 (branch models).");
    }

    eigen_error_flag = 0;

    InitializeLkl_branch(SeqData.omega[0], 0.0, 0.0, SeqData.AAPartition[0]);

    if (return_p(tree.root) != 0) {
        eigen_error_flag = 1;
        return 0.0;
    }

    int zero_patterns = 0;
    int warn_printed  = 0;

    for (l = 0; l < no_pattern; l++) {
        prob = 0.0;
        for (k = 0; k < n; k++) {
            prob += nodes[tree.root].likelihood[l * n + k] * pi[k];
        }

        if (prob < eps) {
            if (warn_printed < 10) {
                printf("⚠️ pattern %d: probability is zero, applying epsilon.\n", l);
                warn_printed++;
            }
            prob = eps;
            zero_patterns++;
        }

        lkl += log(prob) * patterncounts[l];
    }

    if (zero_patterns > 0) {
        double frac_zero = (double)zero_patterns / (double)no_pattern;

        if (frac_zero >= 0.01) {
            fprintf(stderr,
                    "\n[WARN] %d / %d patterns (%.3f) collapsed to epsilon in return_fixGammaP().\n"
                    "       Treating this parameter point as invalid (huge penalty).\n",
                    zero_patterns, no_pattern, frac_zero);
            eigen_error_flag = 1;
            return 0.0;
        }
    }

    return lkl;
}

/**
  * Initialize the FROM61 and FROM64 arrays
  */
int setFROM_61_64(void)
{
    int i, n;
    for (i = 0, n = 0; i < 64; i++) {
        // If GeneticCode[0][i] is -1, it is a stop codon
        if (GeneticCode[0][i] == -1)
            FROM64[i] = -1;
        else {
            // Register in FROM61 and also register conversion in FROM64
            FROM61[n] = i;
            FROM64[i] = n++;
        }
    }
    return (0);
}

/**
  * Frees all dynamically allocated memory
  */
void freeall()
{
    int i;

    // Free the node array
    if (nodes)
        free(nodes);

    // Release all arrays (name, sequence, codon sequence)
    for (i = 0; i < SeqData.No_seq; i++) {
        if (SeqData.name[i])
            free(SeqData.name[i]);
        if (SeqData.seq[i])
            free(SeqData.seq[i]);
        if (SeqData.Codonseq[i])
            free(SeqData.Codonseq[i]);
    }

    // Free arrays for matrices and eigenvectors
    if (PMat)  free(PMat);
    if (U)     free(U);
    if (V)     free(V);
    if (Root)  free(Root);
}

int main(int argc, char *argv[])
{
    FILE *fp;
    int length_label, i, no_param;
    char *optionFName;
    
    //optimization
    double x[NP];
    //int iter;
    double fret;

    start_time();

#ifdef _OPENMP
    {
        int max_threads = omp_get_max_threads();
        int nthreads = (max_threads < 20) ? max_threads : 20;
        omp_set_num_threads(nthreads);
        printf("\n[OpenMP] Using %d threads\n", nthreads);
    }
#endif

    SeqData.ncode = 61;
    setFROM_61_64();

    //no option file was specified, assume it is evonc.dat
    if (argc < 2)
        optionFName = "evonc.dat";
    else
        optionFName = argv[1];

    ReadOptionFile(optionFName);

    fp = gfopen(SeqData.seqFName, "r+");
    Readfile(fp);
    fclose(fp);

    FindPattern();

    PMat = (double *) malloc(NCODON * NCODON * sizeof(double));
    U    = (double *) malloc(NCODON * NCODON * sizeof(double));
    V    = (double *) malloc(NCODON * NCODON * sizeof(double));
    Root = (double *) malloc(NCODON * NCODON * sizeof(double));

    fp = gfopen(SeqData.treeFName, "r+");
    new_tree(SeqData.No_seq * 2 - 1);
    ReadTree(fp, &length_label, 0);
    fclose(fp);

    for (i = 0; i < SeqData.no_parts; i++)
        GetAAPart(i, SeqData.partFNames[i]);

    x[0] = 0;

    if (SeqData.Model == 3) {
        char out_alt[MAXLINE];
        char out_null[MAXLINE];

        strcpy(out_alt, SeqData.outFName);
        make_null_out_name(SeqData.outFName, out_null, sizeof(out_null));

        int saved_fg1 = num_fg1_branches;
        int saved_fg2 = num_fg2_branches;

        /* ---------------------------------------------
         * 1st pass: null (Model 3, but no foreground)
         * --------------------------------------------- */
        use_foreground_branches = 0;
        num_fg1_branches = 0;
        num_fg2_branches = 0;

        strcpy(SeqData.outFName, out_null);

        no_param = InitializeX(x);
        for (i = 0; i < 3; i++)
            dfpmin(x, no_param, GTOL, &iter, &fret, func);

        func(x, no_param);
        double LogL_null = tree.LogL;

        fp = gfopen(SeqData.outFName, "w");
        printX(x, fp);
        fclose(fp);

        printf("\n[MODEL3 driver] Null (Model 3, no foreground) finished. LogL_null = %f\n",
               LogL_null);

        /* ---------------------------------------------
         * 2nd pass: alt (Model 3, with foreground)
         * --------------------------------------------- */
        use_foreground_branches = 1;   /* foreground ON */
        num_fg1_branches = saved_fg1;
        num_fg2_branches = saved_fg2;

        strcpy(SeqData.outFName, out_alt);

        no_param = InitializeX(x);
        for (i = 0; i < 3; i++)
            dfpmin(x, no_param, GTOL, &iter, &fret, func);

        func(x, no_param);

        printf("\n[MODEL3 driver] Alt (Model 3, with foreground) finished.\n");
        printf("LogL_null = %f, LogL_alt = %f\n", LogL_null, tree.LogL);

        fp = gfopen(SeqData.outFName, "w");
        printX(x, fp);
        fclose(fp);

        double elapsed_time = prn_time();
        printf("\nRun Time: %.2f s\n", elapsed_time);
    }
    else {
        /* Model 1, 2, 4 */
        no_param = InitializeX(x);

        for (i = 0; i < 3; i++)
            dfpmin(x, no_param, GTOL, &iter, &fret, func);

        func(x, no_param);

        printf("\nThe log likelihood of the tree is: %f", tree.LogL);

        fp = gfopen(SeqData.outFName, "w");
        printX(x, fp);
        fclose(fp);

        double elapsed_time = prn_time();
        printf("\nRun Time: %.2f s\n", elapsed_time);
    }

    //freeall();

    return 0;
}

/**
  * prints out a matrix, modified from PAML
 **/

int matout1(double x[], int n, int m)
{
    int i, j;
    for (i = 0; i < n; i++) {
    printf("\n");
	FOR(j, m) printf(" %11.6f", x[i * m + j]);
    }
    return (0);
}

/* Newick output */
static void fprintNewickRec(FILE *fp, int node_i, int parent_i)
{
    int i;

    if (nodes[node_i].no_child > 0) {
        fputc('(', fp);
        for (i = 0; i < nodes[node_i].no_child; i++) {
            int child = nodes[node_i].child[i];
            if (i > 0) fputc(',', fp);
            fprintNewickRec(fp, child, node_i);
        }
        fputc(')', fp);
    } else {
        int seq_idx = nodes[node_i].seq_no - 1;
        if (seq_idx >= 0 && seq_idx < SeqData.No_seq)
            fprintf(fp, "%s", SeqData.name[seq_idx]);
        else
            fprintf(fp, "N%d", node_i);
    }

    if (parent_i >= 0) {
        double t = nodes[node_i].time;
        if (t < 0.0) t = 0.0;
        fprintf(fp, ":%g", t);
    }
}

static void fprintNewickTreeWithBL(FILE *fp)
{
    fprintNewickRec(fp, tree.root, -1);
    fprintf(fp, ";\n");
}

/**
  * Output estimated parameters, posterior probabilities, etc. to a file
  */
int printX(double x[], FILE *fp) {
    int i, j, pattern_i, category_i;
    int no_category;
    double prob;

    int num_branch = tree.no_edge;
    int df_model = 0;
    if (SeqData.Model == 1)       // site model alt
        df_model = 8 + num_branch;
    else if (SeqData.Model == 2)  // site model null
        df_model = 5 + num_branch;
    else if (SeqData.Model == 3) { // branch model alt
        int n_gamma_effective = 1;  // background
        if (num_fg1_branches > 0) n_gamma_effective++;
        if (num_fg2_branches > 0) n_gamma_effective++;

        df_model = num_branch      /* branch length */
                 + 1               /* kappa */
                 + 1               /* omega */
                 + n_gamma_effective;
    }
    else if (SeqData.Model == 4)  // branch model null
        df_model = 3 + num_branch;

    // Record the model number and input file
    fprintf(fp, "\n--- MODEL SUMMARY ---\n");
    fprintf(fp, "Model used: %d\n", SeqData.Model);
    fprintf(fp, "Sequence file: %s\n", SeqData.seqFName);
    fprintf(fp, "Tree file: %s\n", SeqData.treeFName);
    fprintf(fp, "Log likelihood: %f\n", tree.LogL);
    fprintf(fp, "np: %d\n", df_model);

    // Record the start time
    char start_str[100];

    struct tm *tm_start = localtime(&tk.begin_time);
    if (tm_start) {
        strftime(start_str, sizeof(start_str), "%Y-%m-%d %H:%M:%S", tm_start);
        fprintf(fp, "Start time: %s\n", start_str);
    } else {
        fprintf(fp, "Start time: (unknown)\n");
    }

    // Output of base frequency
    fprintf(fp, "\n--- BASE FREQUENCIES ---\n");
    fprintf(fp, "T: %f  C: %f  A: %f  G: %f\n",
            SeqData.basepi[0], SeqData.basepi[1],
            SeqData.basepi[2], SeqData.basepi[3]);

    // Branching according to selected model
    if (SeqData.Model == 1 || SeqData.Model == 2) {
        // Site model (model==1 or 2)
        // Calculate the number of combinations for each category
        no_category = SeqData.omegaCats * SeqData.gammaCats[0];

        // Normalize posterior probabilities
        for (pattern_i = 0; pattern_i < SeqData.CNo_pattern; pattern_i++) {
            prob = 0.0;
            for (category_i = 0; category_i < no_category; category_i++) {
                prob += SeqData.postdistGamma[pattern_i * no_category + category_i];
            }
            for (category_i = 0; category_i < no_category; category_i++) {
                SeqData.postdistGamma[pattern_i * no_category + category_i] /= prob;
            }
        }

        // Branch length output
        fprintf(fp, "\n--- BRANCH LENGTHS ---\n");
        for (i = 0; i < tree.no_edge; i++) {
            int parent = tree.edges[i][0];
            int child  = tree.edges[i][1];

            double t = transformXLower(x[i + 1], 1e-4);
            t        = transformXUpper(t, 5.0);
            nodes[child].time = t;

            if (nodes[child].seq_no > 0) {
                /* leaf branch: child is a leaf with a name */
                fprintf(fp, "Branch %2d: %f  (parent node %d -> leaf %s)\n",
                        i, t, parent, SeqData.name[nodes[child].seq_no - 1]);
            } else {
                /* internal branch */
                fprintf(fp, "Branch %2d: %f  (parent node %d -> internal node %d)\n",
                        i, t, parent, child);
            }
        }
        
        // newick tree
        fprintf(fp, "\n--- TREE WITH ESTIMATED BRANCH LENGTHS (NEWICK) ---\n");
        fprintNewickTreeWithBL(fp);

        fprintf(fp, "\n--- ASCII PHYLOGRAM (estimated branch lengths) ---\n");
        DebugPrintAsciiTree(fp);

        // Kappa output
        fprintf(fp, "\n--- PARAMETERS ---\n");
        fprintf(fp, "Kappa: %f\n", SeqData.kappa);

        // Omega output
        for (i = 0; i < SeqData.omegaCats; i++) {
            fprintf(fp, "Omega[%d]: %f\n", i, SeqData.omega[i]);
        }

        // Gamma output
        for (i = 0; i < SeqData.gammaCats[0]; i++) {
            fprintf(fp, "Gamma[0][%d]: %f\n", i, SeqData.gamma[0][i]);
        }

        // Mixture Probability (Category Probability) Output
        fprintf(fp, "\n--- CATEGORY PROBABILITIES ---\n");
        for (i = 0; i < SeqData.omegaCats; i++) {
            for (j = 0; j < SeqData.gammaCats[0]; j++) {
                int idx = i * SeqData.gammaCats[0] + j;
                fprintf(fp, "P(Omega %d, Gamma %d): %f\n", i, j, SeqData.probGamma[idx]);
            }
        }

        // Category posterior probabilities by site
        fprintf(fp, "\n--- POSTERIOR PROBABILITIES PER SITE ---\n");
        for (i = 0; i < SeqData.No_codon; i++) {
            pattern_i = SeqData.CMapPattern[i];
            if (pattern_i >= 0) {
                fprintf(fp, "Site %4d: ", i+1);
                for (j = 0; j < no_category; j++) {
                    double post = SeqData.postdistGamma[pattern_i * no_category + j];
                    fprintf(fp, "C%d: %.3f  ", j+1, post);
                }
                fprintf(fp, "\n");
            }
        }
    }

    else if (SeqData.Model == 3 || SeqData.Model == 4) {
        // Branch model

        // Branch length output
        fprintf(fp, "\n--- BRANCH LENGTHS ---\n");
        for (i = 0; i < tree.no_edge; i++) {
            int parent = tree.edges[i][0];
            int child  = tree.edges[i][1];

            double t = transformXLower(x[i + 1], 1e-4);
            t        = transformXUpper(t, 5.0);
            nodes[child].time = t;

            if (nodes[child].seq_no > 0) {
                /* leaf branch: child is a leaf with a name */
                fprintf(fp, "Branch %2d: %f  (parent node %d -> leaf %s)\n",
                        i, t, parent, SeqData.name[nodes[child].seq_no - 1]);
            } else {
                /* internal branch */
                fprintf(fp, "Branch %2d: %f  (parent node %d -> internal node %d)\n",
                        i, t, parent, child);
            }
        }

        // newick tree
        fprintf(fp, "\n--- TREE WITH ESTIMATED BRANCH LENGTHS (NEWICK) ---\n");
        fprintNewickTreeWithBL(fp);

        fprintf(fp, "\n--- ASCII PHYLOGRAM (estimated branch lengths) ---\n");
        DebugPrintAsciiTree(fp);

        fprintf(fp, "\n--- PARAMETERS ---\n");

        // Kappa output
        fprintf(fp, "Kappa: %f\n", SeqData.kappa);

        // Omega output
        fprintf(fp, "Omega: %f\n", SeqData.omega[0]);

        // Gamma output
        if (SeqData.Model == 3) {
            fprintf(fp, "Gamma (background):    %f\n", SeqData.gamma[0][0]);

            if (num_fg1_branches > 0) {
                fprintf(fp,
                        "Gamma (foreground #1): %f  (branches: %d)\n",
                        SeqData.gamma[0][1], num_fg1_branches);
            }

            if (num_fg2_branches > 0) {
                fprintf(fp,
                        "Gamma (foreground #2): %f  (branches: %d)\n",
                        SeqData.gamma[0][2], num_fg2_branches);
            }

            if (num_fg1_branches == 0 && num_fg2_branches == 0) {
                fprintf(fp,
                        "NOTE: no branches labeled #1 or #2 in the tree "
                        "(only background gamma is used).\n");
            }
        }
        else {
            fprintf(fp, "Gamma (all branches):  %f\n", SeqData.gamma[0][0]);
        }
    }

    // Record the end time
    time_t end_time = time(NULL);
    struct tm *t = localtime(&end_time);
    char end_str[100];
    strftime(end_str, sizeof(end_str), "%Y-%m-%d %H:%M:%S", t);
    fprintf(fp, "\nEnd time: %s\n", end_str);


    return 0;
}

/**
  * Function that records the start time of execution (both user time and wall clock time)
  * Copied from A Book On C
**/
void start_time(void)
{
    /* Initialize begin_time only once per run */
    if (!timer_initialized) {
        tk.begin_time = time(NULL);   /* wall-clock at start */
        timer_initialized = 1;
    }

    tk.begin_clock = tk.save_clock = clock();
    tk.save_time  = time(NULL);       /* latest wall-clock */
}


/**
  * Function to calculate and print execution time (both CPU time and real time)
**/
double prn_time(void)
{
    char s1[MAXSTRING], s2[MAXSTRING];
    int field_width, n1, n2;
    double clocks_per_second = (double) CLOCKS_PER_SEC;
    double user_time, real_time;

    // Calculate the elapsed CPU time (user time) (unit: seconds)
    user_time = (clock() - tk.save_clock) / clocks_per_second;

    // Calculate the elapsed real time (UNIX time) (unit: seconds)
    real_time = difftime(time(NULL), tk.save_time);

    tk.save_clock = clock();
    tk.save_time  = time(NULL);

    n1 = sprintf(s1, "%.1f", user_time);
    n2 = sprintf(s2, "%.1f", real_time);
    field_width = (n1 > n2) ? n1 : n2;

    printf("\n\n%s\n%s%*.1f%s\n%s%*.1f%s\n\n",
        "TIME ELAPSED",
        "User time: ", field_width, user_time, " seconds",
        "Real time: ", field_width, real_time, " seconds");

    return user_time;
}

/**
  * Function that converts bases (T, C, A, G) into their corresponding numeric codes
  *
  * 't' → 0
  * 'c' → 1
  * 'a' → 2
  * 'g' → 3
  *
  * Any other character (such as a gap or N) returns -2.
 **/
int Char2Code(char c)
{
    c = tolower(c);  // Convert to lower case in case the input is upper case

    switch (c) {
        case 't': return 0;  // T → 0
        case 'c': return 1;  // C → 1
        case 'a': return 2;  // A → 2
        case 'g': return 3;  // G → 3
        default:
            return -2;
    }
}


/**
  * Converts a codon to an integer that maps to one of 61 standard codons
  *
  * - If the input is not a valid base ('n', '*', '-', etc.), it returns -2.
  * - If it contains a stop codon, it will terminate with an error.
 **/
int Codon2Code61(char codon[3])
{
    int i, code;

    for (i = 0; i < 3; i++) {
        if (codon[i] != 't' && codon[i] != 'c' && codon[i] != 'a' && codon[i] != 'g')
            return -2;
    }

    // Convert codons to codes 0-63
    code = Codon2Code64(codon);

    // 0-63 code corresponds to 61 valid codons
    code = from64_61[code];

    // If it is -1, it is a stop codon and the program ends with an error
    if (code == -1)
        error("data contains stop codon!");

    return code;
}


/**
  * Function that converts a codon into a number between 0 and 63.
  *
  * Conversion method:
  *   1st character code × 16 + 2nd character code × 4 + 3rd character code
 **/
int Codon2Code64(char codon[3])
{
    int code =
        Char2Code(codon[0]) * 16 + Char2Code(codon[1]) * 4 +
        Char2Code(codon[2]);
    return code;
}

/**
  * Function to calculate codon frequency from base frequency
  *
  * - Input: Frequency of bases (T, C, A, G) stored in SeqData.basepi[]
  * - Output:
  *     ・SeqData.pi[]              :61 codon frequencies before normalization
  *     ・SeqData.normCodon_pi[]    :Normalized frequency of 61 codons (normalized to a sum of 1)
 **/
int calCodonFreqs(void)
{
    int i;
    int codon;
    double nCodonNormFact = 0.0;

    // Stop codon (UAA=34, UAG=35, UGA=50)
    int StopCodons[3] = {34, 35, 50};

    // Calculate the product of the base frequencies of the stop codons and sum them
    for (i = 0; i < 3; i++) {
        nCodonNormFact += SeqData.basepi[(StopCodons[i] / 4) % 4]
                        * SeqData.basepi[StopCodons[i] / 16]
                        * SeqData.basepi[StopCodons[i] % 4];
    }

    // Subtract the frequency of stop codons from the total to obtain a normalization factor
    nCodonNormFact = 1.0 - nCodonNormFact;

    // Calculate the frequency of each codon from the product of base frequencies
    for (i = 0; i < 61; i++) {
        codon = FROM61[i];

        SeqData.pi[i] = SeqData.basepi[codon / 16]
                      * SeqData.basepi[(codon / 4) % 4]
                      * SeqData.basepi[codon % 4];

        SeqData.normCodon_pi[i] = SeqData.pi[i] / nCodonNormFact;
    }

    return 0;
}

/**
  * File opening function
  *
  * - Opens the specified file in the specified mode (e.g. "r" for reading/"w" for writing).
 **/
FILE *gfopen(char *filename, char *mode)
{
    FILE *fp;

    if ((fp = fopen(filename, mode)) == NULL) {
        fprintf(stderr, "Cannot open %s -bye!\n", filename);
        exit(1);
    } else {
        printf("\nFile %s opened successfully\n", filename);
    }

    return fp;
}

/* for model 3 output "result.out" → "result_null.out"
 */
static void make_null_out_name(const char *orig, char *buf, size_t bufsize)
{
    const char *dot = strrchr(orig, '.');
    if (!dot) {
        snprintf(buf, bufsize, "%s_null", orig);
        return;
    }

    size_t base_len = (size_t)(dot - orig);
    if (base_len > bufsize - 6)
        base_len = bufsize - 6;

    memcpy(buf, orig, base_len);
    buf[base_len] = '\0';
    strcat(buf, "_null");
    strncat(buf, dot, bufsize - strlen(buf) - 1);
}

/**
  * Function that reads sequence files and saves the base sequence and codon sequence
 **/
int Readfile(FILE * fp)
{
    int i, j, k, linelength;
    char *line, *p;
    char codon[3];

    linelength = 256;

    // Read the number of sequences and the length of the sequence
    assert(fscanf(fp, "%d %d", &SeqData.No_seq, &SeqData.Seq_length) == 2);
    assert(SeqData.No_seq > 0 && SeqData.Seq_length > 0);

    // Check if SEQLENGTH is large enough
    if (SeqData.Seq_length > SEQLENGTH)
        error("SEQLENGTH is smaller than the actual sequence length, please increase SEQLENGTH!");

    // Calculate the number of codons
    SeqData.No_codon = floor(SeqData.Seq_length / 3);

    // Reserve memory for each sequence's name and base sequence
    for (i = 0; i < SeqData.No_seq; i++) {
        if (SeqData.name[i]) free(SeqData.name[i]);
        if (SeqData.seq[i])  free(SeqData.seq[i]);

        SeqData.name[i] = (char *) malloc((NAMELENGTH + 1) * sizeof(char));
        SeqData.seq[i]  = (char *) malloc(SeqData.Seq_length * sizeof(char));

        for (j = 0; j < NAMELENGTH; j++)
            SeqData.name[i][j] = 0;
    }

    // Temporarily allocate memory for one line to be read
    linelength = max2(linelength, SeqData.Seq_length * 2 + NAMELENGTH + 2);
    if ((line = (char *) malloc(linelength * sizeof(char))) == NULL)
        error("out of memory in creating line");

    for (i = 0; i < SeqData.No_seq; i++) {
        int cc = 0;
        int lp = 0;
        while ((cc == '\n') || (cc == ' ') || (cc == '\0') || (cc == '\r'))
            cc = fgetc(fp);
        while ((cc != '\n') && (cc != ' ') && (cc != '\0') && (cc != '\r') && (cc != EOF) && (lp + 1 <= NAMELENGTH)) {
            line[lp] = cc;
            lp++;
            cc = fgetc(fp);
        }

        if (cc == EOF) {
            printf("\nError: Early EOF encountered reading sequence name %d.\n", i);
            printf("Please check your sequence file.\n");
            exit(-1);
        }

        line[lp] = '\0';

        // Save Name
        strncpy(SeqData.name[i], line, NAMELENGTH);
        SeqData.name[i][NAMELENGTH] = '\0';
        printf("Reading sequence: %s\n", SeqData.name[i]);

        // Get the first line of the sequence
        while (1) {
            p = fgets(line, linelength, fp);
            if (p == NULL)
                error("please check the sequence file");
            else if (*p != '\n')
                break;
        }

        // Read the base characters one by one
        for (j = 0; j < SeqData.Seq_length; j++) {
            while (*p == '\0' || *p == '\n' || *p == '\r') {
                p = fgets(line, linelength, fp);
                if (p == NULL)
                    error("check your sequence file");
            }

            while (*p == ' ') p++;

            // Convert uracil to thymine
            if (*p == 'u') *p = 't';

            // If the base is valid, save it in lower case
            if (*p != '\0' && *p != '\n' && *p != '\r') {
                *p = (char) tolower(*p);
                SeqData.seq[i][j] = *p;
                p++;
            } else {
                j--;
            }
        }

        SeqData.seq[i][j] = '\0';
    }

    // Allocate memory for each sequence's codon sequence (integerized)
    for (i = 0; i < SeqData.No_seq; i++) {
        if (SeqData.Codonseq[i]) free(SeqData.Codonseq[i]);
        if ((SeqData.Codonseq[i] = (int *) malloc((SeqData.Seq_length / 3) * sizeof(int))) == NULL)
            error("out of memory in creating space for codon sequence");
    }

    // Convert every 3 bases into codons
    for (i = 0; i < SeqData.No_codon; i++) {
        for (j = 0; j < SeqData.No_seq; j++) {
            for (k = 0; k < 3; k++)
                codon[k] = SeqData.seq[j][3 * i + k];
            SeqData.Codonseq[j][i] = Codon2Code61(codon);
        }
    }

    return (0);
}

/* --------------------------------------------------------------------
 * ASCII phylogram (vertical, like your C++ example)
 *
 * - Root is printed as "(root)"
 * - For each child:
 *      prefix + connector + [bar of '-' or '='] + name / (*)
 * - connector: '|' for intermediate siblings, '`' for the last one
 * - foreground branches are drawn with '=' and marked "[#1]" or "[#2]"
 * - background branches are drawn with '-'
 * -------------------------------------------------------------------- */

static double max_branch_length_recursive(int node, int is_root)
{
    double m = 0.0;
    int i;

    /* Use the branch length from parent to this node, except for root */
    if (!is_root) {
        double t = nodes[node].time;
        if (t > m) m = t;
    }

    /* Recurse into children */
    for (i = 0; i < nodes[node].no_child; i++) {
        int child = nodes[node].child[i];
        double cmax = max_branch_length_recursive(child, 0);
        if (cmax > m) m = cmax;
    }

    return m;
}

/* Recursive drawing function */
static void ascii_tree_draw(int node, const char *prefix, double maxlen, FILE *fp)
{
    const int MAX_BAR = 40;
    const int MIN_BAR = 2;
    int i;

    int nchild = nodes[node].no_child;
    for (i = 0; i < nchild; i++) {
        int child = nodes[node].child[i];
        int last  = (i == nchild - 1);

        /* branch length from node to child */
        double blen = nodes[child].time;
        if (blen < 0.0) blen = 0.0;

        /* scale bar length relative to maxlen */
        int bar = 0;
        if (maxlen > 0.0)
            bar = (int)((blen / maxlen) * MAX_BAR + 0.5);
        if (bar < MIN_BAR) bar = MIN_BAR;

        /* foreground/background decision (class 0/1/2) */
        int ibranch  = nodes[child].ibranch;
        int class_id = 0;

        if (use_foreground_branches) {
            if (ibranch >= 0 && ibranch < MAXBRANCHES) {
                class_id = branch_gamma_class[ibranch];
            }

            if (class_id == 0 && nodes[child].label > 0) {
                if (nodes[child].label == 1 || nodes[child].label == 2)
                    class_id = nodes[child].label;
            }
        }

        int  is_fg     = (class_id > 0);
        char edge_char = is_fg ? '=' : '-';
        char connector = last ? '`' : '|';

        /* 1) print this line */
        fprintf(fp, "%s%c", prefix, connector);
        for (int b = 0; b < bar; b++)
            fputc(edge_char, fp);
        fputc(' ', fp);

        if (nodes[child].no_child == 0 && nodes[child].seq_no > 0) {
            /* leaf */
            fprintf(fp, "%s", SeqData.name[nodes[child].seq_no - 1]);
        } else {
            /* internal node */
            fprintf(fp, "(*)");
        }

        if (class_id == 1)
            fprintf(fp, " [#1]");
        else if (class_id == 2)
            fprintf(fp, " [#2]");

        fputc('\n', fp);

        /* 2) prepare prefix for the children of this child */
        char next_pref[1024];
        int len = snprintf(next_pref, sizeof(next_pref), "%s%c", prefix, (last ? ' ' : '|'));
        if (len < 0) len = 0;
        if (len >= (int)sizeof(next_pref)) len = (int)sizeof(next_pref) - 1;

        /* add spaces equal to bar length */
        int b;
        for (b = 0; b < bar && len + b < (int)sizeof(next_pref) - 1; b++) {
            next_pref[len + b] = ' ';
        }
        next_pref[len + b] = '\0';

        /* 3) recurse */
        ascii_tree_draw(child, next_pref, maxlen, fp);
    }
}

/* Public function to dump the current tree as ASCII */
void DebugPrintAsciiTree(FILE *fp)
{
    if (tree.no_node <= 0) {
        fprintf(fp, "[ASCII] tree.no_node <= 0, nothing to draw.\n");
        return;
    }

    /* get maximum branch length over the tree */
    double maxlen = max_branch_length_recursive(tree.root, 1);
    if (maxlen <= 0.0)
        maxlen = 1.0;

    fprintf(fp, "\n[DEBUG] ASCII phylogram ('-' = background, '=' = foreground)\n");
    fprintf(fp, "(root)\n");
    ascii_tree_draw(tree.root, "", maxlen, fp);
}

/**
 * Read a user-supplied phylogenetic tree in Newick format.
 *
 * - Based on codeml.c (PAML) style.
 * - This version:
 *     * builds the topology (nodes, edges)
 *     * stores branch lengths into nodes[child].time
 *     * stores branch labels (#1 etc.) into nodes[node].label
 *     * AFTER the whole tree is read, maps node labels to
 *       is_foreground_branch[branch_index] (for models 3 and 4).
 */
int ReadTree(FILE *treeFile, int *length_label, int popline)
{
    int cnode, cparent;
    int inodeb;
    int i, j, level, hasname, haslength, haslabel, ch, lline;
    char check[NS], line[255], delimiters[] = "(),:#$;";

    cparent    = -1;
    inodeb     =  0;
    level      =  0;
    haslength  =  0;
    haslabel   =  0;
    ch         = ' ';
    lline      = 255;

    *length_label = 0;
    tree.no_node  = SeqData.No_seq;
    tree.no_edge  = 0;

    /* Initialise node array */
    for (i = 0; i < 2 * SeqData.No_seq - 1; i++) {
        nodes[i].parent  = -1;
        nodes[i].ibranch = -1;
        nodes[i].no_child = 0;
        nodes[i].label   = 0;     /* node label (for branch marking) */
        nodes[i].time    = 0.0;   /* branch length from parent to this node */
    }

    /* VERY IMPORTANT: clear foreground flags */
    for (i = 0; i < MAXBRANCHES; i++) {
        is_foreground_branch[i] = 0;
        branch_gamma_class[i]   = 0;
    }
    num_fg1_branches = 0;
    num_fg2_branches = 0;

    /* Skip initial whitespace */
    while (isspace(ch))
        ch = fgetc(treeFile);
    ungetc(ch, treeFile);

    /* Check that each species appears exactly once */
    for (i = 0; i < SeqData.No_seq; i++)
        check[i] = 0;

    /* --- Main Newick parser loop --- */
    for (;;) {
        ch = fgetc(treeFile);
        if (ch == EOF)
            return -1;
        else if (!isgraph(ch))
            continue;

        /* Opening parenthesis: start a new internal node */
        else if (ch == '(') {
            level++;
            cnode = tree.no_node++;

            if (cparent >= 0) {
                /* create edge: parent -> cnode */
                nodes[cparent].child[nodes[cparent].no_child++] = cnode;
                nodes[cnode].parent = cparent;

                tree.edges[tree.no_edge][0] = cparent;
                tree.edges[tree.no_edge][1] = cnode;
                nodes[cnode].ibranch        = tree.no_edge++;
            } else {
                /* first internal node is the root */
                tree.root = cnode;
            }
            cparent = cnode;
        }

        /* Closing parenthesis: finished a subtree; cparent is the internal node */
        else if (ch == ')') {
            level--;
            inodeb  = cparent;               /* this internal node may get :length or #label */
            cparent = nodes[cparent].parent; /* move back up the tree */
        }

        /* Colon: branch length for the node stored in inodeb */
        else if (ch == ':') {
            haslength = 1;
            if (fscanf(treeFile, "%lf", &nodes[inodeb].time) != 1) {
                error("Error reading branch length in tree file");
            }
        }

        /* # or $: branch label for the node in inodeb */
        else if (ch == '#' || ch == '$') {
            int lab;
            haslabel = 1;
            if (fscanf(treeFile, "%d", &lab) != 1) {
                error("Error reading branch label in tree file");
            }
            nodes[inodeb].label = lab;
        }

        /* Comma just separates siblings */
        else if (ch == ',') {
            /* nothing to do here */
        }

        /* ';' should only appear after the root is closed (level == 0) */
        else if (ch == ';' && level != 0) {
            error("; in treefile before all parentheses are closed");
        }

        /* Otherwise we are reading a taxon name or a taxon number */
        else {
            /* Read the first two characters (we already consumed ch) */
            line[0] = (char) ch;
            line[1] = (char) fgetc(treeFile);

            /* If number of sequences < 10 and first two chars are digits,
               assume a pure number for the species index. */
            if (SeqData.No_seq < 10 && isdigit((unsigned char)line[0]) &&
                isdigit((unsigned char)line[1])) {
                ungetc(line[1], treeFile);
                line[1] = 0;
            } else {
                /* Otherwise read until a delimiter or EOL/EOF */
                for (i = 1; i < lline; ) {
                    if (strchr(delimiters, line[i]) ||
                        line[i] == (char)EOF ||
                        line[i] == '\n') {
                        ungetc(line[i], treeFile);
                        line[i] = 0;
                        break;
                    }
                    line[++i] = (char) fgetc(treeFile);
                }
            }

            /* Trim trailing non-graphic characters */
            for (j = i - 1; j > 0; j--)
                if (!isgraph((unsigned char)line[j]))
                    line[j] = 0;

            /* Decide if this is a species name or a number */
            for (i = 0, hasname = 0; line[i]; i++) {
                if (!isdigit((unsigned char)line[i])) {
                    hasname = 1;
                    break;
                }
            }

            if (hasname) {
                /* look up name in SeqData.name[] */
                for (i = 0; i < SeqData.No_seq; i++)
                    if (!strcmp(line, SeqData.name[i]))
                        break;
                if ((cnode = i) == SeqData.No_seq) {
                    printf("\nError: species %s not found in sequence data.\n", line);
                    exit(-1);
                }
            } else {
                /* numeric taxon index (1-based in file, 0-based internally) */
                sscanf(line, "%d", &cnode);
                cnode--;
                if (cnode < 0 || cnode >= SeqData.No_seq) {
                    printf("\nError in tree file: species %d not in data.\n", cnode + 1);
                    exit(-1);
                }
            }

            /* Attach this leaf to its parent */
            nodes[cnode].parent = cparent;
            nodes[cnode].seq_no = cnode + 1;

            nodes[cparent].child[nodes[cparent].no_child++] = cnode;

            tree.edges[tree.no_edge][0] = cparent;
            tree.edges[tree.no_edge][1] = cnode;
            nodes[cnode].ibranch        = tree.no_edge++;

            /* mark that this sequence has been seen exactly once */
            check[inodeb = cnode]++;
        }

        if (level == 0)
            break;
    }

    if (popline)
        fgets(line, 254, treeFile);

    if (tree.no_node != tree.no_edge + 1)
        printf("\nnnode %6d != nbranch %6d + 1\n", tree.no_node, tree.no_edge);

    for (i = 0; i < SeqData.No_seq; i++)
        if (check[i] != 1)
            return -1;

    if (tree.no_edge > 2 * SeqData.No_seq - 2) {
        printf("\nnbranch %d", tree.no_edge);
        error("too many branches in tree?");
    }

    *length_label = haslength + 2 * haslabel;

    /* ------------------------------------------------------------
     *  Map node labels (#1, #2) to gamma classes
     *     class 0: background
     *     class 1: foreground #1
     *     class 2: foreground #2
     * ------------------------------------------------------------ */
    if (SeqData.Model == 3 || SeqData.Model == 4) {
        int fg1_count = 0;
        int fg2_count = 0;

        /* reset for actually used branches */
        for (i = 0; i < tree.no_edge; i++) {
            is_foreground_branch[i] = 0;
            branch_gamma_class[i]   = 0;
        }

        for (i = 0; i < tree.no_node; i++) {
            int b   = nodes[i].ibranch;
            int lab = nodes[i].label;

            if (b < 0 || b >= tree.no_edge)
                continue; /* root or invalid */

            if (lab == 1) {
                branch_gamma_class[b]   = 1;
                is_foreground_branch[b] = 1;
                fg1_count++;
            }
            else if (lab == 2) {
                branch_gamma_class[b]   = 2;
                is_foreground_branch[b] = 1;
                fg2_count++;
            }
            else if (lab > 2) {
                fprintf(stderr,
                        "\n[WARN] Node %d has unsupported label #%d "
                        "(only #1 and #2 are recognized for model 3).\n",
                        i, lab);
            }
        }

        num_fg1_branches = fg1_count;
        num_fg2_branches = fg2_count;

        printf("\n[DEBUG] Foreground branches (model %d):\n", SeqData.Model);
        for (i = 0; i < tree.no_edge; i++) {
            if (branch_gamma_class[i] > 0) {
                int parent = tree.edges[i][0];
                int child  = tree.edges[i][1];
                int cls    = branch_gamma_class[i];
                printf("  branch %2d: class=%d  parent=%d, child=%d",
                       i, cls, parent, child);
                if (nodes[child].seq_no > 0)
                    printf(" (leaf %s)", SeqData.name[nodes[child].seq_no - 1]);
                printf("\n");
            }
        }
        printf("[DEBUG] #1 branches = %d, #2 branches = %d\n",
               fg1_count, fg2_count);
    }

    /* ------------------------------------------------------------
     *  Print ASCII phylogram to the terminal
     * ------------------------------------------------------------ */
    *length_label = haslength + 2 * haslabel;

    /* OPTIONAL: debug print of tree and foreground branches */
    DebugPrintAsciiTree(stdout);

    return 0;
}

/**
  * Function that allocates memory for the node array of the tree structure
 **/
int new_tree(int No_nodes)
{
    int i;
    i = (No_nodes) * sizeof(struct node);
    if ((nodes = (struct node *) malloc(i)) == NULL)
        error("out of memory in allocating tree nodes");
    return (0);
}

/**
  * Identical codon sequence patterns are detected at each site
  * to avoid redundant calculations in order to reduce the amount of calculation
  */
int FindPattern()
{
    int i, j, k, n, seq_i;
    int *tempseq[NS];    // Representative sequence of each pattern
    int *seq[NS];        // Actual codon sequence
    int patternmatch;    // Number of perfectly matched sequences
    int patternfound;    // Flag whether it matches an existing pattern
    int patternat;       // Index of the matched pattern
    int Npattern = 0;    // Number of patterns detected so far
    int patterncounts[SEQLENGTH]; // Number of occurrences of each pattern
    int mappattern[SEQLENGTH];    // Which pattern does each position belong to
    int missing_data;    // Flag whether missing data (-2) is included

    n = SeqData.No_codon;  // Number of codons

    // Initialization
    for (i = 0; i < n; i++)
        patterncounts[i] = 0;

    for (i = 0; i < SeqData.No_seq; i++) {
        tempseq[i] = (int *) malloc(n * sizeof(int));
        seq[i] = SeqData.Codonseq[i];
    }

    for (i = 0; i < n; i++) {
        patternfound = 0;
        missing_data = 0;

        for (seq_i = 0; seq_i < SeqData.No_seq; seq_i++) {
            if (seq[seq_i][i] == -2) {
                missing_data = 1;
                break;
            }
        }

        if (missing_data == 0) {
            for (j = 0; j < Npattern; j++) {
                patternmatch = 0;
                for (k = 0; k < SeqData.No_seq; k++) {
                    if (tempseq[k][j] == seq[k][i])
                        patternmatch++;
                }
                if (patternmatch == SeqData.No_seq) {
                    // Exact match → Matches an existing pattern
                    patternfound = 1;
                    patternat = j;
                    break;
                }
            }

            if (patternfound) {
                // Register to an existing pattern
                patterncounts[patternat]++;
                mappattern[i] = patternat;
            } else {
                // Register as a new pattern
                for (k = 0; k < SeqData.No_seq; k++)
                    tempseq[k][Npattern] = seq[k][i];
                mappattern[i] = Npattern;
                patterncounts[Npattern] = 1;
                Npattern++;
            }
        } else {
            // Mark -1 if it contains missing data
            mappattern[i] = -1;
        }
    }

    // Save the number of patterns found
    SeqData.CNo_pattern = Npattern;
    printf("\nFound %d patterns in the coding region.", Npattern);

    // Reduce the codon sequence to just the pattern
    for (i = 0; i < SeqData.No_seq; i++) {
        memcpy(seq[i], tempseq[i], Npattern * sizeof(int));
        seq[i] = (int *) realloc(seq[i], Npattern * sizeof(int));
    }

    // Release the temporary area
    for (i = 0; i < SeqData.No_seq; i++)
        free(tempseq[i]);

    // Save pattern frequency and mapping
    for (i = 0; i < n; i++) {
        SeqData.CPatternCounts[i] = patterncounts[i];
        SeqData.CMapPattern[i] = mappattern[i];
    }

    return 0;
}

/**
  * Initialize the X vector to be optimized
  * It returns the number of free parameters.
  */
int InitializeX(double x[])
{
    int i, j, k;

    /* 1. initialize branch length parameters x[1..no_edge]
    */
    for (i = 0; i < tree.no_edge; i++) {
        int child = tree.edges[i][1];
        double t  = nodes[child].time;

        if (t <= 0.0) {
            t = 0.05;
        }
        x[i + 1] = t;
    }

    /* 2. initial kappa */
    i = tree.no_edge;

    if (SeqData.kappa <= 0.0)
        SeqData.kappa = 4.0;

    x[++i] = SeqData.kappa;

    /* 3. initial base-frequency parameters (unconstrained) */
    {
        double b0 = 0.0, b1 = 0.0, b2 = 0.0;

        if (SeqData.basepi[0] > 0.0 &&
            SeqData.basepi[1] > 0.0 &&
            SeqData.basepi[2] > 0.0 &&
            SeqData.basepi[3] > 0.0) {

            b0 = log(SeqData.basepi[0] / SeqData.basepi[3]);
            b1 = log(SeqData.basepi[1] / SeqData.basepi[3]);
            b2 = log(SeqData.basepi[2] / SeqData.basepi[3]);
        }

        x[++i] = b0; // for T
        x[++i] = b1; // for C
        x[++i] = b2; // for A
    }

    /* 4. model-specific parameters */

    if (SeqData.Model == 1 || SeqData.Model == 2) {
        int n_cat;
        for (j = 0; j < SeqData.omegaCats; j++)
            x[++i] = (j == 0 ? 0.8 : 2.0);

        for (j = 0; j < SeqData.gammaCats[0]; j++)
            x[++i] = (j == 0 ? 0.8 : 2.0);

        n_cat = SeqData.omegaCats * SeqData.gammaCats[0];
        if (n_cat > 1) {
            for (k = 0; k < n_cat - 1; k++)
                x[++i] = 0.5;
        }
    }
    else if (SeqData.Model == 3) {
        /* ---- branch model alternative (Model=3) ----
        */

        /* initial omega */
        double om0 = SeqData.omega[0];
        if (om0 <= 0.0) om0 = 0.9;
        x[++i] = om0;

        /* initial gamma0 (background)：
        */
        double g0 = SeqData.gamma[0][0];
        if (g0 <= 0.0) g0 = 0.8;
        x[++i] = g0;   /* gamma_0 */

        /* initial gamma1 and gamma2 gamma0 estimated by null model (model4) */
        x[++i] = g0;   /* gamma_1 (foreground #1) */
        x[++i] = g0;   /* gamma_2 (foreground #2) */
    }
    else if (SeqData.Model == 4) {
        /* ---- branch model null (Model=4) ---- */

        double om0 = SeqData.omega[0];
        if (om0 <= 0.0) om0 = 0.9;
        x[++i] = om0;

        double g0 = SeqData.gamma[0][0];
        if (g0 <= 0.0) g0 = 0.8;
        x[++i] = g0;   /* gamma (all branches) */
    }

    SeqData.np = i;
    return i;
}

/**
  * Calculates the quadratic penalty for the function
 **/
double quadpenalty(double x[], int n)
{
    int i, j;

    // Defining parameter boundaries
    double tb[] = { 1e-4, 5.0};                    // branch length
    double omegabLessThan1[] = { 1e-4, 0.9999 };    // omega_0
    double omegabGreaterThan1[] = { 1.0001, 25.0 };    // omega_1
    double omegab[] = { 1e-4, 25.0 };                  // omega (for branch model)
    double kappab[] = { 1e-4, 6.50 };                  // kappa
    double gammabLessThan1[] = { 1e-4, 1 };         // gamma_0
    double gammabGreaterThan1[] = { 1.0001, 25.0 };    // gamma_1
    double gammab[] = { 1e-4, 25.0 };                  // gamma (for branch model)
    double pb[] = { 0, 0.9999 };                    // Mixture ratio p

    double xb[NP][2];
    double penalty = 0.0;

    // branch length
    for (i = 1; i <= tree.no_edge; i++)
        for (j = 0; j < 2; j++)
            xb[i][j] = tb[j];

    // kappa
    xb[i][0] = kappab[0];
    xb[i][1] = kappab[1];

    // omega / gamma / mixture ratio
    if (SeqData.Model == 1) { // site model alternative hypothesis
        // omega_0, omega_1
        xb[++i][0] = omegabLessThan1[0]; xb[i][1] = omegabLessThan1[1];
        xb[++i][0] = omegabGreaterThan1[0]; xb[i][1] = omegabGreaterThan1[1];
        // gamma_0, gamma_1
        xb[++i][0] = gammabLessThan1[0]; xb[i][1] = gammabLessThan1[1];
        xb[++i][0] = gammabGreaterThan1[0]; xb[i][1] = gammabGreaterThan1[1];
        // mixture ratio (p1, p2, p3)
        for (j = 0; j < 3; j++) {
            xb[++i][0] = pb[0];
            xb[i][1] = pb[1];
        }
    }
    else if (SeqData.Model == 2) { // site model null hypothesis
        // omega_0, omega_1
        xb[++i][0] = omegabLessThan1[0]; xb[i][1] = omegabLessThan1[1];
        xb[++i][0] = omegabGreaterThan1[0]; xb[i][1] = omegabGreaterThan1[1];
        // only gamma_0
        xb[++i][0] = gammabLessThan1[0]; xb[i][1] = gammabLessThan1[1];
        // mixture ratio (p)
        xb[++i][0] = pb[0]; xb[i][1] = pb[1];
    }
    else if (SeqData.Model == 3) { // branch model alternative hypothesis
        // omega
        xb[++i][0] = omegab[0]; xb[i][1] = omegab[1];
        // gamma_0 (background), gamma_1 (foreground #1), gamma_2 (foreground #2)
        xb[++i][0] = gammab[0]; xb[i][1] = gammab[1];
        xb[++i][0] = gammab[0]; xb[i][1] = gammab[1];
        xb[++i][0] = gammab[0]; xb[i][1] = gammab[1];
    }
    else if (SeqData.Model == 4) { // branch model null hypothesis
        // omega
        xb[++i][0] = omegab[0]; xb[i][1] = omegab[1];
        // gamma
        xb[++i][0] = gammab[0]; xb[i][1] = gammab[1];
    }

    return penalty;
}

/**
  * Function that smoothly corrects values ​​when they fall below a lower limit
 **/
double transformXLower(double x, double lowerbound)
{
    if (x < lowerbound) return lowerbound;
    return x;
}

double transformXUpper(double x, double upperbound)
{
    if (x > upperbound) return upperbound;
    return x;
}

/**
  * Objective function combining log-likelihood and penalty terms
 **/
double func(double x[], int n)
{
    int i, j, k;
    double basepi[3];
    double probGamma[3];
    double lkl;

    tree.LogL = 0;

    eigen_error_flag = 0;

    // --- 1. branch length setting ---
    for (i = 0; i < tree.no_edge; i++) {
        int child = tree.edges[i][1];  /* child node of edge i */
        double t  = transformXLower(x[i + 1], 1e-4);
        t         = transformXUpper(t, 5.0);
        nodes[child].time = t;
    }

    // --- 2. kappa setting ---
    k = tree.no_edge + 1;
    SeqData.kappa = transformXLower(x[k], 1e-4);
    SeqData.kappa = transformXUpper(SeqData.kappa, 25.000);

    // --- 3. Base frequency setting ---
    for (j=0; j<3; j++)
        basepi[j] = x[++k];

    SeqData.basepi[0]=exp(basepi[0])/(exp(basepi[0])+exp(basepi[1])+exp(basepi[2])+1);
    SeqData.basepi[1]=exp(basepi[1])/(exp(basepi[0])+exp(basepi[1])+exp(basepi[2])+1);
    SeqData.basepi[2]=exp(basepi[2])/(exp(basepi[0])+exp(basepi[1])+exp(basepi[2])+1);
    SeqData.basepi[3]=1-SeqData.basepi[0]-SeqData.basepi[1]-SeqData.basepi[2];

    // Codon frequency update
    calCodonFreqs();

    // --- 4. // Read parameters for each model ---
    if (SeqData.Model == 1 || SeqData.Model == 2) {
        SeqData.omega[0] = transformXLower(x[++k], 1e-4);
        SeqData.omega[0] = transformXUpper(SeqData.omega[0], .9999);
        SeqData.omega[1] = transformXLower(x[++k], 1.0001);
        SeqData.omega[1] = transformXUpper(SeqData.omega[1], 25.000);
    }
    else if (SeqData.Model == 3 || SeqData.Model == 4) {
        SeqData.omega[0] = transformXLower(x[++k], 1e-4);
        SeqData.omega[0] = transformXUpper(SeqData.omega[0], 25.000);
    }

    if (SeqData.Model == 1) {
        SeqData.gamma[0][0] = transformXLower(x[++k], 1e-4);
        SeqData.gamma[0][0] = transformXUpper(SeqData.gamma[0][0], .9999);
        SeqData.gamma[0][1] = transformXLower(x[++k], 1.0001);
        SeqData.gamma[0][1] = transformXUpper(SeqData.gamma[0][1], 25.000);

        for (j = 0; j < 3; j++)
            probGamma[j] = x[++k];

        double sum_probGamma = exp(probGamma[0]) + exp(probGamma[1]) + exp(probGamma[2]) + 1;
        SeqData.probGamma[0] = exp(probGamma[0]) / sum_probGamma;
        SeqData.probGamma[1] = exp(probGamma[1]) / sum_probGamma;
        SeqData.probGamma[2] = exp(probGamma[2]) / sum_probGamma;
        SeqData.probGamma[3] = 1 - SeqData.probGamma[0] - SeqData.probGamma[1] - SeqData.probGamma[2];
    }
    else if (SeqData.Model == 2) {
        SeqData.gamma[0][0] = transformXLower(x[++k], 1e-4);
        SeqData.gamma[0][0] = transformXUpper(SeqData.gamma[0][0], .9999);

        SeqData.probGamma[0] = transformXLower(x[++k], 1e-4);
        SeqData.probGamma[0] = transformXUpper(SeqData.probGamma[0], .9999);
        SeqData.probGamma[1] = 1 - SeqData.probGamma[0];
    }
    else if (SeqData.Model == 3) {
        SeqData.gamma[0][0] = transformXLower(x[++k], 1e-4);
        SeqData.gamma[0][0] = transformXUpper(SeqData.gamma[0][0], 25.000);

        SeqData.gamma[0][1] = transformXLower(x[++k], 1e-4);
        SeqData.gamma[0][1] = transformXUpper(SeqData.gamma[0][1], 25.000);

        SeqData.gamma[0][2] = transformXLower(x[++k], 1e-4);
        SeqData.gamma[0][2] = transformXUpper(SeqData.gamma[0][2], 25.000);
    }
    else if (SeqData.Model == 4) {
        SeqData.gamma[0][0] = transformXLower(x[++k], 1e-4);
        SeqData.gamma[0][0] = transformXUpper(SeqData.gamma[0][0], 25.000);
    }

    // --- 5. Likelihood calculation (switching functions depending on the model) ---
    if (SeqData.Model == 3 || SeqData.Model == 4)
        lkl = return_fixGammaP(tree.root);
    else
        lkl = return_gammaPartsP(tree.root);

    // --- 6. eigen error or NaN/inf check---
    if (eigen_error_flag || !isfinite(lkl)) {
        return 1e10;
    }

    tree.LogL = lkl;

    // --- 7. Return the objective function value ---
    return (-1.0 * tree.LogL + quadpenalty(x, n));
}


/**
  * Function that loads the specified AA pair and initializes the AAPartition matrix
  * Stores 0 for amino acid pairs listed in the ZpartAA file, and 1 for others
  */
int GetAAPart (int gamma_i, char* ZpartAAFName)
{
   printf("DEBUG: Trying to open AAPartition file: [%s]\n", ZpartAAFName);
   char line[MAXLINE];
   FILE *fin;
   int n1step=0, i,j,k, iaa,jaa, npair, naa=20;

   // Initialize AAPartition with all 1s (unspecified)
   for(i=0,n1step=0; i<naa; i++)
	   for(j=0; j<i; j++) {
		   SeqData.AAPartition[gamma_i][i*(i-1)/2+j]=1;
		   n1step++;
	   }

   fin=gfopen(ZpartAAFName, "r");

   fgets (line, MAXLINE, fin);

   for (j=0, npair=0; j < MAXLINE-1 && line[j] && line[j]!='\n'; j++) {
       iaa = line[j];
	   if (!isalpha(iaa)) continue;
       jaa = line[++j];
	   if (!isalpha(jaa)) printf("\nerr jaa");
       npair++;

	   iaa = AA2Code((char)iaa);
	   jaa = AA2Code((char)jaa);

	   if(iaa<0||iaa>19||jaa<0||jaa>19)
		   printf("\naa not found");

	   if (iaa < jaa)  {
		   k = jaa; jaa = iaa; iaa = k;
	   }

       printf ("|%c%c (%2d,%2d)| ", AACodes[iaa], AACodes[jaa], iaa, jaa);

       if (iaa == jaa) printf ("\nThis pair has no effect.");

	   if (!SeqData.AAPartition[gamma_i][iaa*(iaa-1)/2+jaa])
		   printf("\nThis pair specified?");

       // Set the specified AA pair to 0
       SeqData.AAPartition[gamma_i][iaa*(iaa-1)/2+jaa] = 0;
   }

   for (i = 0; i < naa; i++) {
	   for(j = 0; j < i; j++)
	       printf ("%3d", SeqData.AAPartition[gamma_i][i*(i-1)/2 + j]);
   }

   return (0);
}

/*
 * Converts one-letter amino acid codes (A, C, D...) into numbers from 0 to 19
 * If not found, returns -1
 */
int AA2Code (char b)
{
   int i, n = 20;

   for (i = 0; i < n; i++)
	   if (b == AACodes[i])
		   return (i);

   printf ("\nerr: strange character '%c' ", b);
   return (-1);
}

/**
 * Read the options file (e.g. evoradical.dat) and initialize the SeqData structure
 */
int ReadOptionFile(char* optionFName) {
	FILE *optionfile;
	char line[MAXLINE], parameter[MAXSTRING], value[MAXSTRING];
	char *p;
	const char *szSeparator = ", ";
	int i;

	optionfile = gfopen(optionFName, "r+");

	// Model initialization
	SeqData.Model = 0;

	while (1) {
		strcpy(line, "");
		p = fgets(line, MAXLINE, optionfile);
		if (p == NULL)
			break;

		if ((*p == '#') || (*p == '\n') || (line[0] == '#') || strstr(line, "=") == NULL)
			continue;

		// Extracting parameter names and values
		sscanf(line, "%s = %s", parameter, value);

		if (strcmp(parameter, "seqfile") == 0)
			strcpy(SeqData.seqFName, value);

		else if (strcmp(parameter, "Model") == 0)
			SeqData.Model = atoi(value);

		else if (strcmp(parameter, "treefile") == 0)
			strcpy(SeqData.treeFName, value);

		else if (strcmp(parameter, "outfile") == 0)
			strcpy(SeqData.outFName, value);

		else if (strcmp(parameter, "GammaFiles") == 0)
			String2Strs(line, SeqData.partFNames, 1, szSeparator);
	}

	if (SeqData.Model == 0) {
		printf("\nModel is not specified\n");
		exit(-1);
	}

	if (SeqData.Model >= 1 && SeqData.Model <= 4)
		SeqData.no_parts = 1;
	else {
		printf("\nUnsupported model number: %d", SeqData.Model);
		exit(-1);
	}

	// Setting the number of gamma categories (varies by model)
	for (i = 0; i < SeqData.no_parts; i++) {
		if (SeqData.Model == 1 || SeqData.Model == 3)  // Alternative hypothesis
			SeqData.gammaCats[i] = 2;
		else  // Null hypothesis
			SeqData.gammaCats[i] = 1;
	}

	// Setting the number of omega categories
	if (SeqData.Model == 1 || SeqData.Model == 2) // Site model
		SeqData.omegaCats = 2;
	else if (SeqData.Model == 3 || SeqData.Model == 4) // Branch model
		SeqData.omegaCats = 1;

	fclose(optionfile);

	printf("\n%s read in successfully:\n", optionFName);
	printf("Number of categories for omega = %d\n", SeqData.omegaCats);
	printf("Number of partitions for gamma = %d\n", SeqData.no_parts);
	for (i = 0; i < SeqData.no_parts; i++)
		printf("Partition %d of gamma: Number of category = %d\n", i + 1, SeqData.gammaCats[i]);

	return 0;
}

int String2Strs(char* line, char* outputStrings[], int no_parts, const char *szSeparator) {
	char *p;
	char value[MAXSTRING];
	const char *strTemp;
	int i;

	p = line;
	while (*p != '=' && *p != '\0')  // Continue until '=' is found
		p++;
	if (*p == '=') p++; // Start after '='

	// Assign to value
	strcpy(value, p);

	// Trim line breaks, carriage returns, and unnecessary symbols
	char *q = strpbrk(value, "\r\n>");
	if (q) *q = '\0';

	// Get first token
	strTemp = strtok(value, szSeparator);
	outputStrings[0] = calloc(strlen(strTemp) + 1, sizeof(char));
	strcpy(outputStrings[0], strTemp);

	// Get remaining tokens
	for (i = 1; i < no_parts; i++) {
		strTemp = strtok(NULL, szSeparator);
		if (strTemp == NULL) {
			printf("\nERROR: Invalid GammaFiles line: %s", line);
			exit(-1);
		}
		outputStrings[i] = calloc(strlen(strTemp) + 1, sizeof(char));
		strcpy(outputStrings[i], strTemp);
	}

	return 0;
}

int String2Ints(char* line, int outputInts[], int no_parts, const char *szSeparator) {
	char *p;
	char value[MAXSTRING];
	const char *strTemp;
	int i;

	// Continue until '=' is found
	p = line;
	while (*p != '=' && *p != '\0')
		p++;
	if (*p == '=') p++;  // Start after '='

	strcpy(value, p);

	// Remove line breaks and garbage
	char *q = strpbrk(value, "\r\n>");
	if (q) *q = '\0';

	// Get first token
	strTemp = strtok(value, szSeparator);
	if (strTemp == NULL) {
		printf("\nERROR: No value found in line: %s", line);
		exit(-1);
	}
	outputInts[0] = atoi(strTemp);

	// Get remaining tokens
	for (i = 1; i < no_parts; i++) {
		strTemp = strtok(NULL, szSeparator);
		if (strTemp == NULL) {
			printf("\nERROR: Not enough values in line: %s", line);
			exit(-1);
		}
		outputInts[i] = atoi(strTemp);
	}

	return 0;
}

int String2Doubles(char* line, double outputDoubles[], int no_parts, const char *szSeparator) {
	char *p;
	char value[MAXSTRING];
	const char *strTemp;
	int i;

	// Find '=' and start reading from the next string
	p = line;
	while (*p != '=' && *p != '\0')
		p++;
	if (*p == '=') p++;

	strcpy(value, p);

	// Remove line breaks and unnecessary symbols
	char *q = strpbrk(value, "\r\n>");
	if (q) *q = '\0';

	// Get first token
	strTemp = strtok(value, szSeparator);
	if (strTemp == NULL) {
		printf("\nERROR: No value found in line: %s", line);
		exit(-1);
	}
	outputDoubles[0] = atof(strTemp);

	// Get remaining tokens
	for (i = 1; i < no_parts; i++) {
		strTemp = strtok(NULL, szSeparator);
		if (strTemp == NULL) {
			printf("\nERROR: Not enough values in line: %s", line);
			exit(-1);
		}
		outputDoubles[i] = atof(strTemp);
	}

	return 0;
}