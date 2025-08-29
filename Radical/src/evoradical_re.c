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
 *  Add branch model and remove some models by Hayate Takeuchi (Final Edit: April 11 2025)
 *************************************/

#include "evoradical_re.h"

// --- Global variables for the Branch model ---
int is_foreground_branch[MAXBRANCHES];  // foreground branch
double gamma0, gamma1;                  // background, foreground gamma

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
    double *pi = SeqData.pi; // Frequency of 61 codons

    // Initialization
    for (i = 0; i < n * n; i++)
        Q[i] = 0;

    // For all codon pairs
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

            // Determine if there is a single base mutation (1 step)
            for (k = 0, ndiff = 0; k < 3; k++) {
                if (from[k] != to[k]) {
                    ndiff++;
                    pos = k;
                }
            }
            if (ndiff != 1) continue; // Skip any changes that are not one step

            Q[i * n + j] = 1.0;

            // If transition, it is kappa times
            if ((from[pos] + to[pos] == 1) || (from[pos] + to[pos] == 5))
                Q[i * n + j] = kappa;

            // Scaled by the frequency of target codon j
            Q[i * n + j] *= pi[j];

            // Obtain amino acid changes
            aa1 = GeneticCode[0][codon1];
            aa2 = GeneticCode[0][codon2];

            // Apply omega or gamma if nonsynonymous substitution
            if (aa1 != aa2) {
                int is_radical = (AAPartition[aa1 * (aa1 - 1) / 2 + aa2] == 0
                               || AAPartition[aa2 * (aa2 - 1) / 2 + aa1] == 0);

                switch (SeqData.Model) {
                    case 1: // site model alternative hypothesis
                        if (is_radical)
                            Q[i * n + j] *= gamma; // Radical substitution: q_ij = pi_j * kappa  * gamma (* omega=1)
                        else
                            Q[i * n + j] *= omega; // Conservative substitutions: q_ij = pi_j * kappa  * omega (* gamma=1)
                        break;

                    case 2: // Site model null hypothesis (fixed gamma_0)
                        Q[i * n + j] *= omega;  // Apply only omega (gamma_0 is embedded later)
                        break;

                    case 3: // branch model alternative hypothesis
                        // Different gamma for foreground and background branches
                        if (is_radical) {
                            double gamma_branch = is_foreground_branch[branch_i] ? gamma1 : gamma0;
                            Q[i * n + j] *= gamma_branch;
                        }
                        Q[i * n + j] *= omega; // omega is common to all branches
                        break;

                    case 4: // branch model null hypothesis
                        if (is_radical)
                            Q[i * n + j] *= gamma; // gamma is common to all branches
                        else
                            Q[i * n + j] *= omega; // omega is common to all branches
                        break;
                }
            }

            // For reversibility, substitute the same value for Q[j][i]
            Q[j * n + i] = Q[i * n + j];
        }
    }

    // Put -∑ Q[i][j] on the diagonal (so that the row sum is 0)
    for (i = 0; i < n; i++)
        Q[i * n + i] = -sum(Q + i * n, n);

    // Perform eigenvalue expansion and decompose the Q matrix into a form that can be used for exp(Qt)
    if (eigen(1, Q, n, Root, ri, U, V, space + n))
        error("eigen function Qmatrix err.");

    xtoy(U, V, n * n);      // U = V^T
    matinv(V, n, n, space); // V = V^{-1}

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
    Qmatrix(Root, U, V, SeqData.kappa, omega, gamma, AAPartition, PMat, n, -1);

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

void return_p(int node_i)
{
    int i, j, k;
    int child_i, pattern_i, temp;
    int n = 61;
    int no_pattern = SeqData.CNo_pattern;
    double t, lkl;

    // If the child node is an internal node (not a leaf node), it is called recursively
    for (i = 0; i < nodes[node_i].no_child; i++) {
        if (nodes[nodes[node_i].child[i]].no_child > 0)
            return_p(nodes[node_i].child[i]);
    }

    // For each child node, propagate the likelihood to the parent node
    for (i = 0; i < nodes[node_i].no_child; i++) {
        child_i = nodes[node_i].child[i];
        t = nodes[child_i].time;

        // In model 3, Q is constructed using different gammas for each branch
        if (SeqData.Model == 3) {
            int branch_i = nodes[child_i].ibranch;
            double gamma_branch = is_foreground_branch[branch_i] ? gamma1 : gamma0;

            // Create the Q matrix for branch and perform the eigendecomposition
            Qmatrix(Root, U, V, SeqData.kappa, SeqData.omega[0],
                    gamma_branch, SeqData.AAPartition[0],
                    PMat, n, branch_i);
        }
        // In model 4, Q is constructed once using a common gamma across all branches
        else if (SeqData.Model == 4) {
            if (child_i == nodes[node_i].child[0]) {
                Qmatrix(Root, U, V, SeqData.kappa, SeqData.omega[0],
                        gamma0, SeqData.AAPartition[0],
                        PMat, n, -1);
            }
        }

        // For model==1,2,4 PMat is already computed
        PMatUVRoot(PMat, t, n, U, V, Root);

        for (pattern_i = 0; pattern_i < no_pattern; pattern_i++) {
            if (nodes[child_i].no_child < 1) {
                for (j = 0; j < n; j++) {
                    temp = SeqData.Codonseq[child_i][pattern_i];
                    if (temp != -2)
                        nodes[node_i].likelihood[n * pattern_i + j] *= PMat[j * n + temp];
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

    // Error except for models 1 and 2
    if (!(SeqData.Model == 1 || SeqData.Model == 2)) {
        error("return_gammaPartsP() is only available for model==1,2 (site model).");
    }

    // Calculate likelihood for each omega × gamma category
    for (i = 0; i < SeqData.omegaCats; i++) {
        for (j = 0; j < SeqData.gammaCats[0]; j++) {
            InitializeLkl_site(SeqData.omega[i], SeqData.gamma[0][j], SeqData.AAPartition[0]);

            return_p(tree.root);  // felsenstein algorithm

            // Store the likelihood (numerator) for each site pattern in the posterior probability storage array
            for (l = 0; l < no_pattern; l++) {
                prob = 0.0;
                for (k = 0; k < n; k++)
                    prob += nodes[tree.root].likelihood[l * n + k] * pi[k];

                SeqData.postdistGamma[
                    l * SeqData.omegaCats * SeqData.gammaCats[0]
                    + i * SeqData.gammaCats[0] + j
                ] = prob * SeqData.probGamma[i * SeqData.gammaCats[0] + j];
            }
        }
    }

    // For each pattern, calculate the log-likelihood by summing the categories
    for (j = 0; j < no_pattern; j++) {
        prob = 0.0;
        for (i = 0; i < SeqData.omegaCats * SeqData.gammaCats[0]; i++)
            prob += SeqData.postdistGamma[j * SeqData.omegaCats * SeqData.gammaCats[0] + i];

        if (prob < eps) {
            printf("⚠️ pattern %d: probability is zero, applying epsilon.\n", j);
            prob = eps;
        }

        lkl += log(prob) * patterncounts[j];
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

    // Model check: model==3 or model==4 only
    if (!(SeqData.Model == 3 || SeqData.Model == 4)) {
        error("return_fixGammaP() is only available for model==3 or 4 (branch models).");
    }

    // Read gamma_0, gamma_1
    double gamma_bg = SeqData.gamma[0][0];  // For Common (model==4) or background (model==3)
    double gamma_fg = (SeqData.Model == 3 && SeqData.gammaCats[0] > 1)
                      ? SeqData.gamma[0][1] // For foreground (model==3)
                      : gamma_bg;

    // Initialization for branch models
    InitializeLkl_branch(SeqData.omega[0], gamma_bg, gamma_fg, SeqData.AAPartition[0]);

    // Likelihood propagation to the root using the Felsenstein algorithm
    return_p(tree.root);

    // Calculate the log-likelihood for each pattern
    for (l = 0; l < no_pattern; l++) {
        prob = 0.0;
        for (k = 0; k < n; k++) {
            prob += nodes[tree.root].likelihood[l * n + k] * pi[k];
        }

        if (prob < eps) {
            printf("⚠️ pattern %d: probability is zero, applying epsilon.\n", l);
            prob = eps;
        }

        lkl += log(prob) * patterncounts[l];
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
    no_param = InitializeX(x);

    for (i = 0; i < 3; i++)
        dfpmin(x, no_param, GTOL, &iter, &fret, func);

    printf("\nThe log likelihood of the tree is: %f", tree.LogL);

    fp = gfopen(SeqData.outFName, "w");
    printX(x, fp);
    fclose(fp);

    double elapsed_time = prn_time();
    printf("\nRun Time: %.2f s\n", elapsed_time);

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
    else if (SeqData.Model == 3)  // branch model alt
        df_model = 4 + num_branch;
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
    strftime(start_str, sizeof(start_str), "%Y-%m-%d %H:%M:%S", localtime(&tk.begin_time));
    fprintf(fp, "Start time: %s\n", start_str);


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
            nodes[i].time = transformXLower(x[i+1], 1e-4);
            fprintf(fp, "Branch %d: %f\n", i, nodes[i].time);
        }

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
            nodes[i].time = transformXLower(x[i+1], 1e-4);
            fprintf(fp, "Branch %d: %f\n", i, nodes[i].time);
        }

        fprintf(fp, "\n--- PARAMETERS ---\n");

        // Kappa output
        fprintf(fp, "Kappa: %f\n", SeqData.kappa);

        // Omega output
        fprintf(fp, "Omega: %f\n", SeqData.omega[0]);

        // Gamma output（if model==3, gamma of foreground and background）
        if (SeqData.Model == 3) {
            fprintf(fp, "Gamma (background): %f\n", SeqData.gamma[0][0]);
            fprintf(fp, "Gamma (foreground): %f\n", SeqData.gamma[0][1]);
        } else {
            fprintf(fp, "Gamma (all branches): %f\n", SeqData.gamma[0][0]);
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
    // Saves current CPU time and wall clock time for record keeping and delta measurement
    tk.begin_clock = tk.save_clock = clock();
    tk.begin_time  = tk.save_time  = time(NULL);
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

/**
  * Function to read user-defined phylogenetic trees (in Newick format)
  * Based on the structure borrowed from PAML's codeml.c, 
  * support for the foreground branch (#1) was added.
  */
int ReadTree(FILE * treeFile, int *length_label, int popline)
{
    int cnode, cparent;
    int inodeb;
    int i, j, level, hasname, haslength, haslabel, ch, lline;
    char check[NS], line[255], delimiters[] = "(),:#$;";

    cparent = -1;
    inodeb = 0;
    level = 0;
    haslength = 0;
    haslabel = 0;
    ch = ' ';
    lline = 255;

    *length_label = 0;
    tree.no_node = SeqData.No_seq;
    tree.no_edge = 0;

    // Initializing the node structure
    for (i = 0; i < 2 * SeqData.No_seq - 1; i++) {
        nodes[i].parent = -1;
        nodes[i].ibranch = -1;
        nodes[i].no_child = 0;
        nodes[i].label = 0;
        nodes[i].time = 0;
    }

    while (isspace(ch)) ch = fgetc(treeFile);
    ungetc(ch, treeFile);

    // Flag array to check if the species occurs only once
    for (i = 0; i < SeqData.No_seq; i++)
        check[i] = 0;

    for (;;) {
        ch = fgetc(treeFile);
        if (ch == EOF) return -1;
        else if (!isgraph(ch)) continue;

        else if (ch == '(') {
            level++;
            cnode = tree.no_node++;
            if (cparent >= 0) {
                nodes[cparent].child[nodes[cparent].no_child++] = cnode;
                nodes[cnode].parent = cparent;
                tree.edges[tree.no_edge][0] = cparent;
                tree.edges[tree.no_edge][1] = cnode;
                nodes[cnode].ibranch = tree.no_edge++;
            } else {
                tree.root = cnode;
            }
            cparent = cnode;
        }

        else if (ch == ')') {
            level--;
            inodeb = cparent;
            cparent = nodes[cparent].parent;
        }

        else if (ch == ':') {
            haslength = 1;
            fscanf(treeFile, "%lf", &nodes[inodeb].time);
        }

        else if (ch == '#' || ch == '$') {
            haslabel = 1;
            fscanf(treeFile, "%d", &nodes[inodeb].label);
        }

        else if (ch == ',') {
        }

        else if (ch == ';' && level != 0) {
            error("; in treefile");
        }

        else {
            // Reading the species name or species number
            line[0] = (char) ch;
            line[1] = (char) fgetc(treeFile);

            if (SeqData.No_seq < 10 && isdigit(line[0]) && isdigit(line[1])) {
                ungetc(line[1], treeFile);
                line[1] = 0;
            } else {
                for (i = 1; i < lline;) {
                    if (strchr(delimiters, line[i]) || line[i] == (char) EOF || line[i] == '\n') {
                        ungetc(line[i], treeFile);
                        line[i] = 0;
                        break;
                    }
                    line[++i] = (char) fgetc(treeFile);
                }
            }

            for (j = i - 1; j > 0; j--)
                if (!isgraph(line[j])) line[j] = 0;

            // Determining whether it is a species name or a number
            for (i = 0, hasname = 0; line[i]; i++)
                if (!isdigit(line[i])) {
                    hasname = 1;
                    break;
                }

            if (hasname) {
                for (i = 0; i < SeqData.No_seq; i++)
                    if (!strcmp(line, SeqData.name[i]))
                        break;
                if ((cnode = i) == SeqData.No_seq) {
                    printf("\nError: species %s not found in sequence data.\n", line);
                    exit(-1);
                }
            } else {
                sscanf(line, "%d", &cnode);
                cnode--;
                if (cnode < 0 || cnode >= SeqData.No_seq) {
                    printf("\nError in tree file: species %d not in data.\n", cnode + 1);
                    exit(-1);
                }
            }

            nodes[cnode].parent = cparent;
            nodes[cnode].seq_no = cnode + 1;
            nodes[cparent].child[nodes[cparent].no_child++] = cnode;
            tree.edges[tree.no_edge][0] = cparent;
            tree.edges[tree.no_edge][1] = cnode;
            nodes[cnode].ibranch = tree.no_edge++;
            check[inodeb = cnode]++;

            // Foreground branch judgment (model==3 or 4)
            if (haslabel && nodes[cnode].ibranch >= 0 && nodes[cnode].ibranch < MAXBRANCHES) {
                if (nodes[cnode].label == 1) {
                    is_foreground_branch[nodes[cnode].ibranch] = 1;
                } else {
                    is_foreground_branch[nodes[cnode].ibranch] = 0;
                }
            }
        }

        if (level == 0) break;
    }

    if (popline) fgets(line, 254, treeFile);

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
    int i = 0, j, k;

    // Initial value of phylogenetic tree branch length
    nodes[nodes[tree.root].child[0]].time = 1e-4; // The branches connecting to the root are short
    for (i = 1; i <= tree.no_edge; i++)
        x[i] = nodes[i - 1].time; // Set the length of each branch to the initial value

    // Initial value of kappa transition/transversion ratio
    x[i] = 4.0;

    // Initialize base frequencies (T, C, A)
    x[++i] = 0.25; // T
    x[++i] = 0.25; // C
    x[++i] = 0.25; // A

    // Branching parameters according to model
    if (SeqData.Model == 1 || SeqData.Model == 2) {
        // Site model

        // Initialize omega (two parameters)
        for (j = 0; j < SeqData.omegaCats; j++)
            x[++i] = (j == 0 ? 0.8 : 2.0); // omega_0 = 0.8, omega_1 = 2.0

        // Initialize gamma (one or two parameters)
        for (j = 0; j < SeqData.gammaCats[0]; j++)
            x[++i] = (j == 0 ? 0.8 : 2.0); // omega_0 = 0.8, omega_1 = 2.0

        // Mixture ratio (number of categories - 1)
        int n_cat = SeqData.omegaCats * SeqData.gammaCats[0];
        if (n_cat > 1) {
            for (k = 0; k < n_cat - 1; k++)
                x[++i] = 0.5;
        }
    }
    else if (SeqData.Model == 3) {
        // Branch model alternative hypothesis

        // Initialize omega (one parameter)
        x[++i] = 0.9;

        // Initialize gamma_0 (background), gamma_1 (foreground)
        x[++i] = 0.8;
        x[++i] = 0.7;
    }
    else if (SeqData.Model == 4) {
        // Branch model null hypothesis

        // Initialize omega (one parameter)
        x[++i] = 0.9;

        // Initialize gamma (one parameter)
        x[++i] = 0.8;
    }

    // Record and return the number of parameters
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
    double tb[] = { 1e-4, 0.3};                    // branch length
    double omegabLessThan1[] = { 1e-4, 0.9999 };    // omega_0
    double omegabGreaterThan1[] = { 1.0001, 5.5 };    // omega_1
    double omegab[] = { 1e-4, 5.5 };                  // omega (for branch model)
    double kappab[] = { 1e-4, 5.5 };                  // kappa
    double gammabLessThan1[] = { 1e-4, 1 };         // gamma_0
    double gammabGreaterThan1[] = { 1.0001, 5.5 };    // gamma_1
    double gammab[] = { 1e-4, 5.5 };                  // gamma (for branch model)
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
        // gamma_0 (background), gamma_1 (foreground)
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
    double xPrime;
    double gap = 1e-4;
    double actualLowerbound;

    if (x >= lowerbound)
        xPrime = x;
    else {
        actualLowerbound = lowerbound - gap;

        xPrime = (-gap + gap * sqrt(lowerbound - x)) / (lowerbound - x - 1) + actualLowerbound;
    }
    return xPrime;
}

/**
  * Function that smoothly corrects when the upper limit is exceeded
 **/
double transformXUpper(double x, double upperbound)
{
    double xPrime;
    double gap = 1e-4;
    double actualUpperbound;

    if (x <= upperbound)
        xPrime = x;
    else {
        actualUpperbound = upperbound + gap;

        xPrime = actualUpperbound - (-gap + gap * sqrt(x - upperbound)) / (-upperbound + x - 1);
    }
    return xPrime;
}

/**
  * Objective function combining log-likelihood and penalty terms
 **/
double func(double x[], int n)
{
    int i, j, k;
    double basepi[3];
    double probGamma[3];
    tree.LogL = 0;

    // Branch length setting
    for (i = 0; i < tree.no_edge; i++)
        nodes[i].time = transformXLower(x[i + 1], 1e-4);
        nodes[i].time = transformXUpper(nodes[i].time, 0.3);

    nodes[nodes[tree.root].child[0]].time = 1e-4;

    // kappa setting
    k = tree.no_edge + 1;
    SeqData.kappa = transformXLower(x[k], 1e-4);
    SeqData.kappa = transformXUpper(SeqData.kappa, 5.500);

    // Base frequency setting
    for (j=0; j<3;j++)
		basepi[j] = x[++k];

	SeqData.basepi[0]=exp(basepi[0])/(exp(basepi[0])+exp(basepi[1])+exp(basepi[2])+1);
	SeqData.basepi[1]=exp(basepi[1])/(exp(basepi[0])+exp(basepi[1])+exp(basepi[2])+1);
	SeqData.basepi[2]=exp(basepi[2])/(exp(basepi[0])+exp(basepi[1])+exp(basepi[2])+1);
	SeqData.basepi[3]=1-SeqData.basepi[0]-SeqData.basepi[1]-SeqData.basepi[2];

    // Codon frequency update
    calCodonFreqs();

    // Read parameters for each model

    if (SeqData.Model == 1 || SeqData.Model == 2) {
        SeqData.omega[0] = transformXLower(x[++k], 1e-4);
        SeqData.omega[0] = transformXUpper(SeqData.omega[0], .9999);
        SeqData.omega[1] = transformXLower(x[++k], 1.0001);
        SeqData.omega[1] = transformXUpper(SeqData.omega[1], 5.500);
    }
    else if (SeqData.Model == 3 || SeqData.Model == 4) {
        SeqData.omega[0] = transformXLower(x[++k], 1e-4);
        SeqData.omega[0] = transformXUpper(SeqData.omega[0], 5.500);
    }

    if (SeqData.Model == 1) {
        SeqData.gamma[0][0] = transformXLower(x[++k], 1e-4);
        SeqData.gamma[0][0] = transformXUpper(SeqData.gamma[0][0], .9999);
        SeqData.gamma[0][1] = transformXLower(x[++k], 1.0001);
        SeqData.gamma[0][1] = transformXUpper(SeqData.gamma[0][1], 5.500);

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
        SeqData.gamma[0][0] = transformXUpper(SeqData.gamma[0][0], 5.500);
        SeqData.gamma[0][1] = transformXLower(x[++k], 1e-4);
        SeqData.gamma[0][1] = transformXUpper(SeqData.gamma[0][1], 5.500);
    }
    else if (SeqData.Model == 4) {
        SeqData.gamma[0][0] = transformXLower(x[++k], 1e-4);
        SeqData.gamma[0][0] = transformXUpper(SeqData.gamma[0][0], 5.500);
    }

    // Likelihood calculation (switching functions depending on the model)
    if (SeqData.Model == 3 || SeqData.Model == 4)
        tree.LogL = return_fixGammaP(tree.root);
    else
        tree.LogL = return_gammaPartsP(tree.root);


    // -Return the objective function value including the penalty
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