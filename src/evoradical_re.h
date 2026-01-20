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
 *  Limitations: 
 *  - Only one tree structure can be input at a time
 *  - Omega constant among sites
 *  - Initial values of kappa and omega fixed
 *
 *  -- Substitution model --
 *  q_{i,j} = pi_j                                      (synonymous, transversion) / 
 *            kappa * pi_j                              (synonymous, transition) / 
 *            gamma_{i,j} * omega_{i, j} * pi_j         (nonsynonymous, transversion) /
 *            gamma_{i,j} * omega_{i, j} * kappa * pi_j (nonsynonymous, transition)
 *
 *              where, gamma_{i,j} = 1 and omega_{i, j} = omega, 
 *                     if the change in physicochemical properties between i and j is conservative, /
 *                     gamma_{i,j} = gamma and omega_{i, j} = 1, 
 *                     if the change in physicochemical properties between i and j is radical.
 *
 *  -- Site model --
 *      -alternative hypothesis (Model==1)-
 *       site categories
 *       (1) omega_0 <= 1, gamma_0 <= 1  (both purifying or neutral)
 *       (2) omega_0 <= 1, gamma_1 > 1   (omega: purifying or neutral, gamma: positive)
 *       (3) omega_1 > 1, gamma_0 <= 1   (omega: positive, gamma: purifying or neutral)
 *       (4) omega_1 > 1, gamma_1 > 1    (both positive)
 *
 *      -null hypothesis (Model==2)-
 *       (1) omega_0 <= 1, gamma_0 <= 1  (both purifying or neutral)
 *       (3) omega_1 > 1, gamma_0 <= 1   (omega: positive, gamma: purifying or neutral)
 *
 *  --branch model--
 *      -alternative hypothesis (Model==3)-
 *       Same omega value for all branches: omega_0
 *       Background branch: gamma_0
 *       Foreground branch: gamma_1
 *
 *      -null hypothesis (Model==4)-
 *       Same omega and gamma values for all branches: omega_0, gamma_0
 *************************************/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "toolsfromPAML.h"
#include "xdfpmin.c"

#define MAXLINE	      256
#define NAMELENGTH    50    //length of sequence names
#define NS            100   //maximum # of sequences
#define NCODON        64    //maximum # of codon
#define SEQLENGTH     10000 //maximum sequence length
#define MAXPARTS      5     //maximum number of partitions
#define MAXCATS	      8     //maximum number of categories
#define NNODE         (NS*2-1)
#define NEDGE         (NS*2-2)
#define NP	          1000  //max # of parameters that needed to be optimized
#define MAXSTRING     1000
#define PENALTYUPPER  1
#define PENALTYLOWER  1
#define GTOL          1.0e-4
#define MAXBRANCHES   NEDGE   //maximum number of foreground branch

struct node {
    char name[NAMELENGTH];
    int code[NNODE];
    double likelihood[SEQLENGTH * NCODON];	/*prob. of all the leaves below the node */
    double time;
    int node_no;
    int seq_no;
    int parent;
    int child[NS];
    int no_child;
    int label;
    int ibranch;
};
typedef struct node NODE;
typedef NODE *SEQUENCE;
SEQUENCE nodes;

// tracking the time
typedef struct {
    clock_t begin_clock, save_clock;
    time_t begin_time, save_time;
} time_keeper;

struct data {
    char *name[NS];
    char *seq[NS];
    int *Codonseq[NS];
    int Seq_length, No_seq, np; //sequence length, number of sequences, number of parameters
    double pi[NCODON];		//stationary distribution of the codons
    double normCodon_pi[NCODON];
    double basepi[4];   //stationary distribution of the nucleotides
    int ncode;
    int No_codon, CPatternCounts[SEQLENGTH/3], CMapPattern[SEQLENGTH/3],CNo_pattern;
    int SeqType; //whether it is only noncoding, only coding, or mixed
    int Model; 	//Model = 1, gamma*omega, both have two categories (alternative hypothesis in site model)
                //(puryfying or neutral (omega or gamma<=1)/ positive (omega or gamma>1))
                //Model = 1: alternative hypothesis in site model
                //           omega_0 <= 1, omega_1 > 1
                //           gamma_0 <= 1, gamma_1 > 1
                //           4 site categories: (omega_0,gamma_0), (omega_0,gamma_1), (omega_1,gamma_0), (omega_1,gamma_1)


                //Model = 2: null hypothesis in site model
                //           omega_0 <= 1, omega_1 > 1
                //           gamma_0 <= 1
                //           2 site categories: (omega_0,gamma_0), (omega_1,gamma_0)


                //Model = 3: alternative hypothesis in branch model
                //           Same omega value for all branches
                //           gamma_bg, gamma_fg


                //Model = 4: null hypothesis in branch model
                //           Same omega and gamma values for all branches

    double omega[MAXCATS], kappa, gamma[MAXPARTS][MAXCATS];
    int omegaCats, gammaCats[MAXPARTS]; //number of categories
    double probGamma[4]; //probability of each category of omega/gamma
    double postdistGamma[SEQLENGTH*4];
    char seqFName[MAXLINE], codingFName[MAXLINE], treeFName[MAXLINE], outFName[MAXLINE]; //file names
    int no_parts; //number of amino acid partitions
    char *partFNames[MAXPARTS]; //file names of the amino acid partitions
    int AAPartition[MAXPARTS][19*20+1]; //holds the indicator variable about whether a AA substition is between partition
} SeqData;


struct treestrct {
    int no_node;
    double LogL;    //loglikelihood of the tree
    int no_edge;
    int root;
    int edges[NEDGE][2];
} tree;

static time_keeper tk;

extern int GeneticCode[][64];
extern int is_foreground_branch[MAXBRANCHES];  // foreground branch flag（1/0）
extern int *site_class;  // [No_codon] each site category（0〜3）
extern double gamma0, gamma1;  // gamma (background, foreground)
extern double omega0, omega1;  // omega (background, foreground)

void start_time(void);
double prn_time(void);

SEQUENCE new_node(void);    //creating a new node
SEQUENCE init_node(int node_no, int seq_no);    //creating a tree
int new_tree(int No_nodes); //allocating memory for a new tree;

int InitializeLkl(double omega, int AAPartition[],double gamma);
int return_p(int node_i);
double return_LH(int node_i); //returns the loglikelihood for one rate models
double return_gammaPartsP(int node_i);
double return_fixGammaP(int node_i);  // likelihood caliculation for  model==3,4

int calc_post(int category); //calculate the posterior distribution
void freeall(void);
int Qmatrix(double Root[], double U[], double V[],
	    double kappa, double omega, double gamma,int AAPartition[], double Q[], int n, int branch_i);


FILE *gfopen(char *filename, char *mode);
int ReadOptionFile (char* optionFName);
int Readfile(FILE * fp);
int ReadOpt(FILE * fp);	//reads in the coding/noncoding regions
int printX(double x[], FILE * fp);
int FindPattern();
int MapPattern();
int GetAAPart (int gamma_i, char* ZpartAAFName); //gets the AA partition from a file
int String2Ints (char* line, int outputInts[], int no_parts, const char *szSeparator);
int String2Strs (char* line, char* outputStrings[], int no_parts, const char *szSeparator);
int String2Doubles (char* line, double outputDoubles[], int no_parts, const char *szSeparator);

int AA2Code (char b);
int Char2Code(char c);
int Codon2Code61(char codon[3]);
int Codon2Code64(char codon[3]);
int calCodonFreqs(void); //calculate the codon frequencies from the base frequencies
void missingData(int baseorcodon);

int setFROM_61_64(void);

/* files from tools.c */
int matout1(double x[], int n, int m);

int from64_61[64] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -1, 10, 11, -1, 12,
    13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
    30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46,
    47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60};
char AACodes[] = "ARNDCQEGHILKMFPSTWYV*-?";

int FROM61[64], FROM64[64];
double *PMat, *U, *V, *Root;
double *BasePMat, *BaseU, *BaseV, *BaseRoot;
int iter; //number of iterations

int ReadTree(FILE * treeFile, int *length_label, int popline);

//optimization stuff
int InitializeX(double x[]);
int SetxBound(int np, double xb[][2]);
double quadpenalty(double x[], int n);
//double transformX(double x, double lowerbound);
double transformXLower(double x, double lowerbound);
double transformXUpper (double x, double upperbound);
double func(double x[], int n);