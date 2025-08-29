####################################
Project Overview
####################################

- This project is designed to provide a tool for identifying and analyzing selective pressures on protein-coding nucleotide sequences, specifically focusing on physicochemical factors. 

- The primary objective of this tool is to evaluate the selective pressure acting on protein-coding sequences using evolutionary models.
- By identifying patterns in amino acid substitution and correlating them with physicochemical properties, the tool helps researchers better understand how these factors contribute to evolutionary processes, particularly under different environmental pressures.

- This tool is based on the original EvoRadical program, which was described in the following paper:

Wong, W.S., Sainudiin, R., & Nielsen, R. (2006). Identification of physicochemical selective pressure on protein encoding nucleotide sequences. BMC Bioinformatics 7, 148. 
https://doi.org/10.1186/1471-2105-7-148

- The original EvoRadical implementation only supported site models, but this version has been expanded to include branch models, offering more flexibility and depth in analyzing evolutionary pressures across phylogenetic trees. 
- Additionally, some of the site models have been compacted and simplified to streamline the analysis while maintaining the ability to capture key evolutionary patterns.

< Models and Algorithms >
- The tool employs a combination of site models and branch models to analyze evolutionary changes in protein-coding sequences
- The tool uses maximum likelihood estimation (MLE) to infer the most probable evolutionary tree and parameter values.

-- Codon substitution model --

q_{i,j} = pi_j                                      (synonymous, transversion) / 
          kappa * pi_j                              (synonymous, transition) / 
          gamma_{i,j} * omega_{i, j} * pi_j         (nonsynonymous, transversion) /
          gamma_{i,j} * omega_{i, j} * kappa * pi_j (nonsynonymous, transition)

            where, gamma_{i,j} = 1 and omega_{i, j} = omega, 
                   if the change in physicochemical properties between i and j is conservative, /
                   gamma_{i,j} = gamma and omega_{i, j} = 1, 
                   if the change in physicochemical properties between i and j is radical.

-- Site Models --

- These models are used to assess evolutionary pressures at individual sites along the sequence.
- They enable the identification of selective forces operating at specific positions within the protein.

    - alternative hypothesis (Model==1) -
       site categories
       (1) omega_0 <= 1, gamma_0 <= 1  (both purifying or neutral)
       (2) omega_0 <= 1, gamma_1 > 1   (omega: purifying or neutral, gamma: positive)
       (3) omega_1 > 1, gamma_0 <= 1   (omega: positive, gamma: purifying or neutral)
       (4) omega_1 > 1, gamma_1 > 1    (both positive)

      - null hypothesis (Model==2) -
       (1) omega_0 <= 1, gamma_0 <= 1  (both purifying or neutral)
       (3) omega_1 > 1, gamma_0 <= 1   (omega: positive, gamma: purifying or neutral)

-- Branch Models --

- Originally, EvoRadical used only site models.
- However, in this project, the code has been modified to incorporate branch models.
- Branch models enable the analysis of selective pressures along different branches of a phylogenetic tree, allowing for the differentiation of pressures acting on different lineages.

    - alternative hypothesis (Model==3) -
       Same omega value for all branches: omega_0
       Background branch: gamma_0
       Foreground branch: gamma_1

    - null hypothesis (Model==4) -
       Same omega and gamma values for all branches: omega_0, gamma_0


####################################
Directory Structure
####################################

- Before typing make, ensure that the directory structure is set up as follows:

Radical/
│
├── ctl/                    # Control files directory
│   └── evoRadical.dat      # Control file for input and model settings
│
├── input/                  # Input files directory
│   ├── compAAParts.dat     # Amino acid pair information for radical substitutions
│   ├── RHO_seq.dat         # Input nucleotide sequence data
│   └── RHO.fa.tree         # Input phylogenetic tree file
│
├── makeRadical/            # Directory of Python scripts for generating radical amino acid pair data
│   └── README.md etc.      # (Refer to makeRadical/README.md for more details)
│
├── src/                    # Source code directory
│   ├── eigen.c             # Eigenvalue and eigenvector calculation functions
│   ├── evoradical_re.c     # Main EvoRadical implementation (modified)
│   ├── evoradical_re.h     # Header file for EvoRadical implementation
│   ├── toolsfromPAML.c     # Tools derived from PAML
│   ├── toolsfromPAML.h     # Header file for PAML tools
│   └── xdfpmin.c           # DFP minimization routine (optimization)
│
├── Makefile                # Makefile to build the project
├── README.md               # This file (README)
└── bin/                    # Directory where compiled executables will be placed
└── obj/                    # Directory where object files will be placed


####################################
Building the Project
####################################

- To compile the project, simply run the make command from the root directory:
- This will compile the source code, generating the necessary object files in the obj directory 
  and the final executable (reEvoRadical) in the bin directory.


####################################
Input Files
####################################

- There are three main input files required to run the program:

- Sequence File (e.g., RHO_seq.dat)
    - The sequence data file must follow a modified FASTA format.
    - The first line should contain two integers: the number of sequences and the length of each sequence.
    - The sequences should not contain any stop codons.
    - Each sequence is preceded by a header line, which starts with a > symbol.

    Example:

    12 1044
    >Hamp
    ATGAACGGGACGGAGGGCCCGCCCCCTCGTCGGCTGGTCCAGGTACATCCCGG
    ...
    >Ppho
    ATGAACGGGACGGAGGGCCTGAACTTCTACGTGCCTTTCTCTAACAAGACAGGC
    ...
    >Dleu
    ATGAACGGGACAGAGGGCCTGAACTTCTACGTGCCTTTCTCTAACAATACAGGCG
    ...

- Phylogenetic Tree File (e.g., seq.tree)
    - The phylogenetic tree file must be in Nexus format.
    - The tree should include branch lengths, which will be used as initial values for the analysis.
    - The system supports trees generated by tools like IQ-TREE with codon substitution models.
    - When using branch models, mark the foreground branches using like "#1" in the tree file.

    Example tree snippet:

    ((A,B)C:0.1,(D,E)F:0.2)G:0.3;

- Amino Acid Pair Information File (e.g., compAAParts.dat)
    - This file contains information about amino acid substitutions that involve 
      radical physicochemical property changes. The pairs should be listed in the following format:

    Example:

    YG YA YP YV YI YL YM YF WG WA WP ...
    
    Each amino acid pair corresponds to a substitution that is considered radical
    (i.e., resulting in a significant physicochemical change in the protein structure).


####################################
Control File (e.g., evoRadical.dat)
####################################

- The control file contains the necessary configurations and parameters to run the analysis.
- Below is the structure and description of the key parameters in the evoRadical.dat file.

    Example of a ctl file:

    Model = 1
    seqfile = input/RHO_seq.dat
    GammaFiles = input/compAAParts.dat
    treefile = input/RHO.fa.tree
    outfile = Output/Site/H1/mlradSH1_2.out

    Parameters:
        Model: This specifies the model to use for the analysis:
                1: Site model (alternative hypothesis)
                2: Site model (null hypothesis)
                3: Branch model (alternative hypothesis)
                4: Branch model (null hypothesis)

                The choice of model depends on whether you want to analyze site-specific 
                or branch-specific evolutionary pressures.
        
        seqfile: Path to the input nucleotide sequence file (e.g., RHO_seq.dat).
                This file must follow the format described earlier.
        
        GammaFiles: Path to the file containing the radical amino acid substitution pairs 
                (e.g., compAAParts.dat).
                This file should list pairs of amino acids involved in radical physicochemical changes.
        
        treefile: Path to the input phylogenetic tree file (e.g., RHO.fa.tree).
                This tree file should be in Nexus format, and it must include branch lengths.
        
        outfile: Path to the output file where the results of the analysis will be written 
                (e.g., mlradSH1_2.out).
                This file will contain the log-likelihood, parameter estimates, 
                and other results from the analysis.


####################################
Running the Program
####################################

- To run the program, use the following command in your terminal:

    reEvoRadical ctl/evoRadical.dat


####################################
Output Description
####################################

- The output file generated by reEvoRadical contains detailed information about the evolutionary analysis.
- It includes model-specific parameters, branch lengths, category probabilities, 
  and posterior probabilities for each site.
- The content of the output is divided into several sections as follows:

    - Model Summary
        - Model used: Indicates the model type that was applied.
        - Sequence file: The file containing the input sequences used for the analysis.
        - Tree file: The file containing the input phylogenetic tree used for the analysis.
        - Log likelihood: The log likelihood value for the tree under the given model and parameters.
        - np: The number of free parameters used in the model.

    - Start time: The date and time when the analysis started.

    - Base Frequencies: The base frequencies for each of the four nucleotide bases (T, C, A, G) are given.
                        These frequencies represent the relative abundance of each base across the sequences.
    
    - Branch Lengths: The branch lengths of the phylogenetic tree are listed for each branch.
                      Each branch length corresponds to the amount of evolutionary change along that branch.
    
    - Parameters
        - Kappa: The transition/transversion rate ratio.
        - Omega[0]: The omega value for the first category
        - Omega[1]: The omega value for the second category.
        - Gamma[0][0]: The gamma value for the first category.
        - Gamma[0][1]: The gamma value for the second category.
    
    - Category Probabilities (only site model)
      The category probabilities represent the proportion of time each category (omega and gamma) 
      is used for the evolutionary model at each site.
      These probabilities are calculated based on the model and the data and reflect the likelihood of each category occurring in the sequences.
    
    - Posterior Probabilities Per Site (only site model)
      For each site in the sequence data, posterior probabilities for each category are calculated.
      These values indicate the likelihood of each category being the most likely for that specific site.
      The output shows the probabilities for each of the 4 or 2 categories at each site in the alignment.

    - End time: The time when the analysis finished.

####################################
Credits
####################################

This software is based on the original EvoRadical code (Wong, Sainudiin, & Nielsen 2006), 
with modifications and improvements by Hayate Takeuchi. 
The original implementation was developed by Wendy Wong in 2004.

Eigenvalue computation for real nonsymmetric matrices
The source code for this part was written by Tianlin Wang at the University of Illinois.
The underlying algorithm for eigenvalue computation follows the work in Handbook for 
Automatic Computation, vol 2 by Wilkinson and Reinsch (1971).
Most of the source code is based on a public domain software package called MATCALC.

Local optimization
The code for local optimization was adopted from Numerical Recipes in C (Press et al., 1992).
(C) Copr. 1986-92, Numerical Recipes Software.

Program organization
The structure of the code was heavily inspired by the codeml.c source code from the PAML package
(Yang, 1997-2002).

Branch model implementation
A branch model has been added to this version, along with the removal of some original models,
by Hayate Takeuchi.
