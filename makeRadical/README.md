# makeRadical.py

## Overview

`makeRadical.py` is a command-line tool to generate amino acid substitution partitions used in the **EvoRadical (Wong, Sainudiin, & Nielsen 2006)** evolutionary framework.

It identifies "radical" versus "conservative" amino acid changes based on a given AAindex (amino acid property values), restricted only to **nonsynonymous single nucleotide mutations** according to a specified codon table.

This allows for biologically meaningful classification of radical substitutions based on molecular evolution principles.

The output format is compatible with **EvoRadical**, which expects a flat list of radical amino acid substitution pairs.

## Features

- Loads and validates a 20-amino-acid AAindex file.
- Loads codon table definitions from an external JSON file.
- Constructs all single-nucleotide nonsynonymous substitutions based on the selected codon table.
- Calculates property differences between amino acids (absolute difference).
- Splits substitutions into "conservative" and "radical" using quantile-based threshold.
- Writes output in EvoRadical-compatible format.
- Logs all operations and errors in detail.

---

## Usage

```bash
python3 makeRadical.py -i compressibility.dat -c universal -s 0.4 -o compAAParts.dat
```

## Arguments

-i, --input         :Input AAindex file (format: A 10.0, one per line for 20 amino acids)
-c, --codon_table   :Codon table name defined in genetic_code.json (e.g., universal, yeast_mt, etc.)
-s, --split         :Split threshold (e.g., 0.4 for top 50% to be labeled as radical substitution)
-o, --output        :Output file to save radical amino acid pairs
-j, --json          :Codon table JSON file (default: genetic_code.json)

## Input File Formats

1. AAindex File (e.g. compressibility.dat)

    - Should contain exactly 20 unique amino acid single-letter codes with numeric values.

    - Example:

        A 12.5
        R 8.3
        N 9.1
        ...

2. Codon Table JSON (default: genetic_code.json)
    
    - Must include:
        - amino_acid_indices :mapping from single-letter amino acids to integers
        - codon_order        :array of 64 codons in standard order
        - genetic_codes      :dictionary mapping codon table names to lists of 64 values

## Output

A space-separated list of radical substitution amino acid pairs, such as:

    AE PQ EG HP LP PT AD PS CR

This is the format expected by EvoRadical for identifying radical amino acid transitions.

## Logging

- All processing steps and errors are logged to makeRadical.log. 

- The log includes:

    - Start and end timestamps
    - Input command
    - AAindex parsing result
    - Codon table loading status
    - Identified nonsynonymous substitutions
    - Amino acid pair differences

## Author

- Hayate Takeuchi
- Created: March 31, 2025
- License: MIT License

This code is open-source and free to use, modify, and distribute under the terms of the MIT License.