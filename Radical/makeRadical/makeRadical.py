import argparse
import logging
import sys
import json
from datetime import datetime

EXPECTED_AAS = set("ARNDCQEGHILKMFPSTWYV")


def setup_logger(logfile='makeRadical.log'):
    logging.basicConfig(
        filename=logfile,
        filemode='w',
        level=logging.INFO,
        format='%(message)s'
    )
    start_time = datetime.now()
    logging.info(f"# Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"# Command: {' '.join(sys.argv)}\n")


def finalize_logger():
    end_time = datetime.now()
    logging.info(f"\n# End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")


def read_aa_index(file_path):
    aa_dict = {}
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) != 2:
                    logging.warning(f"Ignored invalid line -> {line}")
                    continue
                aa, value = parts
                if aa in aa_dict:
                    logging.error(f"Duplicate entry for amino acid: {aa}")
                    raise SystemExit(1)
                try:
                    aa_dict[aa] = float(value)
                except ValueError:
                    logging.error(f"Value conversion failed -> {line}")
                    raise SystemExit(1)

        logging.info(f"{len(aa_dict)} amino acid values loaded.")
        for aa, val in aa_dict.items():
            logging.info(f"{aa}: {val}")

        input_set = set(aa_dict.keys())
        missing = EXPECTED_AAS - input_set
        extra = input_set - EXPECTED_AAS

        if missing:
            logging.error(f"Missing amino acids: {' '.join(sorted(missing))}")
        if extra:
            logging.error(f"Unexpected amino acids: {' '.join(sorted(extra))}")

        if missing or extra:
            logging.error("Input does not contain exactly 20 standard amino acids. Exiting.")
            raise SystemExit(1)

    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
        raise SystemExit(1)
    return aa_dict


def load_codon_table(json_path, codon_table_name):
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)

        if codon_table_name not in data["genetic_codes"]:
            logging.error(f"Codon table '{codon_table_name}' not found in {json_path}")
            raise SystemExit(1)

        codon_table = data["genetic_codes"][codon_table_name]
        codon_order = data["codon_order"]
        aa_indices = data["amino_acid_indices"]

        if len(codon_table) != 64 or len(codon_order) != 64:
            logging.error("Codon table or codon order does not contain 64 entries.")
            raise SystemExit(1)

        logging.info(f"Codon table loaded: {codon_table_name} ({json_path})")
        return codon_order, codon_table, aa_indices

    except FileNotFoundError:
        logging.error(f"Genetic code file not found: {json_path}")
        raise SystemExit(1)
    except json.JSONDecodeError:
        logging.error("Failed to parse JSON file.")
        raise SystemExit(1)


def generate_nonsynonymous_pairs(codon_order, codon_table):
    bases = ['T', 'C', 'A', 'G']
    pairs = set()

    for idx, codon in enumerate(codon_order):
        aa1_index = codon_table[idx]
        if aa1_index == -1:
            continue

        for i in range(3):
            for b in bases:
                if b == codon[i]:
                    continue
                mutated = codon[:i] + b + codon[i+1:]
                if mutated in codon_order:
                    m_idx = codon_order.index(mutated)
                    aa2_index = codon_table[m_idx]
                    if aa2_index == -1 or aa2_index == aa1_index:
                        continue
                    aa1 = get_aa_symbol_by_index(aa1_index)
                    aa2 = get_aa_symbol_by_index(aa2_index)
                    pairs.add(tuple(sorted((aa1, aa2))))

    return sorted(pairs)


def get_aa_symbol_by_index(index):
    aa_order = "ARNDCQEGHILKMFPSTWYV"
    return aa_order[index]


def main():
    parser = argparse.ArgumentParser(description="Generate radical AA substitution list.")
    parser.add_argument("-i", "--input", required=True, help="AAindex file (e.g., AAindex.dat)")
    parser.add_argument("-c", "--codon_table", required=True, help="Codon table name (e.g., universal)")
    parser.add_argument("-j", "--json", default="genetic_code.json", help="Genetic code JSON file")
    parser.add_argument("-s", "--split", type=float, default=0.3, help="Upper percentile for radical substitutions")
    parser.add_argument("-o", "--output", required=True, help="Output radical pair file")
    args = parser.parse_args()

    setup_logger()
    logging.info("Started reading AAindex")
    aa_values = read_aa_index(args.input)
    logging.info("Finished reading AAindex\n")

    logging.info("Started loading codon table")
    codon_order, codon_table, aa_indices = load_codon_table(args.json, args.codon_table)
    logging.info("Finished loading codon table\n")

    logging.info("Generating nonsynonymous amino acid substitution pairs")
    aa_pairs = generate_nonsynonymous_pairs(codon_order, codon_table)
    logging.info(f"Identified {len(aa_pairs)} unique nonsynonymous amino acid pairs.")

    changes = [abs(aa_values[a] - aa_values[b]) for a, b in aa_pairs]
    logging.info("AA_pair={" + ", ".join(f"{a}{b}" for a, b in aa_pairs) + "}")
    logging.info("AAindex_change={" + ", ".join(f"{v:.3f}" for v in changes) + "}")

    # sort and determine threshold
    zipped = sorted(zip(changes, aa_pairs), reverse=True)
    threshold_index = int(len(zipped) * args.split)
    radical_pairs = [pair for _, pair in zipped[:threshold_index]]

    with open(args.output, 'w') as f:
        for a1, a2 in radical_pairs:
            f.write(f"{a1}{a2} ")

    logging.info(f"Radical substitutions written to {args.output}")
    finalize_logger()


if __name__ == "__main__":
    main()
