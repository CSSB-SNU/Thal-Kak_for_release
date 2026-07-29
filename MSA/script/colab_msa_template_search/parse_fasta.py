import sys, os, glob, re
import argparse

"""make query fasta with casp fasta file considering proper stoichiometry"""

NUCLEIC_ACIDS = set("ACGTUN")


def is_nucleic_acid(seq):
    return set(seq.upper()).issubset(NUCLEIC_ACIDS)


def get_na_type(seq):
    """Return 'dna' or 'rna' based on sequence content."""
    if "U" in seq.upper():
        return "rna"
    return "dna"


def parse_fasta(fasta, stoi, out_dir):

    with open(fasta, "r") as f:
        data = f.read()

    k = os.path.basename(fasta).split(".")[0]
    chains = data.split("\n")[1::2]
    stoi_parsed = re.findall(r"([A-Z])(\d+|n)", stoi)
    protein_entities = []  # list of (seq, n)
    na_chains = []
    for seq, (chain_chr, n) in zip(chains, stoi_parsed):
        n = 1 if n == "n" else int(n)
        if is_nucleic_acid(seq):
            na_type = get_na_type(seq)
            na_chains.append({"sequence": seq, "type": na_type, "copy": n})
            print(f"  Found {na_type.upper()} chain {chain_chr} (copy: {n})")
            continue
        protein_entities.append((seq, n))

    # Interleave copies: round-robin so entity 1 gets chains A, D, ... entity 2 gets B, E, ...
    final = []
    max_copies = max((n for _, n in protein_entities), default=0)
    for r in range(max_copies):
        for seq, n in protein_entities:
            if r < n:
                final.append(seq)
    parsed_path = f"{out_dir}/{k}_parsed.fa"
    with open(parsed_path, "w") as fasta_file:
        fasta_file.write(f">{k}\n")
        fasta_file.write(":".join(final))

    return parsed_path, na_chains

def main(args):
    fasta = args.fasta
    stoi = args.stoi
    out_dir = args.out_dir
    
    out_fa = parse_fasta(fasta, stoi, out_dir)
    print(out_fa)

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    argparser.add_argument("--fasta", default=None, help="Path to input fasta")
    argparser.add_argument("--stoi", default="A1", help="stoichiometry information")
    argparser.add_argument("--out_dir", default=None, help="Path to output directory")
    args = argparser.parse_args()    
    main(args)
