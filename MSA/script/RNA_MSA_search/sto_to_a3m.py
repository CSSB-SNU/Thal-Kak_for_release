import os
import glob
import argparse
from Bio import AlignIO


def sto_to_a3m(sto_path, output_path):
    """
    Convert Stockholm (.sto) to A3M, redefining Match/Insert states
    relative to the query (first sequence).
    """
    try:
        alignment = AlignIO.read(sto_path, "stockholm")
        num_seqs = len(alignment)
        if num_seqs == 0:
            return False, 0

        # Match columns: alignment columns where the query has a residue
        # (i.e. not a gap). These become the A3M match states.
        query_record = alignment[0]
        query_seq_str = str(query_record.seq)
        match_cols_set = {
            i for i, char in enumerate(query_seq_str) if char not in ("-", ".")
        }

        with open(output_path, "w") as f_out:
            for record in alignment:
                header = f">{record.id}"
                if record.description and record.description != "<unknown description>":
                    header += f" {record.description}"
                f_out.write(f"{header}\n")

                # A3M rules:
                #  - match columns  -> uppercase (gap kept as '-')
                #  - insert columns -> lowercase (gap dropped entirely)
                raw_seq = str(record.seq)
                new_seq = []
                for i, char in enumerate(raw_seq):
                    if i in match_cols_set:
                        new_seq.append("-" if char in ("-", ".") else char.upper())
                    elif char not in ("-", "."):
                        new_seq.append(char.lower())

                f_out.write("".join(new_seq) + "\n")

        print(f"[Converted] {sto_path} -> {output_path} ({num_seqs} seqs)")
        return True, num_seqs

    except Exception as e:
        print(f"Error converting {sto_path}: {e}")
        return False, 0


def convert(input_path, output_path):
    """Convert .sto to .a3m. Input can be a single file."""
    if os.path.isfile(input_path):
        sto_to_a3m(input_path, output_path)
    else:
        print(f"Input {input_path} is not a file. Please provide a valid .sto file.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help=".sto file or directory")
    parser.add_argument("--output", required=True, help=".a3m file or directory")
    args = parser.parse_args()
    convert(args.input, args.output)
