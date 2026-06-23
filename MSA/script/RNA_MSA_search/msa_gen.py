import os
import subprocess
import argparse
import sys
from Bio import SeqIO


def run_command(cmd):
    """Run shell command silently unless error occurs"""
    try:
        subprocess.run(
            cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        # Print only on failure.
        print(
            f"Command failed: {' '.join(cmd)}\nError: {e.stderr.decode()}",
            file=sys.stderr,
        )
        raise e


def get_query_length(query_path):
    record = next(SeqIO.parse(query_path, "fasta"))
    return len(record.seq)


def run_nhmmer(query_path, db_prefix, output_tbl, query_len, num_cpu=4):
    # AF3 Params
    f3_val = "0.02" if query_len < 50 else "0.00005"
    cmd = [
        "nhmmer",
        "-E",
        "0.001",
        "--incE",
        "0.001",
        "--rna",
        "--watson",
        "--F3",
        f3_val,
        "--tblout",
        output_tbl,
        "-o",
        "/dev/null",
        "--cpu",
        str(num_cpu),
        query_path,
        db_prefix,
    ]
    run_command(cmd)


def extract_hits(
    tbl_file, fasta_db, output_hits_fasta, temp_prefix, max_sequences=10000
):
    hits = []
    if not os.path.exists(tbl_file):
        return False

    with open(tbl_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 1:
                target_name = parts[0]
                if target_name not in hits:
                    hits.append(target_name)

    if not hits:
        return False

    # Cutoff applied
    if len(hits) > max_sequences:
        # print(f"    [Info] Cutoff applied: {len(hits)} -> {max_sequences} hits")
        hits = hits[:max_sequences]

    # Unique temp list file name
    temp_list = f"{temp_prefix}.list"

    with open(temp_list, "w") as f:
        for h in hits:
            f.write(h + "\n")

    # esl-sfetch
    cmd = ["esl-sfetch", "-f", "-o", output_hits_fasta, fasta_db, temp_list]
    run_command(cmd)

    if os.path.exists(temp_list):
        os.remove(temp_list)
    return True


def remove_duplicates(fasta_file):
    if not os.path.exists(fasta_file):
        return

    seen_sequences = set()
    unique_records = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_str = str(record.seq).upper()
        if seq_str not in seen_sequences:
            seen_sequences.add(seq_str)
            unique_records.append(record)

    # Overwrite if duplicates found
    SeqIO.write(unique_records, fasta_file, "fasta")


def realign_hits(query_path, hits_fasta, output_sto, temp_hmm):
    # 1. Build Query Profile (Unique name)
    run_command(["hmmbuild", "--rna", temp_hmm, query_path])

    # 2. Align Hits
    run_command(
        [
            "hmmalign",
            "--rna",
            "--mapali",
            query_path,
            "-o",
            output_sto,
            temp_hmm,
            hits_fasta,
        ]
    )

    if os.path.exists(temp_hmm):
        os.remove(temp_hmm)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--query", required=True)
    parser.add_argument("--db_dir", required=True)
    parser.add_argument("--output", default="final_msa.sto")
    args = parser.parse_args()

    # Per-query unique id so parallel runs do not clash on temp files.
    # Derive the id from the query filename (e.g. R1203.fasta -> R1203).
    query_id = os.path.splitext(os.path.basename(args.query))[0]

    # Prefix temp files with query_id to keep them unique.
    temp_combined_hits = f"temp_{query_id}_combined.fasta"
    temp_hmm_file = f"temp_{query_id}_query.hmm"

    dbs = [
        {
            "name": "Rfam",
            "hmm": "rfam/rfam_v_latest.mdf",
            "seq": "rfam/rfam_clust_rep_seq.fasta",
        },
        {
            "name": "RNAcentral",
            "hmm": "rnacentral/rnacentral_v_latest.mdf",
            "seq": "rnacentral/rnacentral_clust_rep_seq.fasta",
        },
    ]

    try:
        query_len = get_query_length(args.query)

        # Init combined file
        open(temp_combined_hits, "w").close()
        found_any = False

        for db in dbs:
            # print(f"--> Searching {db['name']}...")
            hmm_path = os.path.join(args.db_dir, db["hmm"])
            seq_path = os.path.join(args.db_dir, db["seq"])

            # Unique temp filenames.
            tbl_out = f"temp_{query_id}_{db['name']}.tbl"
            hits_out = f"temp_{query_id}_{db['name']}.fasta"
            # Prefix for the esl-sfetch list.
            list_prefix = f"temp_{query_id}_{db['name']}_list"

            try:
                # 1 CPU per search in parallel mode to save resources
                run_nhmmer(args.query, hmm_path, tbl_out, query_len, num_cpu=1)

                if extract_hits(tbl_out, seq_path, hits_out, list_prefix):
                    with (
                        open(hits_out, "r") as infile,
                        open(temp_combined_hits, "a") as outfile,
                    ):
                        outfile.write(infile.read())
                    found_any = True
            except Exception as e:
                print(f"[Warning] Error in {db['name']} for {query_id}: {e}")

            # Cleanup DB specific temps
            if os.path.exists(tbl_out):
                os.remove(tbl_out)
            if os.path.exists(hits_out):
                os.remove(hits_out)

        if found_any:
            # print("--> Removing duplicates...")
            remove_duplicates(temp_combined_hits)
            # print("--> Generating Final MSA...")
            realign_hits(args.query, temp_combined_hits, args.output, temp_hmm_file)
        else:
            print(f"No hits found for {query_id}")

    finally:
        # Final cleanup of temp files.
        if os.path.exists(temp_combined_hits):
            os.remove(temp_combined_hits)
        if os.path.exists(temp_hmm_file):
            os.remove(temp_hmm_file)


if __name__ == "__main__":
    main()
