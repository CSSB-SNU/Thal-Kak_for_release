import argparse, subprocess, yaml, os, glob
from MSA.script.colab_msa_template_search.parse_fasta import parse_fasta
from MSA.script.colab_msa_template_search.colab_a3m_to_yaml import (
    split_colab_a3m_write_yaml,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
with open(os.path.join(ROOT, "config.yaml")) as _f:
    CONFIG = yaml.safe_load(_f)


def msa_generation(args):
    base_dir = args.output_dir or os.path.dirname(os.path.abspath(args.seq))
    target_name = os.path.basename(args.seq).split(".")[0]

    print("Generating MSA...")
    parsed_fasta_path, na_chains = parse_fasta(args.seq, args.stoi, base_dir)

    msa_dir = os.path.join(base_dir, "msa", args.msa)
    os.makedirs(msa_dir, exist_ok=True)

    # Prior run's MSA params (read before method_log.yaml is rewritten below) so
    # the RNA MSA search can be skipped when seq/stoi/method are unchanged — e.g.
    # running several structure models on the same target.
    _prior_path = os.path.join(msa_dir, "method_log.yaml")
    _prior = {}
    if os.path.exists(_prior_path):
        with open(_prior_path) as f:
            _prior = yaml.safe_load(f) or {}
    params_match = (
        _prior.get("msa") == args.msa
        and _prior.get("seq") == args.seq
        and _prior.get("stoi") == args.stoi
    )

    # Check if there are protein chains to search
    with open(parsed_fasta_path, "r") as f:
        protein_seq = f.read().split("\n", 1)[1].strip()
    has_protein = len(protein_seq) > 0

    if has_protein:
        method_log_path = os.path.join(msa_dir, "method_log.yaml")
        existing_a3m = glob.glob(os.path.join(msa_dir, "*.a3m"))
        skip_msa = False
        if os.path.exists(method_log_path) and existing_a3m:
            with open(method_log_path) as f:
                method_log = yaml.safe_load(f)
            if (
                method_log.get("msa") == args.msa
                and method_log.get("seq") == args.seq
                and method_log.get("stoi") == args.stoi
            ):
                skip_msa = True

        if skip_msa:
            print(
                "MSA already generated with the same parameters, skipping MSA generation."
            )
            output_msa = existing_a3m[0]
            output_yaml = split_colab_a3m_write_yaml(output_msa)
        else:
            match args.msa:
                case "colab":
                    print("Running MSA generation using Colab...")
                    subprocess.run(
                        f"colabfold_batch --msa-only --templates "
                        f"{parsed_fasta_path} {msa_dir}",
                        shell=True,
                        check=True,
                    )
                    output_msa = glob.glob(os.path.join(msa_dir, "*.a3m"))[0]
                    output_yaml = split_colab_a3m_write_yaml(output_msa)
                    print(f"MSA generated at: {output_msa}")

        with open(output_yaml, "r") as f:
            yaml_content = yaml.safe_load(f)
    else:
        print("No protein chains found, skipping protein MSA generation.")
        yaml_content = {"a3m": []}

    # Write method log
    method_log = {"msa": args.msa, "seq": args.seq, "stoi": args.stoi}
    if yaml_content.get("templates"):
        method_log["templates"] = yaml_content["templates"]
    with open(os.path.join(msa_dir, "method_log.yaml"), "w") as f:
        yaml.dump(method_log, f)

    # Add NA chains to data yaml
    yaml_content["method_log"] = os.path.join(msa_dir, "method_log.yaml")
    if na_chains:
        from MSA.script.RNA_MSA_search.sto_to_a3m import convert as sto_to_a3m

        rna_msa_script = os.path.join(ROOT, "MSA", "script", "RNA_MSA_search")
        rna_db_dir = CONFIG["rna_msa_db_dir"]
        for i, na in enumerate(na_chains):
            na_fa_path = os.path.join(msa_dir, f"{target_name}_na_{i}.fa")
            with open(na_fa_path, "w") as f:
                f.write(f">{target_name}_na_{i}\n{na['sequence']}\n")

            # Run RNA MSA search and convert to a3m
            unpaired_path = os.path.abspath(na_fa_path)
            if na["type"] == "rna":
                na_sto_path = os.path.join(msa_dir, f"{target_name}_na_{i}.sto")
                na_a3m_path = os.path.join(msa_dir, f"{target_name}_na_{i}.a3m")
                # Reuse an existing RNA MSA when the same (msa, seq, stoi) already
                # produced a non-empty a3m (e.g. a prior structure-model run on
                # this target). The search is deterministic, so this only skips
                # redundant work.
                if (
                    params_match
                    and os.path.exists(na_a3m_path)
                    and os.path.getsize(na_a3m_path) > 0
                ):
                    print(
                        f"  RNA MSA already generated for chain {i} with the same "
                        f"parameters, skipping RNA MSA search."
                    )
                else:
                    print(f"  Running RNA MSA search for chain {i}...")
                    subprocess.run(
                        f"python {rna_msa_script}/msa_gen.py "
                        f"--query {na_fa_path} "
                        f"--db_dir {rna_db_dir} "
                        f"--output {na_sto_path}",
                        shell=True,
                        check=True,
                    )
                    sto_to_a3m(na_sto_path, na_a3m_path)
                if os.path.exists(na_a3m_path):
                    unpaired_path = os.path.abspath(na_a3m_path)
                else:
                    unpaired_path = (
                        na_fa_path  # Fallback to fasta if a3m generation fails
                    )

            yaml_content["a3m"].append(
                {
                    "paired_path": None,
                    "unpaired_path": unpaired_path,
                    "copy": na["copy"],
                    "type": na["type"],
                }
            )
    data_yaml = f"{base_dir}/{target_name}.yaml"
    with open(data_yaml, "w") as f:
        f.write("# Fill in the following fields before running structure prediction\n")
        f.write("# job_name:\n")
        f.write("# output_dir:\n")
        f.write("# seed:\n\n")
        yaml.dump(yaml_content, f, indent=2)

    return data_yaml


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--msa",
        type=str,
        required=True,
        choices=["colab"],
        help="MSA generation method to use",
    )
    parser.add_argument(
        "--seq",
        type=str,
        required=True,
        help="Path to the input sequence file in FASTA format",
    )
    parser.add_argument(
        "--stoi",
        type=str,
        required=True,
        help="Stoichiometry information, e.g. 'A1' for one chain A",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default=None,
        help="Output directory (default: same as input FASTA directory)",
    )
    args = parser.parse_args()
    msa_generation(args)
