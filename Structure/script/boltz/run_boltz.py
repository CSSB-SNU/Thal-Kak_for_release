import yaml, subprocess, os, json, argparse, shutil, time, sys
from datetime import datetime
from pathlib import Path
import pandas as pd
from collections import Counter
import shlex

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
COMMON_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "common")
for _p in (SCRIPT_DIR, COMMON_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from chain_utils import assign_chain_indices
from template_cleaner import clean_template_for_boltz


def auth_to_label(cif_path, auth_chain):
    # Convert auth_chain to label_chain using cif file
    tags = []
    i_auth = i_label = None
    in_atom_loop = False
    cnt = Counter()

    with open(cif_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            if s == "loop_":
                tags = []
                in_atom_loop = False
                i_auth = i_label = None
                continue

            # collect atom_site tags
            if s.startswith("_atom_site."):
                tags.append(s.split()[0])
                continue

            # first data line after atom_site tags
            if tags and not in_atom_loop:
                try:
                    i_auth = tags.index("_atom_site.auth_asym_id")
                    i_label = tags.index("_atom_site.label_asym_id")
                    in_atom_loop = True
                except ValueError:
                    tags = []  # not the loop we want
                    continue

            if in_atom_loop:
                # end of atom_site loop
                if (
                    s.startswith("_")
                    or s == "loop_"
                    or s.startswith("data_")
                    or s.startswith("#")
                ):
                    break
                row = shlex.split(s)
                if len(row) > max(i_auth, i_label) and row[i_auth] == auth_chain:
                    cnt[row[i_label]] += 1

    if not cnt:
        raise ValueError(f"auth chain {auth_chain} not found")
    return cnt.most_common(1)[0][0]


def fix_pdb_chain_ids(pdb_path):
    """Fix 2-character chain IDs (e.g., AA, BA) to lowercase single chars (a, b, ...)."""
    with open(pdb_path, "r") as f:
        lines = f.readlines()

    # Collect 2-char chain IDs in order of first appearance
    seen = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and line[22].isalpha():
            chain_id = line[21:23]
            if chain_id not in seen:
                seen.append(chain_id)

    if not seen:
        return

    chain_map = {cid: chr(97 + i) for i, cid in enumerate(seen)}

    fixed_lines = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and line[22].isalpha():
            chain_id = line[21:23]
            line = line[:21] + chain_map[chain_id] + line[23:]
        fixed_lines.append(line)

    with open(pdb_path, "w") as f:
        f.writelines(fixed_lines)


def reorder_pdb_chains(pdb_path):
    """Reorder atom records so chains appear alphabetically (A, B, C, D).

    boltz groups chain copies by entity, so an A2B2 complex comes out as
    A, C, B, D, whereas the other models emit A, B, C, D. Reorder the chain
    blocks to match and renumber atom serials so they stay monotonic. Chain
    labels are unchanged -- only record order and serials.
    """
    with open(pdb_path, "r") as f:
        lines = f.readlines()

    blocks, order, header, trailer = {}, [], [], []
    seen_atom = False
    for line in lines:
        if line.startswith(("ATOM", "HETATM", "ANISOU", "TER")):
            seen_atom = True
            chain_id = line[21:22]
            if chain_id not in blocks:
                blocks[chain_id] = []
                order.append(chain_id)
            blocks[chain_id].append(line)
        elif not seen_atom:
            header.append(line)
        else:
            trailer.append(line)

    if order == sorted(order):
        return  # already alphabetical

    out = list(header)
    serial = 0
    for chain_id in sorted(blocks):
        for line in blocks[chain_id]:
            if line.startswith("ANISOU"):
                s = serial  # ANISOU shares the serial of its preceding atom
            else:
                serial += 1
                s = serial
            out.append(f"{line[:6]}{s:>5}{line[11:]}")
    out.extend(trailer)

    with open(pdb_path, "w") as f:
        f.writelines(out)


def rename_output_dir(output_dir):
    if os.path.exists(f"{output_dir}"):
        output_dir += datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")
    return output_dir


def mv_output_dir(temp_dir, output_dir):
    os.makedirs(output_dir)
    os.rename(temp_dir, f"{output_dir}")
    os.makedirs(f"{output_dir}/common", exist_ok=True)
    try:
        os.rmdir("./temp/")
    except:
        pass


def main(data_yaml, boltz2_yaml):
    with open(data_yaml, "r") as file:
        data_config = yaml.safe_load(file)
    with open(boltz2_yaml, "r") as file:
        boltz_config = yaml.safe_load(file)

    ### Make input.yaml and csv files
    ## initialize
    name = os.path.basename(data_yaml).split(".")[0]
    temp_dir = f"temp/boltz_results_{name}"
    # Preserve leftover state from a previous failed run -- stale msa CSVs,
    # cleaned templates, or input.yaml under temp_dir would otherwise leak
    # into this job. Rename with a timestamp so debug evidence is kept.
    if os.path.exists(temp_dir):
        failed_dir = f"{temp_dir}_failed_{datetime.now():%Y%m%d_%H%M%S}"
        os.rename(temp_dir, failed_dir)
        print(f"[boltz] preserved previous temp at {failed_dir}")
    os.makedirs(f"{temp_dir}/msa")
    yaml_output = "version: 1\nsequences:\n"  # Initialize YAML output
    n_chain = 0

    ## a3m parsing
    chain_len, chain_copy = [], []
    protein_chain_ids = set()
    a3m_chain_indices = assign_chain_indices(
        [e["copy"] for e in data_config["a3m"]],
        [0 if e.get("type", "protein") == "protein" else 1 for e in data_config["a3m"]],
    )
    n_chain = sum(e["copy"] for e in data_config["a3m"])
    for i, entity in enumerate(data_config["a3m"]):
        paired_a3m, unpaired_a3m = [], []
        if entity["paired_path"] != None:
            with open(entity["paired_path"], "r") as f:
                paired_a3m = f.readlines()
        if entity["unpaired_path"] != None:
            with open(entity["unpaired_path"], "r") as f:
                unpaired_a3m = f.readlines()
        sequence = unpaired_a3m[1].strip()
        chain_len.append(len(sequence))
        chain_copy.append(entity["copy"])

        if entity["type"] == "protein":
            ## Write msa.csv files
            with open(f"{temp_dir}/msa/{name}_msa_{i}.csv", "w") as csv_file:
                write_data = "key,sequence\n"
                for j, seq in enumerate(paired_a3m[1::2]):
                    write_data += f"{j},{seq.strip()}\n"
                for seq in unpaired_a3m[1::2]:
                    write_data += f"-1,{seq.strip()}\n"
                csv_file.write(write_data[:-1])  # remove last newline

        ## Append input.yaml
        chain_ids = [
            chr(65 + k) if k < 26 else chr(65 + k % 26) + chr(64 + k // 26)
            for k in a3m_chain_indices[i]
        ]
        yaml_output += f"  - {entity['type']}:\n"
        yaml_output += f"      id: [{','.join(chain_ids)}]\n"
        yaml_output += f"      sequence: {sequence}\n"
        if entity["type"] == "protein":
            protein_chain_ids.update(chain_ids)
            yaml_output += (
                f"      msa: {Path.cwd() / Path(temp_dir)}/msa/{name}_msa_{i}.csv\n"
            )

    ## ligand parsing
    if "ligand" in data_config and data_config["ligand"]:
        for entity in data_config["ligand"]:
            yaml_output += f"  - ligand:\n"
            yaml_output += f"      id: [{','.join([chr(65 + k) if k < 26 else chr(65 + k%26)+chr(64 + k//26) for k in range(n_chain, n_chain + entity['copy'])])}]\n"
            if "smiles" in entity and entity["smiles"]:
                yaml_output += f"      smiles: '{entity['smiles']}'\n"
            elif "ccd" in entity and entity["ccd"]:
                yaml_output += f"      ccd: '{entity['ccd']}'\n"
            n_chain += entity["copy"]

    ## template parsing
    if "templates" in data_config and data_config["templates"]:
        cleaned_tmpl_dir = Path(temp_dir) / "cleaned_templates"
        cleaned_tmpl_dir.mkdir(parents=True, exist_ok=True)
        template_yaml = ""
        for template in data_config["templates"]:
            # boltz only supports templates on protein chains. Drop any
            # non-protein (RNA/DNA) query chains, else boltz raises
            # "Chain X assigned for template ... is not one of the protein chains!"
            query_chains = [
                c for c in template["chain_query"] if c in protein_chain_ids
            ]
            if not query_chains:
                print(
                    f"[boltz] skipping non-protein template "
                    f"{Path(template['path']).name} "
                    f"(query chains {list(template['chain_query'])})"
                )
                continue
            src = Path(template["path"])
            # Assume chain_template has only one unique value
            chain = template["chain_template"][0]
            cleaned = cleaned_tmpl_dir / f"{src.stem}_{chain}.cif"
            label_chain, n_renamed = clean_template_for_boltz(src, cleaned, chain)
            if n_renamed:
                print(
                    f"[boltz] {src.name}/{chain}: replaced {n_renamed} unknown "
                    f"residue(s) with UNK/N/DN"
                )
            template_yaml += f"  - cif: {cleaned}\n"
            template_yaml += f"    template_id: [{','.join([label_chain] * len(query_chains))}]\n"
            template_yaml += f"    chain_id: [{','.join(query_chains)}]\n"
        if template_yaml:
            yaml_output += "templates:\n" + template_yaml

    ## constraint parsing
    if "constraints" in boltz_config and boltz_config["constraints"]:
        yaml_output += "constraints:\n  "
        yaml_dump = yaml.dump(boltz_config["constraints"], default_flow_style=False)
        yaml_output += yaml_dump.replace("\n", "\n  ")

    with open(f"{temp_dir}/{name}.yaml", "w") as input_file:
        input_file.write(yaml_output)

    ### Run boltz2

    if boltz_config['no_kernels']:
        command = f"boltz predict {temp_dir}/{name}.yaml --out_dir {Path(temp_dir).parents[0]} --output_format {boltz_config['output_format']} --diffusion_samples {boltz_config['n_samples']} --recycling_steps {boltz_config['recycling_steps']} --sampling_steps {boltz_config['sampling_steps']} --no_kernels"
    else:
        command = f"boltz predict {temp_dir}/{name}.yaml --out_dir {Path(temp_dir).parents[0]} --output_format {boltz_config['output_format']} --diffusion_samples {boltz_config['n_samples']} --recycling_steps {boltz_config['recycling_steps']} --sampling_steps {boltz_config['sampling_steps']}"
    # boltz's native MSA subsampler (random num_subsampled_msa=1024 rows by
    # default). Off/absent = full MSA.
    if boltz_config.get("subsample_msa", False):
        command += " --subsample_msa"
    output_dir = (
        f"{data_config['output_dir']}/boltz_results_{name}_{data_config['job_name']}"
    )
    if data_config["seed"] != None:
        base_command = command
        for seed in data_config["seed"]:
            subprocess.run(
                f"{base_command} --seed {seed}",
                shell=True,
                check=True,
            )
            os.rename(
                f"{temp_dir}/predictions/{name}/",
                f"{temp_dir}/predictions/{name}_seed_{seed}/",
            )
    else:
        subprocess.run(
            command,
            shell=True,
            check=True,
        )
    output_dir = rename_output_dir(output_dir)
    mv_output_dir(temp_dir, output_dir)

    ### Organize output metrics
    ## Create metrics csv
    csv_rows = []
    for seed in data_config["seed"] if data_config["seed"] != None else [None]:
        for i in range(boltz_config["n_samples"]):
            if seed is not None:
                prediction_dir = f"{output_dir}/predictions/{name}_seed_{seed}/"
            else:
                prediction_dir = f"{output_dir}/predictions/{name}/"

            with open(
                f"{prediction_dir}/confidence_{name}_model_{i}.json", "r"
            ) as json_file:
                confidence_json = json.load(json_file)
            csv_rows.append(
                {
                    "target": name,
                    "option": data_config["job_name"],
                    "model": "boltz",
                    "seed-sample": f"seed_{seed}_sample_{i}",
                    "mean_plddt": f"{confidence_json['complex_plddt']*100:.3f}",
                    "ptm": f"{confidence_json['ptm']:.3f}",
                    "iptm": f"{confidence_json['iptm']:.3f}",
                    "ranking_score": f"{confidence_json['confidence_score']:.3f}",
                }
            )

            ## Copy output pdb files to common directory with seed info in filename
            shutil.copy(
                f"{prediction_dir}/{name}_model_{i}.pdb",
                f"{output_dir}/common/{name}_seed_{seed}_sample_{i}.pdb",
            )
            fix_pdb_chain_ids(f"{output_dir}/common/{name}_seed_{seed}_sample_{i}.pdb")
            reorder_pdb_chains(f"{output_dir}/common/{name}_seed_{seed}_sample_{i}.pdb")
    df = pd.DataFrame(csv_rows)
    df.sort_values(by=["ranking_score"], inplace=True, ascending=False)
    df.to_csv(f"{output_dir}/common/{name}_results_summary.csv", index=False)

    # for plotting
    chain_info = {"chain_len": chain_len, "chain_copy": chain_copy}
    with open(f"{output_dir}/chain_info.json", "w") as json_file:
        json.dump(chain_info, json_file)

    print(f"RESULT_DIR:{output_dir}")
    return output_dir


if __name__ == "__main__":
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_yaml", type=str, required=True)
    parser.add_argument("--boltz2_yaml", type=str, required=True)
    args = parser.parse_args()

    main(args.data_yaml, args.boltz2_yaml)

    end = time.time()
    print(f"Total time: {end - start:.2f} seconds")
