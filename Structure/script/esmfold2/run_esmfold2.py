"""ESMFold2 inference runner for the Thal-Kak pipeline.

Reads a Thal-Kak data yaml (the one produced by MSA generation) plus an
ESMFold2 config yaml, builds an ``esm`` ``StructurePredictionInput`` from the
per-chain a3m files (paired + unpaired merged with ``key=`` taxonomy markers
for cross-chain pairing), then loops over ``data_yaml['seed']`` and writes
each diffusion sample as ``{target}_seed_{N}_sample_{M}.pdb`` into the
result directory. A ``{target}_results_summary.csv`` mirroring the schema of
the other model runners is also written, and a ``RESULT_DIR:<path>`` line is
emitted on stdout so the orchestrator can pick it up.
"""

import argparse
import csv
import os
import string
import sys
from datetime import datetime
from pathlib import Path

import gemmi
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

# Silence rdkit's pickle-version warning spammed while loading ccd.pkl
# (built with a newer rdkit than the pinned 2024.9.6).
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

from esm.models.esmfold2.processor import ESMFold2InputBuilder
from esm.utils.msa.msa import MSA
from esm.utils.parsing import FastaEntry, read_sequences
from esm.utils.structure.input_builder import (
    DNAInput,
    LigandInput,
    ProteinInput,
    RNAInput,
    StructurePredictionInput,
)
from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
COMMON_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "common")
if COMMON_DIR not in sys.path:
    sys.path.insert(0, COMMON_DIR)
from chain_utils import assign_chain_indices


def _chain_letter(index: int) -> str:
    """Convert 0-based chain index to a chain id (A..Z, then AA, BA, ...)."""
    if index < 26:
        return chr(ord("A") + index)
    return chr(ord("A") + index % 26) + chr(ord("A") + index // 26)


def _chain_breaks(chain_ids):
    """Token positions where the chain id transitions, matching the convention
    used by the other model wrappers' pAE/pLDDT plots."""
    arr = np.asarray(chain_ids)
    return [i + 1 for i in range(len(arr) - 1) if arr[i] != arr[i + 1]]


def _plot_pae_matrix(ax, pae, breaks, title):
    pae = np.asarray(pae)
    num_res = pae.shape[0]
    im = ax.imshow(pae, cmap="bwr", vmin=0, vmax=30, extent=(0, num_res, num_res, 0))
    ax.set_xlim(0, num_res)
    ax.set_ylim(num_res, 0)
    ax.set_aspect("equal")
    for b in breaks:
        ax.axvline(x=b, color="black", linewidth=1.0)
        ax.axhline(y=b, color="black", linewidth=1.0)
    chain_boundaries = [0] + list(breaks) + [num_res]
    ytick_positions, ytick_labels = [], []
    for i in range(len(chain_boundaries) - 1):
        ytick_positions.append((chain_boundaries[i] + chain_boundaries[i + 1]) / 2)
        ytick_labels.append(string.ascii_uppercase[i])
    ax.set_yticks(ytick_positions)
    ax.set_yticklabels(ytick_labels, fontsize=10, fontweight="bold")
    ax.set_title(title, fontsize=9)
    return im


def _save_confidence_plots(common_dir, target, plot_data):
    """Save PAE matrix grid + pLDDT comparison PNGs for the top-5 samples."""
    if not plot_data:
        return
    plot_data = sorted(plot_data, key=lambda d: d["rank_key"], reverse=True)[:5]

    fig_pae, axes = plt.subplots(1, len(plot_data), figsize=(4 * len(plot_data) + 1, 4))
    if len(plot_data) == 1:
        axes = [axes]
    im = None
    for ax, d in zip(axes, plot_data):
        if d["pae"] is None:
            ax.set_axis_off()
            continue
        im = _plot_pae_matrix(ax, d["pae"], d["breaks"], f"{d['rank_label']}\n{d['stem']}")
    if im is not None:
        fig_pae.subplots_adjust(right=0.85)
        cbar_ax = fig_pae.add_axes([0.88, 0.15, 0.02, 0.7])
        fig_pae.colorbar(im, cax=cbar_ax, label="predicted Expected Error")
        fig_pae.savefig(
            os.path.join(common_dir, f"{target}_pae_rankings.png"),
            dpi=300,
            bbox_inches="tight",
        )
    plt.close(fig_pae)

    plt.figure(figsize=(12, 5))
    drew = False
    for d in plot_data:
        if d["plddt"] is None:
            continue
        plt.plot(d["plddt"], label=d["rank_label"], alpha=0.8)
        if d["rank_label"] == "rank_001":
            for b in d["breaks"]:
                plt.axvline(x=b, color="black", linestyle="--", linewidth=0.7, alpha=0.4)
        drew = True
    if drew:
        plt.title(f"{target} pLDDT Comparison")
        plt.xlabel("Residue Index")
        plt.ylabel("pLDDT")
        plt.ylim(0, 100)
        plt.legend()
        plt.grid(True, axis="y", alpha=0.2)
        plt.savefig(
            os.path.join(common_dir, f"{target}_plddt_comparison.png"),
            dpi=300,
            bbox_inches="tight",
        )
    plt.close("all")


def _build_chain_msa(paired_path, unpaired_path) -> MSA:
    """Merge our split paired/unpaired a3m files into one ESMFold2 MSA with
    ``key=`` taxonomy markers. Both source files include the query at row 0,
    so we keep it once and drop the duplicate."""
    paired_rows = list(read_sequences(paired_path)) if paired_path else []
    unpaired_rows = list(read_sequences(unpaired_path)) if unpaired_path else []

    if not paired_rows and not unpaired_rows:
        raise ValueError("Both paired_path and unpaired_path are empty")

    entries: list[FastaEntry] = []
    query_header, query_seq = (paired_rows or unpaired_rows)[0]
    entries.append(FastaEntry(query_header, query_seq))

    # Paired rows past the query: same row index across chains → cross-chain pair
    for i, (header, seq) in enumerate(paired_rows[1:], start=1):
        entries.append(FastaEntry(f"{header} key={i}", seq))

    # Unpaired rows past the query: key=-1 → block-diagonal unpaired region
    for header, seq in unpaired_rows[1:]:
        entries.append(FastaEntry(f"{header} key=-1", seq))

    return MSA(entries)


def _make_polymer_input(entry, chain_ids, use_msa):
    """Build a ProteinInput / RNAInput / DNAInput from one a3m entry."""
    mol_type = entry.get("type", "protein")
    paired = entry.get("paired_path")
    unpaired = entry.get("unpaired_path")

    # Sequence: first row of either a3m (both have query at row 0)
    src = paired or unpaired
    if not src:
        raise ValueError(f"a3m entry has no paired_path or unpaired_path: {entry}")
    _, query_seq = next(iter(read_sequences(src)))

    if mol_type == "protein":
        msa = _build_chain_msa(paired, unpaired) if use_msa else None
        return ProteinInput(id=chain_ids, sequence=query_seq, msa=msa)
    if mol_type == "rna":
        return RNAInput(id=chain_ids, sequence=query_seq)
    if mol_type == "dna":
        return DNAInput(id=chain_ids, sequence=query_seq)
    raise ValueError(f"Unsupported a3m entry type: {mol_type!r}")


def _make_ligand_input(entry, chain_ids):
    if entry.get("smiles"):
        return LigandInput(id=chain_ids, smiles=entry["smiles"])
    if entry.get("ccd"):
        ccd = entry["ccd"]
        return LigandInput(id=chain_ids, ccd=ccd if isinstance(ccd, list) else [ccd])
    raise ValueError(f"Ligand entry has neither smiles nor ccd: {entry}")


def _result_root(output_dir, target, job_name):
    root = os.path.join(output_dir, f"esmfold2_results_{target}_{job_name}")
    if os.path.exists(root):
        root += datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")
    os.makedirs(os.path.join(root, "common"))
    return root


def main(data_yaml_path, esm_yaml_path):
    with open(data_yaml_path) as f:
        data_cfg = yaml.safe_load(f)
    with open(esm_yaml_path) as f:
        esm_cfg = yaml.safe_load(f)

    target = data_cfg.get("name") or Path(data_yaml_path).stem
    job_name = data_cfg["job_name"]
    output_dir = data_cfg["output_dir"]
    seeds = data_cfg["seed"]
    if not isinstance(seeds, list):
        seeds = [seeds]

    a3m_entries = data_cfg.get("a3m", [])
    ligand_entries = data_cfg.get("ligand", [])

    model_variant = esm_cfg.get("model_variant", "biohub/ESMFold2")
    num_loops = int(esm_cfg.get("num_loops", 3))
    num_sampling_steps = int(esm_cfg.get("num_sampling_steps", 200))
    num_diffusion_samples = int(esm_cfg.get("num_diffusion_samples", 5))
    # MSA is only used by the full ESMFold2 variant; the -Fast variant skips it.
    is_fast = model_variant.endswith("-Fast")
    use_msa = esm_cfg.get("use_msa", True) and not is_fast

    # Chain id assignment: protein round-robin first, NA next, ligands last
    # (matches other models; keeps template chain_query (protein-only at MSA
    # time) aligned with the runner's chain pool).
    polymer_copies = [entry.get("copy", 1) for entry in a3m_entries]
    polymer_types = [
        0 if entry.get("type", "protein") == "protein" else 1
        for entry in a3m_entries
    ]
    chains_per_entity = assign_chain_indices(polymer_copies, polymer_types)
    n_polymer = sum(polymer_copies)
    chains_per_ligand = []
    k = n_polymer
    for entry in ligand_entries:
        c = entry.get("copy", 1)
        chains_per_ligand.append(list(range(k, k + c)))
        k += c

    sequences = []
    for idx, entry in enumerate(a3m_entries):
        ids = [_chain_letter(i) for i in chains_per_entity[idx]]
        sequences.append(_make_polymer_input(entry, ids, use_msa))
    for j, entry in enumerate(ligand_entries):
        ids = [_chain_letter(i) for i in chains_per_ligand[j]]
        sequences.append(_make_ligand_input(entry, ids))

    spi = StructurePredictionInput(sequences=sequences)

    print(f"[esmfold2] loading {model_variant} ...", flush=True)
    model = ESMFold2Model.from_pretrained(model_variant).cuda().eval()
    builder = ESMFold2InputBuilder()

    result_root = _result_root(output_dir, target, job_name)
    common = os.path.join(result_root, "common")
    cif_dir = os.path.join(result_root, "cif_outputs")
    os.makedirs(cif_dir, exist_ok=True)

    summary_rows = []
    plot_data = []
    for seed in seeds:
        print(f"[esmfold2] seed={seed} folding ...", flush=True)
        results = builder.fold(
            model,
            spi,
            num_loops=num_loops,
            num_sampling_steps=num_sampling_steps,
            num_diffusion_samples=num_diffusion_samples,
            seed=int(seed),
            complex_id=target,
        )
        if not isinstance(results, list):
            results = [results]

        for m, res in enumerate(results):
            stem = f"{target}_seed_{seed}_sample_{m}"
            cif_path = os.path.join(cif_dir, f"{stem}.cif")
            pdb_path = os.path.join(common, f"{stem}.pdb")
            with open(cif_path, "w") as f:
                f.write(res.complex.to_mmcif())
            gemmi.read_structure(cif_path).write_pdb(pdb_path)
            # ESMFold2 emits plddt in [0, 1]; rescale to the 0-100 convention
            # the rest of the pipeline (and the comparison plot) expects.
            plddt_arr = (
                res.plddt.detach().cpu().numpy() * 100.0
                if res.plddt is not None
                else None
            )
            pae_arr = res.pae.detach().cpu().numpy() if res.pae is not None else None
            breaks = _chain_breaks(res.complex.chain_id)
            mean_plddt = float(plddt_arr.mean()) if plddt_arr is not None else None
            ptm = float(res.ptm) if res.ptm is not None else None
            iptm = float(res.iptm) if res.iptm is not None else None
            # Ranking score: AF3-style 0.8*iptm + 0.2*ptm when available, else ptm
            ranking = (
                0.8 * iptm + 0.2 * ptm
                if (iptm is not None and ptm is not None)
                else (ptm if ptm is not None else mean_plddt)
            )
            summary_rows.append(
                {
                    "target": target,
                    "option": job_name,
                    "model": "esmfold2",
                    "seed-sample": f"seed_{seed}_sample_{m}",
                    "ranking_score": (
                        round(ranking, 4) if ranking is not None else ""
                    ),
                    "mean_plddt": round(mean_plddt, 3) if mean_plddt is not None else "",
                    "ptm": round(ptm, 3) if ptm is not None else "",
                    "iptm": round(iptm, 3) if iptm is not None else "",
                }
            )
            plot_data.append(
                {
                    "stem": stem,
                    "rank_key": ranking if ranking is not None else float("-inf"),
                    "plddt": plddt_arr,
                    "pae": pae_arr,
                    "breaks": breaks,
                }
            )

    summary_path = os.path.join(common, f"{target}_results_summary.csv")
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "target",
                "option",
                "model",
                "seed-sample",
                "ranking_score",
                "mean_plddt",
                "ptm",
                "iptm",
            ],
        )
        writer.writeheader()
        # Sort by ranking_score desc so consumers that read top-N get best first
        summary_rows.sort(
            key=lambda r: float(r["ranking_score"]) if r["ranking_score"] != "" else -1,
            reverse=True,
        )
        writer.writerows(summary_rows)

    plot_data.sort(key=lambda d: d["rank_key"], reverse=True)
    for i, d in enumerate(plot_data, 1):
        d["rank_label"] = f"rank_{i:03d}"
    _save_confidence_plots(common, target, plot_data)

    print(f"RESULT_DIR:{result_root}")
    return result_root


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_yaml", type=str, required=True)
    parser.add_argument("--esmfold2_yaml", type=str, required=True)
    args = parser.parse_args()
    main(args.data_yaml, args.esmfold2_yaml)
