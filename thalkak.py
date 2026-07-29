import argparse, yaml, os, csv, re, glob, shutil, logging, subprocess, sys, io, contextlib
from collections import namedtuple

ROOT = os.path.dirname(os.path.abspath(__file__))


def _select_top5_for_job(decoy_dir, top5_dir):
    """Pick the top-5 structure-validated PDBs for one prediction job, ranked
    by the model's own confidence (``*_results_summary.csv`` written next to the
    decoys), copy them as model_[1-5].pdb into ``top5_dir`` and write a
    method_log.yaml summarising the picks."""
    from Structure.script.common.structure_validation import validation

    summary_csvs = glob.glob(os.path.join(decoy_dir, "*_results_summary.csv"))
    if not summary_csvs:
        raise FileNotFoundError(f"No *_results_summary.csv found in {decoy_dir}")
    with open(summary_csvs[0]) as f:
        rows = sorted(
            csv.DictReader(f), key=lambda r: float(r["ranking_score"]), reverse=True
        )

    os.makedirs(top5_dir, exist_ok=True)

    picked = []
    for row in rows:
        if len(picked) >= 5:
            break
        candidates = glob.glob(os.path.join(decoy_dir, f"*{row['seed-sample']}*.pdb"))
        if not candidates:
            continue
        candidate = candidates[0]
        if not validation(candidate):
            continue
        picked.append(candidate)

    method_log = {"models": {}}
    for i, src in enumerate(picked, 1):
        shutil.copy(src, os.path.join(top5_dir, f"model_{i}.pdb"))
        entry = {}
        ml_src = os.path.join(os.path.dirname(src), "method_log.yaml")
        if os.path.exists(ml_src):
            with open(ml_src) as f:
                entry.update(yaml.safe_load(f) or {})
        m = re.search(r"_seed_(\d+)_sample_(\d+)\.pdb$", os.path.basename(src))
        if m:
            entry["seed"] = int(m.group(1))
            entry["sample"] = int(m.group(2))
        method_log["models"][f"model_{i}"] = entry

    with open(os.path.join(top5_dir, "method_log.yaml"), "w") as f:
        yaml.dump(method_log, f, sort_keys=False)

    if len(picked) < 5:
        print(
            f"  WARN: only {len(picked)} valid candidate(s) for "
            f"{os.path.basename(top5_dir)} (wanted 5)"
        )
    return picked


def run_full(args):
    from MSA.msa_generation import msa_generation
    from Structure.structure_prediction import structure_prediction
    from Relax.relaxation import relaxation

    print(
        f"Running Thal-Kak with MSA method: {args.msa}, structure prediction "
        f"model: {args.structure}, input sequence file: {args.seq}"
    )
    base_dir = args.base_dir if args.base_dir else os.path.dirname(args.seq)
    os.makedirs(base_dir, exist_ok=True)

    # MSA generation
    if args.stoi == "UNK":
        with open(args.seq, "r") as f:
            num_chains = len(f.read().strip().split("\n")[1::2])
        args.stoi = "".join(f"{chr(65 + i)}1" for i in range(num_chains))
        print(f"STOI is set to UNK, inferred as {args.stoi}")
    msa_args = namedtuple("MsaArgs", ["msa", "seq", "stoi", "output_dir"])(
        msa=args.msa,
        seq=args.seq,
        stoi=args.stoi,
        output_dir=base_dir,
    )
    data_yaml = msa_generation(msa_args)

    # Edit data yaml for structure prediction
    with open(data_yaml, "r") as f:
        yaml_content = yaml.safe_load(f)
    yaml_content["job_name"] = f"{args.msa}_{args.structure}"
    yaml_content["output_dir"] = os.path.join(base_dir, "structure")
    yaml_content["seed"] = list(range(args.seed_start, args.seed_start + args.n_seed))
    if args.ligand_yaml:
        with open(args.ligand_yaml, "r") as f:
            ligand_cfg = yaml.safe_load(f) or {}
        if ligand_cfg.get("ligand"):
            yaml_content["ligand"] = ligand_cfg["ligand"]
            print(
                f"Merged {len(ligand_cfg['ligand'])} ligand entry/entries "
                f"from {args.ligand_yaml}"
            )
    with open(data_yaml, "w") as f:
        yaml.dump(yaml_content, f, indent=2)

    # Structure prediction
    print("Running structure prediction...")
    sp_args = namedtuple("SpArgs", ["model", "data_config", "model_config"])(
        model=args.structure,
        data_config=data_yaml,
        model_config=(
            args.model_config
            if args.model_config
            else os.path.join(ROOT, "examples", f"{args.structure}.yaml")
        ),
    )
    result_dir = structure_prediction(sp_args)
    print(f"Structure prediction results saved at: {result_dir}")

    # Per-job top-5 selection by the model's own confidence. Job folder name
    # mirrors the structure result dir (timestamp suffix preserved if any).
    decoy_dir = os.path.join(result_dir, "common")
    job_name = os.path.basename(result_dir)
    top5_dir = os.path.join(base_dir, "top5", job_name)
    print(f"Selecting top-5 for {job_name}...")
    _select_top5_for_job(decoy_dir, top5_dir)
    print(f"Top-5 saved at: {top5_dir}")

    # Relax the top-5 in place.
    print(f"Running relaxation ({args.relax}) on top-5...")
    relax_dir = relaxation(top5_dir, args.relax, args.relax_config)
    print(f"Relaxation complete: {relax_dir}")


def cli():
    setup_logging()
    parser = argparse.ArgumentParser(description="Thal-Kak structure prediction pipeline")
    subparsers = parser.add_subparsers(dest="mode", required=True, help="Pipeline mode")

    # full
    p_full = subparsers.add_parser("full", help="Run the full pipeline")
    p_full.add_argument(
        "--msa", type=str, required=True, choices=["colab"],
        help="MSA generation method",
    )
    p_full.add_argument(
        "--structure", type=str, required=True,
        choices=["boltz2", "chai1", "protenix", "esmfold2"],
        help="Structure prediction model",
    )
    p_full.add_argument(
        "--relax", type=str, required=True, choices=["none", "openmm"],
        help="Relaxation method",
    )
    p_full.add_argument(
        "--seq", type=str, required=True, help="Path to input FASTA file"
    )
    p_full.add_argument(
        "--stoi", type=str, required=True, help="Stoichiometry, e.g. 'A1'"
    )
    p_full.add_argument(
        "--model_config", type=str, default=None,
        help="Model config yaml (default: examples/{structure}.yaml)",
    )
    p_full.add_argument(
        "--base_dir", type=str, default=None,
        help="Output directory (default: same as input FASTA directory)",
    )
    p_full.add_argument(
        "--n_seed", type=int, default=5,
        help="Number of seeds for structure prediction",
    )
    p_full.add_argument(
        "--seed_start", type=int, default=1,
        help="First seed value (seeds are seed_start .. seed_start + n_seed - 1)",
    )
    p_full.add_argument(
        "--ligand_yaml", type=str, default=None,
        help="Optional YAML file with a 'ligand' list to merge into the data config",
    )
    p_full.add_argument(
        "--relax_config", type=str, default=None,
        help="Relax config yaml (default: examples/{relax}.yaml)",
    )

    # msa
    p_msa = subparsers.add_parser("msa", help="Run MSA generation only")
    p_msa.add_argument(
        "--msa", type=str, required=True, choices=["colab"],
        help="MSA generation method",
    )
    p_msa.add_argument(
        "--seq", type=str, required=True, help="Path to input FASTA file"
    )
    p_msa.add_argument(
        "--stoi", type=str, required=True, help="Stoichiometry, e.g. 'A1'"
    )
    p_msa.add_argument(
        "--output_dir", type=str, default=None,
        help="Output directory (default: same as input FASTA directory)",
    )

    # structure
    p_sp = subparsers.add_parser("structure", help="Run structure prediction only")
    p_sp.add_argument(
        "--model", type=str, required=True,
        choices=["boltz2", "chai1", "protenix", "esmfold2"],
        help="Structure prediction model",
    )
    p_sp.add_argument(
        "--data_config", type=str, required=True, help="Path to data config yaml"
    )
    p_sp.add_argument(
        "--model_config", type=str, default=None,
        help="Model config yaml (default: examples/{model}.yaml)",
    )

    # relax
    p_relax = subparsers.add_parser("relax", help="Run relaxation only")
    p_relax.add_argument(
        "--decoy_dir", type=str, required=True, help="Directory containing PDB files"
    )
    p_relax.add_argument(
        "--relax", type=str, required=True, choices=["none", "openmm"],
        help="Relaxation method",
    )
    p_relax.add_argument(
        "--relax_config", type=str, default=None,
        help="Relax config yaml (default: examples/{relax}.yaml)",
    )

    args = parser.parse_args()

    match args.mode:
        case "full":
            run_full(args)
        case "msa":
            from MSA.msa_generation import msa_generation
            msa_generation(args)
        case "structure":
            from Structure.structure_prediction import structure_prediction
            if not args.model_config:
                args.model_config = os.path.join(ROOT, "examples", f"{args.model}.yaml")
            structure_prediction(args)
        case "relax":
            from Relax.relaxation import relaxation
            relaxation(args.decoy_dir, args.relax, args.relax_config)


# =============================== logging ===============================
# Only the thalkak.* tree uses this setup, so third-party INFO (numexpr, jax, ...)
# stays out. Tag column = stage for INFO, level for WARNING+. External tools are
# funnelled through the same logger via run_logged (subprocess) / log_stream.

_PKG = "thalkak"
_SUB = "  │ "
_DATEFMT = "%Y-%m-%d %H:%M:%S"


class _Fmt(logging.Formatter):
    def format(self, record):
        record.tag = (record.levelname if record.levelno >= logging.WARNING
                      else record.name.rsplit(".", 1)[-1])
        return super().format(record)


def setup_logging(level=logging.INFO):
    """Configure the thalkak.* logger to stdout. Idempotent."""
    os.environ.setdefault("PYTHONUNBUFFERED", "1")  # child tools stream line-by-line
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except (AttributeError, ValueError):
        pass
    pkg = logging.getLogger(_PKG)
    if pkg.handlers:
        return pkg
    pkg.setLevel(level)
    pkg.propagate = False
    h = logging.StreamHandler(sys.stdout)
    h.setFormatter(_Fmt("%(asctime)s | %(tag)-9s | %(message)s", _DATEFMT))
    pkg.addHandler(h)
    return pkg


def get_logger(name):
    """Stage logger; `name` (e.g. 'msa') fills the tag column."""
    return logging.getLogger(f"{_PKG}.{name}")


def section(log, title, width=60):
    log.info(f" {title} ".center(width, "="))


def _emit(log, level, line):
    """Log one line of external-tool output: keep only the final \\r-overwrite
    (tqdm bars), strip trailing space, indent under _SUB, skip if blank."""
    line = line.rsplit("\r", 1)[-1].rstrip()
    if line:
        log.log(level, "%s%s", _SUB, line)


def run_logged(cmd, log=None, check=True, **kw):
    """subprocess.run-like, but stream the merged stdout+stderr through `log` one
    timestamped line at a time. Returns exit code; raises on non-zero if `check`."""
    log = log or get_logger("subprocess")
    shell = isinstance(cmd, str)
    log.info("$ %s", cmd if shell else " ".join(map(str, cmd)))
    proc = subprocess.Popen(cmd, shell=shell, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, bufsize=1, **kw)
    for line in proc.stdout:
        _emit(log, logging.INFO, line)
    code = proc.wait()
    if code:
        log.error("↳ command exited with status %d", code)
        if check:
            raise subprocess.CalledProcessError(code, cmd)
    return code


def log_lines(text, log=None, level=logging.INFO):
    """Emit already-captured text line-by-line (parallel jobs, kept contiguous)."""
    log = log or get_logger("subprocess")
    for line in (text or "").splitlines():
        _emit(log, level, line)


class _LineWriter(io.TextIOBase):
    """File-like shim: forward complete lines to a logger via _emit; delegate
    fileno/isatty to the real stream so libraries probing stdout don't crash."""
    def __init__(self, log, level, stream):
        self._log, self._level, self._stream, self._buf = log, level, stream, ""

    def write(self, s):
        self._buf += s
        while "\n" in self._buf:
            line, self._buf = self._buf.split("\n", 1)
            _emit(self._log, self._level, line)
        return len(s)

    def flush(self):
        _emit(self._log, self._level, self._buf)
        self._buf = ""

    def fileno(self):
        return self._stream.fileno()

    def isatty(self):
        return False


@contextlib.contextmanager
def log_stream(log=None, level=logging.INFO):
    """Route sys.stdout/stderr through `log` inside the block (imported in-process
    steps: boltz/chai/esmfold/protenix). C-level fd writes still print raw."""
    log = log or get_logger("subprocess")
    out, err = _LineWriter(log, level, sys.stdout), _LineWriter(log, level, sys.stderr)
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            yield
        finally:
            out.flush()
            err.flush()

# =======================================================================


if __name__ == "__main__":
    cli()
