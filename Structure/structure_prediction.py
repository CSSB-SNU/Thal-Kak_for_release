import argparse, yaml, os, glob, sys
from datetime import datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def structure_prediction(args):
    with open(args.data_config) as f:
        data_yaml = yaml.safe_load(f)
    output_dir = data_yaml["output_dir"]
    target_name = os.path.basename(args.data_config).split(".")[0]
    job_name = data_yaml["job_name"]

    match args.model:
        case "boltz2":
            print("Running inference with Boltz2...")
            from Structure.script.boltz.run_boltz import main as run_boltz
            from Structure.script.boltz.boltz_confidence import main as boltz_confidence
            result_root = run_boltz(args.data_config, args.model_config)
            boltz_confidence(result_root, target_name)

        case "chai1":
            print("Running inference with Chai-1...")
            from Structure.script.chai.convert_yaml_to_json import convert_yaml_to_json
            from Structure.script.chai.run_chai import main as run_chai
            os.makedirs("temp/", exist_ok=True)

            data_json_path = f"temp/{target_name}.json"
            convert_yaml_to_json(args.data_config, data_json_path)

            model_name = os.path.basename(args.model_config).split(".")[0]
            model_json_path = f"temp/{target_name}_{model_name}.json"
            convert_yaml_to_json(args.model_config, model_json_path)

            result_root = run_chai(data_json_path, model_json_path)

            # move json files to result directory and clean up temp
            os.rename(data_json_path, f"{result_root}/{target_name}.json")
            os.rename(model_json_path, f"{result_root}/{model_name}.json")
            try:
                os.rmdir("temp/")
            except:
                pass

        case "esmfold2":
            print("Running inference with ESMFold2...")
            from Structure.script.esmfold2.run_esmfold2 import main as run_esmfold2
            result_root = run_esmfold2(args.data_config, args.model_config)

        case "protenix":
            print("Running inference with Protenix...")
            protenix_root = f"{ROOT}/Structure/submodules/protenix"
            seed = data_yaml["seed"]
            seed = ",".join(map(str, seed if isinstance(seed, list) else [seed]))

            result_root = f"{output_dir}/protenix_results_{target_name}_{job_name}"
            if os.path.exists(result_root):
                result_root += datetime.now().strftime("_%Y_%m_%d_%H_%M_%S")
            common_dir = f"{result_root}/common"
            os.makedirs(common_dir)

            os.environ["PROTENIX_ROOT_DIR"] = protenix_root
            os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
            os.environ["TQDM_DISABLE"] = "1"
            os.environ["LAYERNORM_TYPE"] = "torch"

            from argparse import Namespace
            from Structure.script.protenix.process_msa_to_json import main as protenix_msa_to_json
            from Structure.script.protenix.protenix_confidence import process_protenix_results

            protenix_msa_to_json(
                Namespace(
                    data=args.data_config,
                    protenix=args.model_config,
                    save_path=result_root,
                    name=target_name,
                )
            )

            if protenix_root not in sys.path:
                sys.path.insert(0, protenix_root)
            from runner.inference import run as protenix_run

            with open(args.model_config) as f:
                protenix_yaml = yaml.safe_load(f)
            inference_argv = [
                "protenix_inference",
                "--model_name", protenix_yaml["model_name"],
                "--seeds", seed,
                "--dump_dir", result_root,
                "--input_json_path", f"{result_root}/input.json",
                "--model.N_cycle", str(protenix_yaml["N_cycle"]),
                "--sample_diffusion.N_sample", str(protenix_yaml["N_sample"]),
                "--sample_diffusion.N_step", str(protenix_yaml["N_step"]),
                "--triangle_attention", "triattention",
                "--triangle_multiplicative", "cuequivariance",
                "--use_rna_msa", "true",
                "--use_template", "true",
            ]

            min_size_test = protenix_yaml.get("data.msa.min_size.test")
            if min_size_test is not None:
                inference_argv += ["--data.msa.min_size.test", str(min_size_test)]

            # TFG's VinaStericPotential crashes on single-chain inputs:
            # potentials.py:1206 calls a closure with 2 positional args that
            # is defined to take 1. Skip TFG for monomers until upstream fixes.
            total_chains = sum(e.get("copy", 1) for e in data_yaml.get("a3m") or [])
            total_chains += sum(l.get("copy", 1) for l in data_yaml.get("ligand") or [])
            is_monomer = total_chains == 1

            if protenix_yaml["use_tfg_guidance"]:
                if is_monomer:
                    print("Skipping TFG: monomer input triggers protenix VinaSteric bug.")
                else:
                    inference_argv += ["--sample_diffusion.guidance.enable", "true"]

            old_argv = sys.argv
            sys.argv = inference_argv
            try:
                protenix_run()
            finally:
                sys.argv = old_argv

            # Confidence scoring
            protenix_output = f"{result_root}/{target_name}"
            process_protenix_results(protenix_output, job_name)

            # copy to common
            for file in glob.glob(f"{protenix_output}/seed_*/predictions/*.pdb"):
                os.system(f"cp {file} {result_root}/common/")
            for file in glob.glob(f"{protenix_output}/*.png"):
                os.system(f"mv {file} {result_root}/common/")
            for file in glob.glob(f"{protenix_output}/*.csv"):
                os.system(f"mv {file} {result_root}/common/")

    # Write method log (inherit from MSA)
    method_log_path = data_yaml.get("method_log")
    if method_log_path and os.path.exists(method_log_path):
        with open(method_log_path) as f:
            method_log = yaml.safe_load(f)
    else:
        method_log = {"msa": None}
    method_log["structure"] = args.model
    with open(os.path.join(result_root, "common", "method_log.yaml"), "w") as f:
        yaml.dump(method_log, f)

    return result_root


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model",
        type=str,
        required=True,
        choices=["boltz2", "chai1", "protenix", "esmfold2"],
        help="The model to use for inference.",
    )
    parser.add_argument(
        "--data_config",
        type=str,
        required=True,
        help="Path to the data configuration yaml file.",
    )
    parser.add_argument(
        "--model_config",
        type=str,
        required=True,
        help="Path to the model configuration yaml file.",
    )

    args = parser.parse_args()
    structure_prediction(args)
