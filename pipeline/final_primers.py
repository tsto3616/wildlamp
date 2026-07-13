import csv

from alignment.load_seq import load_alignment, clean_alignment
from utils.seq import rc
from primers.optimise_primer import optimise_primer
from thermodynamics.scoring import thermo_score
from primers.build_full_lset import build_full_lset
from thermodynamics.vienna_eval import validate_lamp_set_vienna
from perturbation.analysis import perturbation_analysis
from perturbation.robustness import compute_perturbation_robustness_scores
from primers.primer_from_csv import load_primers_from_csv


def final_primers(alignment_path, primer_csv, output_csv="final_primers_output.csv"):
    """
    Run full optimisation + Vienna + perturbation pipeline and write results to CSV.
    """

    # ------------------------------------------------------------
    # LOAD ALIGNMENT
    # ------------------------------------------------------------
    names, raw = load_alignment(alignment_path)
    aln = clean_alignment(raw)
    ref = aln[0]

    # ------------------------------------------------------------
    # LOAD PRIMERS FROM CSV
    # ------------------------------------------------------------
    primers = load_primers_from_csv(primer_csv)

    required = ["F3","F2","F1","B3","B2","B1"]
    for r in required:
        if r not in primers:
            raise ValueError(f"Primer {r} missing from CSV.")

    optimised = {}

    # ------------------------------------------------------------
    # OPTIMISE EACH PRIMER
    # ------------------------------------------------------------
    for role, seq in primers.items():
        if role not in required:
            continue

        opt_seq, opt_pos, opt_score, opt_cov = optimise_primer(
            aln, names, ref,
            seq,
            window_shift=10,
            tm_min=55.0,
            tm_max=62.0,
            max_mutations=2,
            no_mutate_3prime=5,
            cost_max=4.0,
            min_binding_transcripts=5
        )

        if opt_seq is None:
            continue

        tm_final, _ = thermo_score(opt_seq)

        optimised[role] = {
            "seq": opt_seq,
            "pos": opt_pos,
            "tm": tm_final,
            "score": opt_score,
            "species_matched": opt_cov
        }

    # ------------------------------------------------------------
    # BUILD FIP / BIP (CSV overrides optimisation)
    # ------------------------------------------------------------
    if "FIP" in primers:
        FIP = primers["FIP"]
    elif "F1" in optimised and "F2" in optimised:
        FIP = rc(optimised["F1"]["seq"]) + optimised["F2"]["seq"]
    else:
        FIP = None

    if "BIP" in primers:
        BIP = primers["BIP"]
    elif "B1" in optimised and "B2" in optimised:
        BIP = rc(optimised["B1"]["seq"]) + optimised["B2"]["seq"]
    else:
        BIP = None

    # ------------------------------------------------------------
    # BUILD FULL LAMP SET FOR VIENNA + PERTURBATION
    # ------------------------------------------------------------
    lset = build_full_lset(primers, optimised, ref)

    # ------------------------------------------------------------
    # VIENNARNA FOLDING
    # ------------------------------------------------------------
    lset = validate_lamp_set_vienna(lset)

    # ------------------------------------------------------------
    # PERTURBATION ROBUSTNESS
    # ------------------------------------------------------------
    pert = perturbation_analysis(lset, aln, names)
    robust = compute_perturbation_robustness_scores(pert)

    # ------------------------------------------------------------
    # WRITE CSV OUTPUT
    # ------------------------------------------------------------
    fieldnames = [
        "role", "original_seq", "optimised_seq", "position",
        "tm", "score", "species_matched",
        "mfe", "structure", "vienna_pass",
        "FIP", "BIP",
        "robust_mean", "robust_min", "robust_max"
    ]

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for role in required:
            orig = primers[role]
            opt = optimised.get(role, {})
            vienna = lset["vienna"][role]

            writer.writerow({
                "role": role,
                "original_seq": orig,
                "optimised_seq": opt.get("seq"),
                "position": opt.get("pos"),
                "tm": opt.get("tm"),
                "score": opt.get("score"),
                "species_matched": opt.get("species_matched"),
                "mfe": vienna["mfe"],
                "structure": vienna["structure"],
                "vienna_pass": vienna["pass"],
                "FIP": FIP,
                "BIP": BIP,
                "robust_mean": robust["mean"],
                "robust_min": robust["min"],
                "robust_max": robust["max"]
            })

    # ------------------------------------------------------------
    # RETURN RESULTS
    # ------------------------------------------------------------
    return {
        "optimised": optimised,
        "FIP": FIP,
        "BIP": BIP,
        "vienna": lset["vienna"],
        "robustness": robust,
        "csv": output_csv
    }
