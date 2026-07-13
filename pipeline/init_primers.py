#!/usr/bin/env python3

import csv

from alignment.load_seq import load_alignment, clean_alignment
from alignment.conservation import compute_conservation

from primers.scan import scan_primers
from primers.species import filter_primers_by_species, primer_species_coverage
from primers.clash import build_clash_map
from primers.detect_lamp import detect_lamp_sets

from primers.amplicons import amplicon_heuristic, amplicon_score
from thermodynamics.primer3_eval import evaluate_primer3
from thermodynamics.vienna_eval import validate_lamp_set_vienna

from primers.lamp_score import lamp_score, rtlamp_score

from perturbation.analysis import perturbation_analysis
from perturbation.robustness import (
    build_robustness_matrix,
    compute_perturbation_robustness_scores
)

from failures.debug import debug_no_matches
from failures.failure_mode import analyze_failure_modes


def init_primers(alignment_path, output):
    """
    Standalone pipeline command.
    Users call:
        init_primers("alignment.fasta", "output_prefix")

    Produces:
        <prefix>_validated_top20.csv
        <prefix>_failure_modes.csv
        <prefix>_perturbation.csv
    """

    # ------------------------------------------------------------
    # LOAD + CLEAN ALIGNMENT
    # ------------------------------------------------------------
    names, raw = load_alignment(alignment_path)
    aln = clean_alignment(raw)
    ref = aln[0]

    # ------------------------------------------------------------
    # CONSERVATION + PRIMER SCANNING
    # ------------------------------------------------------------
    cons_array = compute_conservation(aln)

    primers = scan_primers(ref)
    primers = filter_primers_by_species(primers, aln, names)

    clash = build_clash_map(primers)
    lamp_sets = detect_lamp_sets(ref, primers, aln, names, clash)

    scored = []
    N_species = len(names)

    # ------------------------------------------------------------
    # SCORING LOOP
    # ------------------------------------------------------------
    for s in lamp_sets:

        species = primer_species_coverage(s, aln, names)
        s["species_matched"] = len(species)
        s["species_list"] = species

        amp_start = s["F2_pos"]
        amp_end   = s["B2_pos"] + len(s["B2"])
        amplicon  = ref[amp_start:amp_end]

        amp_info = amplicon_heuristic(amplicon)
        s["amp_gc"]      = amp_info["gc"]
        s["amp_penalty"] = amp_info["penalty"]
        s["amp_score"]   = amplicon_score(amp_info)

        s["primer3_eval"] = evaluate_primer3(s)
        s = validate_lamp_set_vienna(s)

        s["final_score"]  = lamp_score(s, cons_array, N_species)
        s["rtlamp_score"] = rtlamp_score(s)

        scored.append(s)

    print("Primers scanned:", len(primers))
    print("LAMP sets detected:", len(lamp_sets))
    print("Scored sets:", len(scored))
    print("Vienna pass:", sum(1 for s in scored if s["vienna_pass"]))

    # ------------------------------------------------------------
    # ROBUSTNESS
    # ------------------------------------------------------------
    mat, labels = build_robustness_matrix(scored, aln, names)
    compute_perturbation_robustness_scores(mat, scored)

    for s in scored:
        s["combined_score"] = 0.7 * s["rtlamp_score"] + 0.3 * s["perturb_score"]

    validated = [s for s in scored if s["vienna_pass"]]

    if len(validated) == 0:
        debug_no_matches(scored, mat, labels)
        print("No validated sets found.")
        return

    validated.sort(key=lambda x: x["combined_score"], reverse=True)

    # ------------------------------------------------------------
    # OUTPUT 1: TOP 20 VALIDATED SETS
    # ------------------------------------------------------------
    validated_file = f"{output}_validated_top20.csv"

    fieldnames = [
        "combined_score","rtlamp_score","perturb_score",
        "perturb_mean","perturb_min","perturb_frac",
        "final_score","species_matched",
        "amplicon_len","amp_gc","amp_penalty","amp_score",
        "F3","F3_pos","F2","F2_pos","F1","F1_pos",
        "B3","B3_pos","B2","B2_pos","B1","B1_pos",
        "FIP","BIP","species_list",
        "F3_tm","F3_hairpin","F3_dimer",
        "F2_tm","F2_hairpin","F2_dimer",
        "F1_tm","F1_hairpin","F1_dimer",
        "B3_tm","B3_hairpin","B3_dimer",
        "B2_tm","B2_hairpin","B2_dimer",
        "B1_tm","B1_hairpin","B1_dimer",
        "F3_mfe","F3_struct",
        "F2_mfe","F2_struct",
        "F1_mfe","F1_struct",
        "B3_mfe","B3_struct",
        "B2_mfe","B2_struct",
        "B1_mfe","B1_struct"
    ]

    with open(validated_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for s in validated[:20]:
            row = {
                "combined_score": s["combined_score"],
                "rtlamp_score": s["rtlamp_score"],
                "perturb_score": s["perturb_score"],
                "perturb_mean": s["perturb_mean"],
                "perturb_min": s["perturb_min"],
                "perturb_frac": s["perturb_frac"],
                "final_score": s["final_score"],
                "species_matched": s["species_matched"],
                "amplicon_len": s["amplicon_len"],
                "amp_gc": s["amp_gc"],
                "amp_penalty": s["amp_penalty"],
                "amp_score": s["amp_score"],
                "F3": s["F3"], "F3_pos": s["F3_pos"],
                "F2": s["F2"], "F2_pos": s["F2_pos"],
                "F1": s["F1"], "F1_pos": s["F1_pos"],
                "B3": s["B3"], "B3_pos": s["B3_pos"],
                "B2": s["B2"], "B2_pos": s["B2_pos"],
                "B1": s["B1"], "B1_pos": s["B1_pos"],
                "FIP": s["FIP"], "BIP": s["BIP"],
                "species_list": ",".join(s["species_list"]),
                "F3_tm": s["primer3_eval"]["F3"]["tm"],
                "F3_hairpin": s["primer3_eval"]["F3"]["hairpin_tm"],
                "F3_dimer": s["primer3_eval"]["F3"]["dimer_tm"],
                "F2_tm": s["primer3_eval"]["F2"]["tm"],
                "F2_hairpin": s["primer3_eval"]["F2"]["hairpin_tm"],
                "F2_dimer": s["primer3_eval"]["F2"]["dimer_tm"],
                "F1_tm": s["primer3_eval"]["F1"]["tm"],
                "F1_hairpin": s["primer3_eval"]["F1"]["hairpin_tm"],
                "F1_dimer": s["primer3_eval"]["F1"]["dimer_tm"],
                "B3_tm": s["primer3_eval"]["B3"]["tm"],
                "B3_hairpin": s["primer3_eval"]["B3"]["hairpin_tm"],
                "B3_dimer": s["primer3_eval"]["B3"]["dimer_tm"],
                "B2_tm": s["primer3_eval"]["B2"]["tm"],
                "B2_hairpin": s["primer3_eval"]["B2"]["hairpin_tm"],
                "B2_dimer": s["primer3_eval"]["B2"]["dimer_tm"],
                "B1_tm": s["primer3_eval"]["B1"]["tm"],
                "B1_hairpin": s["primer3_eval"]["B1"]["hairpin_tm"],
                "B1_dimer": s["primer3_eval"]["B1"]["dimer_tm"],
                "F3_mfe": s["vienna_validation"]["F3"]["mfe"],
                "F3_struct": s["vienna_validation"]["F3"]["structure"],
                "F2_mfe": s["vienna_validation"]["F2"]["mfe"],
                "F2_struct": s["vienna_validation"]["F2"]["structure"],
                "F1_mfe": s["vienna_validation"]["F1"]["mfe"],
                "F1_struct": s["vienna_validation"]["F1"]["structure"],
                "B3_mfe": s["vienna_validation"]["B3"]["mfe"],
                "B3_struct": s["vienna_validation"]["B3"]["structure"],
                "B2_mfe": s["vienna_validation"]["B2"]["mfe"],
                "B2_struct": s["vienna_validation"]["B2"]["structure"],
                "B1_mfe": s["vienna_validation"]["B1"]["mfe"],
                "B1_struct": s["vienna_validation"]["B1"]["structure"],
            }
            writer.writerow(row)

    print(f"Validated sets written to {validated_file}")

    # ------------------------------------------------------------
    # OUTPUT 2: FAILURE MODES
    # ------------------------------------------------------------
    failure_file = f"{output}_failure_modes.csv"

    fail_fields = [
        "species","role","primer_pos","primer_index",
        "original_primer","ref_base","alt_bases",
        "suggested_primer","mismatch_cost"
    ]

    with open(failure_file, "w", newline="") as f2:
        writer2 = csv.DictWriter(f2, fieldnames=fail_fields)
        writer2.writeheader()

        for s in validated[:20]:
            failures = analyze_failure_modes(s, aln, names)
            for fail in failures:
                writer2.writerow(fail)

    print(f"Failure modes written to {failure_file}")

    # ------------------------------------------------------------
    # OUTPUT 3: PERTURBATION TABLE
    # ------------------------------------------------------------
    perturb_file = f"{output}_perturbation.csv"

    pert_fields = [
        "set_index","role","original_primer","perturbed_primer",
        "primer_pos","species_matched","species_failed",
        "tm","hairpin_tm","dimer_tm","mfe","structure",
        "perturb_score","perturb_mean","perturb_min","perturb_frac"
    ]

    with open(perturb_file, "w", newline="") as f4:
        writer4 = csv.DictWriter(f4, fieldnames=pert_fields)
        writer4.writeheader()

        for idx, s in enumerate(scored):
            pert = perturbation_analysis(s, aln, names)
            for row in pert:
                row["set_index"] = idx
                row["perturb_score"] = s["perturb_score"]
                row["perturb_mean"]  = s["perturb_mean"]
                row["perturb_min"]   = s["perturb_min"]
                row["perturb_frac"]  = s["perturb_frac"]
                writer4.writerow(row)

    print(f"Perturbation analysis written to {perturb_file}")
    print("Pipeline complete.")
