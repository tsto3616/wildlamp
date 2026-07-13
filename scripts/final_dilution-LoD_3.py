#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
from Bio import SeqIO
import primer3.bindings as p3
import csv

############################################################
# MANUFACTURER BASELINE RGB (DEFAULT START COLOUR)
############################################################

HNB_BASELINE = (30, 60, 180)     # Blue
PR_BASELINE  = (220, 80, 120)    # Pink

############################################################
# USER END RGB (MEASURED COLOUR)
############################################################

USER_END_RGB  = (240, 100, 140)   # (R,G,B)
USER_END_TIME = 60   # minutes

############################################################
# GENE CONFIGURATION
############################################################

GENES = {
    "DIO2": {
        "fasta": "DIO2_N_fur.fasta",
        "F3": "AGGGGCTGACTTCCTGTT",
        "F2": "AGGGGCTGACTTCCTGTT",
        "F1": "CAGGAAGACGGATGTGCA",
        "B3": "GCCTTTGAACCTGTGTGT",
        "B2": "GCCCCTTCTCCTACAACC",
        "B1": "ACGTTGGCTGGAGAAGAA"
    },
    "DIO1": {
        "fasta": "DIO1_N_fur.fasta",
        "F3": "TTCTGGAAAAGCTCCACA",
        "F2": "GTAGTCTGGAAAGATTCC",
        "F1": "CTTGACTCTCCTCCCAGC",
        "B3": "GACTCCCGTCCCCTGCTGA",
        "B2": "AGCTCGGATTGGTCTGCT",
        "B1": "AAAGCCGCAAGTTAGCTGA"
    },
    "DIO3": {
        "fasta": "DIO3_N_fur.fasta",
        "F3": "ACTCTCCCTCCCCTCCCT",
        "F2": "GCTGGCTCGACATCCCTT",
        "F1": "GAGTGTGGATGTGATGCC",
        "B3": "CCAGACGCAGAGAGACTT",
        "B2": "GCCTGCCGCTCGCCCTAGTT",
        "B1": "TTGTTCCTGAGTTGGGGGTG"
    },
    "SLC25A43": {
        "fasta": "slc25a43_N_fur.fasta",
        "F3": "GGCTGGGAGTCTTGCAGG",
        "F2": "TGTGCAGAACCTGCTGGA",
        "F1": "CCTCCATGCTTTTTCTAC",
        "B3": "GAAACTCTGGAATGGACG",
        "B2": "CCTCCAGAACTTCGCCAA",
        "B1": "CGTAACCCAGACGCTCTC"
    },
    "HSP70": {
        "fasta": "HSP70_N_fur.fasta",
        "F3": "TCTGCTTGGATTTCCGTG",
        "F2": "ATGTCACAATAAAGCCATC",
        "F1": "AAGACCTTTGCAGACACT",
        "B3": "GTAGCTTGCATTTGAACT",
        "B2": "CCCCATATTAGTTTTGAATG",
        "B1": "ACTATCATTTTACAGTAGTCACACTTGG"
    }
}

############################################################
# BASIC UTILITIES
############################################################

def rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def build_FIP(F1: str, F2: str) -> str:
    return rc(F1) + F2

def build_BIP(B1: str, B2: str) -> str:
    return rc(B1) + B2

def load_gene(path: str) -> str:
    rec = next(SeqIO.parse(path, "fasta"))
    return str(rec.seq).upper()

############################################################
# MISMATCH / BINDING APPROXIMATION
############################################################

def best_match_mismatches(primer: str, template: str) -> int:
    Lp = len(primer)
    Lt = len(template)
    min_mm = Lp
    for i in range(Lt - Lp + 1):
        window = template[i:i+Lp]
        mm = sum(1 for a, b in zip(primer, window) if a != b)
        if mm < min_mm:
            min_mm = mm
            if min_mm == 0:
                break
    return min_mm

############################################################
# THERMODYNAMICS (LAMP-CORRECTED)
############################################################

def compute_hairpin(seq: str) -> float:
    return p3.calc_hairpin(seq).dg

def compute_homodimer(seq: str) -> float:
    return p3.calc_homodimer(seq).dg

def compute_heterodimer(seqA: str, seqB: str) -> float:
    return p3.calc_heterodimer(seqA, seqB).dg

def has_3prime_dimer(seqA: str, seqB: str, min_match: int = 4) -> bool:
    rcB = rc(seqB)
    tailA = seqA[-8:]
    tailB = rcB[-8:]
    matches = sum(1 for a, b in zip(tailA, tailB) if a == b)
    return matches >= min_match

def compute_primer_features(seq: str,
                            template: str,
                            all_primers: dict) -> dict:
    tm = p3.calc_tm(seq)
    hairpin_dG = compute_hairpin(seq)
    homodimer_dG = compute_homodimer(seq)
    mismatches = best_match_mismatches(seq, template)

    hetero = {}
    dangerous_3prime = {}
    for nameB, seqB in all_primers.items():
        if seqB == seq:
            continue
        dG = compute_heterodimer(seq, seqB)
        hetero[nameB] = dG
        dangerous_3prime[nameB] = has_3prime_dimer(seq, seqB)

    return {
        "tm": tm,
        "hairpin_dG": hairpin_dG,
        "homodimer_dG": homodimer_dG,
        "heterodimer_dG": hetero,
        "dangerous_3prime": dangerous_3prime,
        "mismatches": mismatches
    }

############################################################
# BIOPHYSICAL LAMP PRIMER SCORING
############################################################

def fraction_sequestered(dG_struct: float, scale: float = 4.0) -> float:
    dG = dG_struct / 100.0
    return 1.0 / (1.0 + math.exp(dG / scale))

def fraction_free(dG_struct: float, scale: float = 4.0) -> float:
    return 1.0 - fraction_sequestered(dG_struct, scale=scale)

def binding_efficiency(mismatches: int, per_mismatch_drop: float = 0.25) -> float:
    return max(0.0, (1.0 - per_mismatch_drop) ** mismatches)

def dimer_sink_factor(dangerous_3prime_flags: dict, per_dimer_drop: float = 0.5) -> float:
    n = sum(1 for v in dangerous_3prime_flags.values() if v)
    return (1.0 - per_dimer_drop) ** n

def lamp_tm_score(tm: float) -> float:
    return max(0.0, 1.0 - abs(tm - 60.0) / 12.0)

def lamp_primer_biophysical_score(feat: dict) -> float:
    f_free_struct = (
        fraction_free(feat["hairpin_dG"]) *
        fraction_free(feat["homodimer_dG"])
    )
    f_bind = binding_efficiency(feat["mismatches"])
    f_dimer = dimer_sink_factor(feat["dangerous_3prime"])
    tm_term = lamp_tm_score(feat["tm"])
    score = tm_term * f_free_struct * f_bind * f_dimer
    return max(0.0, score)

############################################################
# RT-LAMP SIMULATION
############################################################

def rt_lamp_rate(FIP_feat: dict,
                 BIP_feat: dict,
                 outer_feats: list) -> float:
    inner = lamp_primer_biophysical_score(FIP_feat) + lamp_primer_biophysical_score(BIP_feat)
    outer = sum(lamp_primer_biophysical_score(f) for f in outer_feats)
    combined = 0.7 * inner + 0.3 * outer
    base_k = 0.15
    scale_k = 0.40
    return base_k + scale_k * combined

def simulate_rt_qlamp(FIP_seq: str,
                      BIP_seq: str,
                      outer_seqs: list,
                      template_seq: str,
                      N0: float = 800,
                      t_max: float = 60,
                      dt: float = 0.25):

    all_primers = {
        "FIP": FIP_seq,
        "BIP": BIP_seq,
        "F3": outer_seqs[0],
        "F2": outer_seqs[1],
        "B3": outer_seqs[2]
    }

    FIP_feat = compute_primer_features(FIP_seq, template_seq, all_primers)
    BIP_feat = compute_primer_features(BIP_seq, template_seq, all_primers)
    outer_feats = [
        compute_primer_features(p, template_seq, all_primers)
        for p in outer_seqs
    ]

    k_eff = rt_lamp_rate(FIP_feat, BIP_feat, outer_feats)

    threshold_copies = 1e6
    plateau = 9e8
    t_rt = 8.0
    rt_eff = 0.45

    times, copies = [], []
    tt = None
    cDNA = 0.0
    t = 0.0

    while t <= t_max:
        if t < t_rt:
            cDNA = N0 * (t / t_rt) * rt_eff
            N = cDNA
        else:
            N = cDNA * math.exp(k_eff * (t - t_rt))

        if N > plateau:
            N = plateau

        times.append(t)
        copies.append(N)

        if tt is None and N >= threshold_copies:
            tt = t

        t += dt

    return np.array(times), np.array(copies), tt

############################################################
# DYE MODELS (HNB + PHENOL RED)
############################################################

def color_hnb_from_copies(N: float):
    Mg0 = 6.0
    beta = 1.1e-8
    Mg_free = max(Mg0 - beta * N, 0.0)
    x = 1 / (1 + math.exp(-(Mg_free - 3.5) * 2))
    r0, g0, b0 = HNB_BASELINE
    r1, g1, b1 = 120, 220, 220
    return (
        int(r0 + (1 - x) * (r1 - r0)),
        int(g0 + (1 - x) * (g1 - g0)),
        int(b0 + (1 - x) * (b1 - b0))
    )

def color_pr_from_copies(N: float):
    pH0 = 8.3
    alpha = 0.12
    pH = pH0 - alpha * math.log10(max(N, 1))
    x = 1 / (1 + math.exp(-(pH - 7.4) * 4))
    r0, g0, b0 = PR_BASELINE
    r1, g1, b1 = 240, 220, 60
    return (
        int(r0 + (1 - x) * (r1 - r0)),
        int(g0 + (1 - x) * (g1 - g0)),
        int(b0 + (1 - x) * (b1 - b0))
    )

############################################################
# ΔRGB + HYBRID QUANTIFICATION
############################################################

def delta_rgb(R: int, G: int, B: int,
              R0: int, G0: int, B0: int) -> int:
    return abs(R - R0) + abs(G - G0) + abs(B - B0)

QUANT_CONC = [1e6, 3e5, 1e5, 3e4, 1e4, 3e3, 1e3, 3e2, 1e2]

def simulate_quant_curve(FIP: str,
                         BIP: str,
                         outer: list,
                         template_seq: str,
                         dye_func,
                         baseline):

    logN_list = []
    tt_list = []
    tt_colour_list = []

    R0, G0, B0 = baseline

    for N0 in QUANT_CONC:
        times, copies, tt = simulate_rt_qlamp(
            FIP,
            BIP,
            outer,
            template_seq,
            N0=N0
        )

        tt_colour = None
        for t, N in zip(times, copies):
            R, G, B = dye_func(N)
            dRGB = delta_rgb(R, G, B, R0, G0, B0)
            if dRGB >= 40:
                tt_colour = t
                break

        logN_list.append(np.log10(N0))
        tt_list.append(tt)
        tt_colour_list.append(tt_colour)

    return np.array(logN_list), np.array(tt_list), np.array(tt_colour_list)

def fit_standard_curve(logN, tt):
    tt_arr = np.array([np.nan if (v is None) else v for v in tt], dtype=float)
    mask = ~np.isnan(tt_arr)
    slope, intercept = np.polyfit(logN[mask], tt_arr[mask], 1)
    return intercept, slope

def copies_from_Tt_colour(tt_colour: float,
                          a_colour: float,
                          b_colour: float) -> float:
    logN = (tt_colour - a_colour) / b_colour
    return 10 ** logN

############################################################
# PLOTTING
############################################################

GENE_COLOURS = {
    "DIO2": "blue",
    "DIO1": "red",
    "DIO3": "green",
    "SLC25A43": "purple",
    "HSP70": "orange"
}

def plot_multi_gene_curve(results: dict,
                          key: str,
                          title: str,
                          ylabel: str):
    plt.figure(figsize=(10, 7))
    for gene, data in results.items():
        logN = data["logN"]
        y = data[key]
        a, b = data[f"a_{key}"], data[f"b_{key}"]
        plt.scatter(logN, y, color=GENE_COLOURS[gene])
        xfit = np.linspace(min(logN), max(logN), 200)
        plt.plot(xfit, a + b * xfit,
                 color=GENE_COLOURS[gene],
                 linewidth=2,
                 label=gene)
    plt.xlabel("log10(Input copies)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_delta_rgb():
    copies = np.logspace(2, 8, 300)
    R0_h, G0_h, B0_h = HNB_BASELINE
    R0_p, G0_p, B0_p = PR_BASELINE
    d_hnb = [delta_rgb(*color_hnb_from_copies(N), R0_h, G0_h, B0_h)
             for N in copies]
    d_pr = [delta_rgb(*color_pr_from_copies(N), R0_p, G0_p, B0_p)
            for N in copies]
    plt.figure(figsize=(8, 6))
    plt.plot(np.log10(copies), d_hnb,
             color="cyan", linewidth=2, label="HNB")
    plt.plot(np.log10(copies), d_pr,
             color="magenta", linewidth=2, label="Phenol Red")
    plt.axhline(40, linestyle="--", color="black",
                label="ΔRGB threshold (40)")
    plt.xlabel("log10(copies)")
    plt.ylabel("ΔRGB")
    plt.title("ΔRGB vs log10(copies)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

############################################################
# THERMODYNAMICS TABLE + CSV
############################################################

def generate_all_thermo_tables(GENES: dict):
    rows = []
    for gene_name, cfg in GENES.items():
        gene_seq = load_gene(cfg["fasta"])
        primers = {
            "F3": cfg["F3"],
            "F2": cfg["F2"],
            "F1": cfg["F1"],
            "B3": cfg["B3"],
            "B2": cfg["B2"],
            "B1": cfg["B1"]
        }

        print(f"\n=== Thermodynamics Table for {gene_name} ===")
        print(f"{'Primer':<6} {'Tm':>6} {'Hairpin_dG':>12} "
              f"{'Homodimer_dG':>14} {'WorstHet_dG':>14} "
              f"{'Mismatches':>12}")
        print("-" * 70)

        for name, seq in primers.items():
            feat = compute_primer_features(seq, gene_seq, primers)
            worst_het = min(feat["heterodimer_dG"].values()) \
                if feat["heterodimer_dG"] else 0.0

            print(f"{name:<6} "
                  f"{feat['tm']:>6.2f} "
                  f"{feat['hairpin_dG']:>12.2f} "
                  f"{feat['homodimer_dG']:>14.2f} "
                  f"{worst_het:>14.2f} "
                  f"{feat['mismatches']:>12}")

            rows.append({
                "Gene": gene_name,
                "Primer": name,
                "Sequence": seq,
                "Tm": feat["tm"],
                "Hairpin_dG": feat["hairpin_dG"],
                "Homodimer_dG": feat["homodimer_dG"],
                "WorstHeterodimer_dG": worst_het,
                "Mismatches": feat["mismatches"]
            })
    return rows

def write_thermo_csv(rows, filename: str = "primer_thermodynamics.csv"):
    fieldnames = ["Gene", "Primer", "Sequence", "Tm",
                  "Hairpin_dG", "Homodimer_dG",
                  "WorstHeterodimer_dG", "Mismatches"]
    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nThermodynamics CSV written to: {filename}")

def generate_standard_curve_summary(results: dict,
                                    filename: str = "standard_curve_summary.csv"):
    fieldnames = ["Gene", "Slope", "Intercept", "R2"]
    with open(filename, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for gene, data in results.items():
            logN = data["logN"]
            tt = data["tt"]
            tt_arr = np.array([np.nan if v is None else v for v in tt],
                              dtype=float)
            mask = ~np.isnan(tt_arr)
            slope, intercept = np.polyfit(logN[mask], tt_arr[mask], 1)
            y_pred = intercept + slope * logN[mask]
            ss_res = np.sum((tt_arr[mask] - y_pred) ** 2)
            ss_tot = np.sum((tt_arr[mask] - np.mean(tt_arr[mask])) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
            writer.writerow({
                "Gene": gene,
                "Slope": slope,
                "Intercept": intercept,
                "R2": r2
            })
    print(f"\nStandard curve summary CSV written to: {filename}")

############################################################
# MAIN
############################################################

def main():
    results_hnb = {}
    results_pr = {}

    for gene_name, cfg in GENES.items():
        gene_seq = load_gene(cfg["fasta"])
        FIP = build_FIP(cfg["F1"], cfg["F2"])
        BIP = build_BIP(cfg["B1"], cfg["B2"])
        outer = [cfg["F3"], cfg["F2"], cfg["B3"]]

        logN_hnb, tt_hnb, ttc_hnb = simulate_quant_curve(
            FIP,
            BIP,
            outer,
            gene_seq,
            color_hnb_from_copies,
            HNB_BASELINE
        )
        a_c_hnb, b_c_hnb = fit_standard_curve(logN_hnb, tt_hnb)
        a_col_hnb, b_col_hnb = fit_standard_curve(logN_hnb, ttc_hnb)

        results_hnb[gene_name] = {
            "logN": logN_hnb,
            "tt": tt_hnb,
            "tt_colour": ttc_hnb,
            "a_tt": a_c_hnb, "b_tt": b_c_hnb,
            "a_tt_colour": a_col_hnb, "b_tt_colour": b_col_hnb
        }

        logN_pr, tt_pr, ttc_pr = simulate_quant_curve(
            FIP,
            BIP,
            outer,
            gene_seq,
            color_pr_from_copies,
            PR_BASELINE
        )
        a_c_pr, b_c_pr = fit_standard_curve(logN_pr, tt_pr)
        a_col_pr, b_col_pr = fit_standard_curve(logN_pr, ttc_pr)

        results_pr[gene_name] = {
            "logN": logN_pr,
            "tt": tt_pr,
            "tt_colour": ttc_pr,
            "a_tt": a_c_pr, "b_tt": b_c_pr,
            "a_tt_colour": a_col_pr, "b_tt_colour": b_col_pr
        }

    plot_multi_gene_curve(results_hnb, "tt",
                          "Dilution Curve (HNB)", "Tt")
    plot_multi_gene_curve(results_hnb, "tt_colour",
                          "Colour Curve (HNB)", "Tt_colour")

    plot_multi_gene_curve(results_pr, "tt",
                          "Dilution Curve (Phenol Red)", "Tt")
    plot_multi_gene_curve(results_pr, "tt_colour",
                          "Colour Curve (Phenol Red)", "Tt_colour")

    plot_delta_rgb()

    if USER_END_TIME is not None:
        print("User endpoint time:", USER_END_TIME, "minutes")
        gene_ref = "DIO2"
        a_col = results_pr[gene_ref]["a_tt_colour"]
        b_col = results_pr[gene_ref]["b_tt_colour"]
        N_est = copies_from_Tt_colour(USER_END_TIME, a_col, b_col)
        print(f"Estimated PR copies from endpoint Tt ({gene_ref}):", N_est)

    thermo_rows = generate_all_thermo_tables(GENES)
    write_thermo_csv(thermo_rows)
    generate_standard_curve_summary(results_pr, "standard_curve_PR.csv")

if __name__ == "__main__":
    main()
