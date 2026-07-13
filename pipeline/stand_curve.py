import csv
import matplotlib.pyplot as plt
import numpy as np 

from ..utils.load_seq import load_fasta
from ..utils.build_primers import build_BIP, build_FIP
from ..simulations.sim_curve import simulate_quant_curve_pr
from ..simulations.curve_fit import fit_standard_curve
from ..utils.gene_primer import load_gene_primer_csv

def stand_curves(primer_gene_csv, output="standard_curve.csv"):
    """
    EXACT same structure as your original function,
    but reading from CSV instead of GENES dict.
    """

    GENES = load_gene_primer_csv(primer_gene_csv)
    results = {}

    # ---------------------------------------------------------
    # LOOP OVER GENES (identical to your original)
    # ---------------------------------------------------------
    for gene_name, cfg in GENES.items():

        # Load FASTA
        template_seq = load_fasta(cfg["fasta"])

        # Build FIP/BIP
        FIP = build_FIP(cfg["F1"], cfg["F2"])
        BIP = build_BIP(cfg["B1"], cfg["B2"])

        # Outer primers (unchanged)
        outer = [cfg["F3"], cfg["F2"], cfg["B3"]]

        # Run simulation (Phenol Red only)
        logN, tt, ttc = simulate_quant_curve_pr(FIP, BIP, outer, template_seq)

        # Fit curves
        a_tt, b_tt = fit_standard_curve(logN, tt)
        a_tc, b_tc = fit_standard_curve(logN, ttc)

        # Store results
        results[gene_name] = {
            "logN": logN,
            "tt": tt,
            "tt_colour": ttc,
            "a_tt": a_tt, "b_tt": b_tt,
            "a_tt_colour": a_tc, "b_tt_colour": b_tc
        }

    # ---------------------------------------------------------
    # PLOT Tt (unchanged)
    # ---------------------------------------------------------
    plt.figure(figsize=(10, 7))
    for gene, data in results.items():
        plt.scatter(data["logN"], data["tt"], label=f"{gene} Tt")
        xfit = np.linspace(min(data["logN"]), max(data["logN"]), 200)
        plt.plot(xfit, data["a_tt"] + data["b_tt"] * xfit)
    plt.xlabel("log10(Input copies)")
    plt.ylabel("Tt (min)")
    plt.title("Phenol Red Standard Curve (Tt)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # ---------------------------------------------------------
    # PLOT Tt_colour (unchanged)
    # ---------------------------------------------------------
    plt.figure(figsize=(10, 7))
    for gene, data in results.items():
        plt.scatter(data["logN"], data["tt_colour"], label=f"{gene} Tt_colour")
        xfit = np.linspace(min(data["logN"]), max(data["logN"]), 200)
        plt.plot(xfit, data["a_tt_colour"] + data["b_tt_colour"] * xfit)
    plt.xlabel("log10(Input copies)")
    plt.ylabel("Tt_colour (min)")
    plt.title("Phenol Red Colour Threshold Curve")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # ---------------------------------------------------------
    # WRITE CSV (unchanged)
    # ---------------------------------------------------------
    with open(output, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "Gene",
            "Slope_Tt", "Intercept_Tt",
            "Slope_Tt_colour", "Intercept_Tt_colour"
        ])
        for gene, data in results.items():
            w.writerow([
                gene,
                data["b_tt"], data["a_tt"],
                data["b_tt_colour"], data["a_tt_colour"]
            ])

    return results
