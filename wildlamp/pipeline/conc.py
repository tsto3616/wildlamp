import csv
from simulations.dye import delta_rgb, PR_BASELINE
from simulations.curve_fit import copies_from_Tt


def estimate_copies(user_R, user_G, user_B,
                    gene_name, curve_csv,
                    output="final_csv.csv",
                    dilution_factor=1):
    """
    Estimate copies using:
    - gene name (for selecting curve parameters)
    - user-specified Phenol Red colour (R,G,B)
    - slope/intercept from standard curve CSV
    - dilution factor (default = 1, no dilution)
    - write result to output CSV
    """

    # ------------------------------------------------------------
    # Load curve parameters
    # ------------------------------------------------------------
    with open(curve_csv) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    row = next((r for r in rows if r["Gene"] == gene_name), None)
    if row is None:
        raise ValueError(f"Gene {gene_name} not found in {curve_csv}")

    intercept = float(row["Intercept_Tt_colour"])
    slope     = float(row["Slope_Tt_colour"])

    # ------------------------------------------------------------
    # Compute ΔRGB from user colour
    # ------------------------------------------------------------
    R0, G0, B0 = PR_BASELINE
    dRGB_user = delta_rgb(user_R, user_G, user_B, R0, G0, B0)

    # ------------------------------------------------------------
    # Convert ΔRGB → equivalent Tt_colour
    # ------------------------------------------------------------
    Tt_colour_user = intercept + (dRGB_user - 40)

    # ------------------------------------------------------------
    # Convert Tt_colour → copies
    # ------------------------------------------------------------
    N_est_raw = copies_from_Tt(Tt_colour_user, intercept, slope)

    # ------------------------------------------------------------
    # Apply dilution factor
    # ------------------------------------------------------------
    N_est = N_est_raw * dilution_factor

    # ------------------------------------------------------------
    # Write output CSV
    # ------------------------------------------------------------
    with open(output, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "Gene", "User_R", "User_G", "User_B",
            "DeltaRGB", "Tt_colour_est",
            "Copies_est_raw", "Dilution_factor", "Copies_est_corrected"
        ])
        w.writerow([
            gene_name,
            user_R, user_G, user_B,
            dRGB_user,
            Tt_colour_user,
            N_est_raw,
            dilution_factor,
            N_est
        ])

    return N_est
