"""
Diagnostic reporting for cases where no robust LAMP sets are found.

Provides:
- debug_no_matches(): print detailed diagnostics for robustness, species
  coverage, thermodynamics, ViennaRNA, and geometry failures.
"""

import numpy as np


def debug_no_matches(scored, mat, labels, roles=["F3","F2","F1","B1","B2","B3"]):
    """
    Print a detailed diagnostic report explaining why no robust LAMP sets passed.

    Parameters
    ----------
    scored : list[dict]
        List of scored LAMP sets.
    mat : np.ndarray
        Robustness matrix (perturbation robustness).
    labels : list[str]
        Column labels for robustness matrix.
    roles : list[str]
        Primer roles to inspect.

    Returns
    -------
    None
        Prints diagnostic information to stdout.
    """

    print("\n==============================")
    print(" DEBUG REPORT: No robust LAMP sets found")
    print("==============================\n")

    # 1. Global robustness
    row_means = np.nanmean(mat, axis=1)
    print(f"Global mean robustness across all sets: {np.nanmean(row_means):.3f}")
    print(f"Best set robustness: {np.nanmax(row_means):.3f}")
    print(f"Worst set robustness: {np.nanmin(row_means):.3f}\n")

    # 2. Role-specific fragility
    print("Role-specific fragility analysis:")
    for r in roles:
        role_cols = [i for i, lab in enumerate(labels) if lab.startswith(r)]
        role_vals = mat[:, role_cols]

        role_mean = np.nanmean(role_vals)
        role_min  = np.nanmin(role_vals)

        print(f"  {r}: mean={role_mean:.3f}, min={role_min:.3f}, cols={len(role_cols)}")
    print("\n")

    # 3. Species coverage failures
    print("Species coverage failures:")
    for i, s in enumerate(scored):
        if s["species_matched"] < 3:  # adjustable threshold
            print(f"  Set {i+1}: weak species coverage ({s['species_matched']})")
    print("\n")

    # 4. Thermodynamic failures
    print("Thermodynamic failures:")
    for i, s in enumerate(scored):
        for r in roles:
            ev = s["primer3_eval"][r]
            if ev["hairpin_tm"] > 40 or ev["dimer_tm"] > 40:
                print(f"  Set {i+1}, {r}: hairpin={ev['hairpin_tm']}, dimer={ev['dimer_tm']}")
    print("\n")

    # 5. ViennaRNA failures
    print("ViennaRNA failures:")
    for i, s in enumerate(scored):
        for r in roles:
            v = s["vienna_validation"][r]
            if not v["pass"]:
                print(f"  Set {i+1}, {r}: mfe={v['mfe']} structure={v['structure']}")
    print("\n")

    # 6. Geometry failures
    print("Geometry failures:")
    for i, s in enumerate(scored):
        if not s["vienna_pass"]:
            print(f"  Set {i+1}: geometry or spacing failure")

    print("\n==============================")
    print(" END DEBUG REPORT")
    print("==============================\n")
