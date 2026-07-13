"""
wildlamp: A computational pipeline for designing conserved RT-LAMP primers
across genomically uncharacterised wildlife species.

This top-level API exposes all major functions so users can simply:
    import wildlamp as wl
and call:
    wl.init_primer(...)
    wl.load_alignment(...)
    wl.mismatch_score(...)
"""

# -----------------------------
# Alignment utilities
# -----------------------------
from .alignment.load_seq import (
    load_alignment,
    clean_alignment
)

from .alignment.conservation import compute_conservation

# -----------------------------
# Primer utilities
# -----------------------------
from .primers.scan import (
    primer_ok, 
    scan_primers
)   

from .primers.mismatch import mismatch_cost
from .primers.species import (
    primer_species_coverage, 
    filter_primers_by_species
)

from .primers.clash import (
    primers_clash,
    build_clash_map
)

from .primers.detect_lamp import detect_lamp_sets

from .primers.amplicons import amplicon_heuristic, amplicon_score

from .primers.lamp_score import (primer_basic_score, 
                                 primer_thermo_score, 
                                 primer_struct_score, 
                                 primer_conservation_score,
                                 lamp_score,
                                 rtlamp_score
)

# -----------------------------
# Perturbation
# -----------------------------
from .perturbation.analysis import perturbation_analysis
from .perturbation.mutate import mutate_base, perturb_primer

# -----------------------------
# debug
# -----------------------------
from .failures.debug import debug_no_matches
from .failures.suggest import suggest_degenerate
from .failures.failure_mode import analyze_failure_modes

# -----------------------------
# thermodynamics
# -----------------------------
from .thermodynamics.primer3_eval import evaluate_primer3
from .thermodynamics.vienna_eval import vienna_fold_dna, validate_primer_vienna, validate_lamp_set_vienna

# -----------------------------
# Main pipeline entry point
# -----------------------------
from .pipeline.init_primers import init_primers

# -----------------------------
# Utilities
# -----------------------------
from .utils.seq import rc, gc, bad_run

# -----------------------------
# Public API
# -----------------------------
__all__ = [
    # Alignment
    "load_alignment", "clean_alignment",
    "compute_conservation",

    # Primers
    "scan_primers", "primer_ok", "mismatch_cost",
    "primer_species_coverage", "filter_primers_by_species",
    "primers_clash", "build_clash_map",
    "detect_lamp_sets",
    "amplicon_heuristic", "amplicon_score",
    "primer_basic_score", "primer_thermo_score", "primer_struct_score", "primer_conservation_score",
    "lamp_score", "rtlamp_score", 

    # perturbation
    "perturbation_analysis",
    "mutate_base", "perturb_primer",
    
    # debug
    "debug_no_matches",
    "suggest_degenerate", "analyze_failure_modes", "validate_primer_vienna", "validate_lamp_set_vienna",
    
    # thermodynamics
    "evaluate_primer3", "vienna_fold_dna",
    # utilities
    "gc", "rc", "bad_run",
    
    # Pipeline
    "init_primers",
]
