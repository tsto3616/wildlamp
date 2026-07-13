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
# Primer_hits
# -----------------------------
from .primer_hits.crossreactive import blast_sequence, biologically_meaningful

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
# Species
# -----------------------------
from .species.accession import extract_accession, extract_accessions_from_csv
from .species.ncbi_lookup import accession_to_gene_id, get_species_from_accession
from .species.table import write_species_lookup

# -----------------------------
# Main pipeline entry point
# -----------------------------
from .pipeline.init_primers import init_primers
from .pipeline.species_fetch import species_fetch
from .pipeline.primer2hits import primers2hits
from .pipeline.csv2hits import csv2hits

# -----------------------------
# Utilities
# -----------------------------
from .utils.seq import gc, bad_run, rc

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

    # primer_hits
    "blast_sequence", "biologically_meaningful",
    
    # perturbation
    "perturbation_analysis",
    "mutate_base", "perturb_primer",
    
    # debug
    "debug_no_matches",
    "suggest_degenerate", "analyze_failure_modes", "validate_primer_vienna", "validate_lamp_set_vienna",
    
    # thermodynamics
    "evaluate_primer3", "vienna_fold_dna",
    
    # species
    "extract_accession", "extract_accessions_from_csv",
    "accession_to_gene_id", "get_species_from_accession",
    "write_species_lookup",
    
    # utilities
    "gc", "bad_run", "rc",
    
    # Pipeline
    "init_primers",
    "species_fetch",
    "primers2hits",
    "csv2hits",
]
