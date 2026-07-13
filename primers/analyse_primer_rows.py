from utils.gc_content import gc_content
from perturbation.mismatch_counter import best_match
from thermodynamics.primer_access import binding_accessibility
from thermodynamics.simple_thermo import thermo

def analyse_row(row, template_seq, template_structure, template_mfe):
    primer_name = row["Primer"]
    primer_seq = row["Sequence"]

    mm, pos = best_match(primer_seq, template_seq)

    if pos != -1:
        acc = binding_accessibility(template_structure, pos, pos + len(primer_seq))
    else:
        acc = 0.0

    hp, hd = thermo(primer_seq)

    return {
        "Primer": primer_name,
        "Sequence": primer_seq,
        "GC%": gc_content(primer_seq),
        "Mismatches": mm,
        "Binding_position": pos,
        "Accessibility_%": acc,
        "Hairpin_dG": hp,
        "Homodimer_dG": hd,
        "Transcript_MFE": template_mfe
    }
