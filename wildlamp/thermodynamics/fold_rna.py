import RNA

def fold_rna(sequence):
    fc = RNA.fold_compound(sequence.replace("T", "U"))
    structure, mfe = fc.mfe()
    return structure, mfe
