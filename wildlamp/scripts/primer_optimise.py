#!/usr/bin/env python

from Bio import SeqIO
import primer3
import RNA

###############################################
# ALIGNMENT + UTILITIES
###############################################

def load_alignment(path):
    names, seqs = [], []
    for rec in SeqIO.parse(path, "fasta"):
        names.append(rec.id)
        seqs.append(str(rec.seq).upper().replace("U","T"))
    return names, seqs

def clean_alignment(aln):
    max_len = max(len(s) for s in aln)
    aln = [s.ljust(max_len, "-") for s in aln]
    cols = list(zip(*aln))
    clean_cols = [col for col in cols if any(b in "ACGT" for b in col)]
    return ["".join(col[i] for col in clean_cols) for i in range(len(aln))]

def gc(seq):
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0

def bad_run(seq, max_run=6):
    run = 1
    for a, b in zip(seq, seq[1:]):
        if a == b:
            run += 1
            if run > max_run:
                return True
        else:
            run = 1
    return False

def rc(seq):
    return seq.translate(str.maketrans("ACGT","TGCA"))[::-1]

###############################################
# COST-BASED MISMATCH + COVERAGE
###############################################

def mismatch_cost(primer, seg):
    cost = 0.0
    L = len(primer)
    for i, (a, b) in enumerate(zip(primer, seg), start=1):
        if a == b:
            continue
        if b == "-":
            c_type = 3
        elif (a, b) in [("A","G"),("G","A"),("C","T"),("T","C")]:
            c_type = 1
        else:
            c_type = 2
        w_pos = 2.0 if i > L-5 else 1.0
        cost += w_pos * c_type
    return cost

def primer_species_coverage(primer, pos, aln, names, cost_max=4.0):
    ok = []
    for i, seq in enumerate(aln):
        seg = seq[pos:pos+len(primer)]
        if len(seg) != len(primer):
            continue
        if mismatch_cost(primer, seg) <= cost_max:
            ok.append(names[i])
    return ok

###############################################
# THERMODYNAMICS
###############################################

def thermo_score(seq):
    ta = primer3.thermoanalysis.ThermoAnalysis()
    tm = ta.calc_tm(seq)
    hp = ta.calc_hairpin(seq)
    hd = ta.calc_homodimer(seq)
    hairpin_tm = hp.tm if hp.structure_found else 0.0
    dimer_tm   = hd.tm if hd.structure_found else 0.0
    g = gc(seq)
    gc_pen = 0
    if g < 0.35:
        gc_pen = (0.35 - g) * 100
    elif g > 0.70:
        gc_pen = (g - 0.70) * 100
    score = tm - hairpin_tm - dimer_tm - gc_pen
    return tm, score

###############################################
# GENERIC PRIMER FILTER
###############################################

def primer_ok(p, min_len=18, max_len=30):
    if not (min_len <= len(p) <= max_len): return False
    if any(b not in "ACGT" for b in p): return False
    if not (0.35 <= gc(p) <= 0.70): return False
    if bad_run(p): return False
    return True

###############################################
# LOCAL OPTIMISER FOR ONE PRIMER
###############################################

def optimise_primer(
    aln, names, ref,
    original_primer,
    window_shift=10,
    tm_min=55.0,
    tm_max=62.0,
    max_mutations=2,
    no_mutate_3prime=5,
    cost_max=4.0,
    min_binding_transcripts=5
):
    L = len(original_primer)
    pos0 = ref.find(original_primer)
    if pos0 == -1:
        print(f"[WARN] Primer {original_primer} not found in cleaned ref.")
        return None, None, None, None

    best_seq = None
    best_score = float("-inf")
    best_pos = None
    best_cov = 0

    for shift in range(-window_shift, window_shift+1):
        pos = pos0 + shift
        if pos < 0 or pos+L > len(ref):
            continue

        window = ref[pos:pos+L]
        if any(b not in "ACGT" for b in window):
            continue

        base = window
        if not primer_ok(base, min_len=L, max_len=L):
            continue

        cov_names = primer_species_coverage(base, pos, aln, names, cost_max=cost_max)
        if len(cov_names) < min_binding_transcripts:
            continue

        tm, thermo = thermo_score(base)
        if not (tm_min <= tm <= tm_max):
            continue

        best_local_seq = base
        best_local_score = thermo + 10.0 * len(cov_names)
        best_local_cov = len(cov_names)

        for i in range(L):
            if i >= L - no_mutate_3prime:
                continue
            for b in "ACGT":
                if b == base[i]:
                    continue
                cand = base[:i] + b + base[i+1:]
                muts = sum(1 for a, c in zip(base, cand) if a != c)
                if muts > max_mutations:
                    continue

                if not primer_ok(cand, min_len=L, max_len=L):
                    continue

                cov_names_cand = primer_species_coverage(cand, pos, aln, names, cost_max=cost_max)
                cov_count = len(cov_names_cand)
                if cov_count < min_binding_transcripts:
                    continue

                tm_c, thermo_c = thermo_score(cand)
                if not (tm_min <= tm_c <= tm_max):
                    continue

                score_c = thermo_c + 10.0 * cov_count
                if score_c > best_local_score:
                    best_local_seq = cand
                    best_local_score = score_c
                    best_local_cov = cov_count

        if best_local_score > best_score:
            best_seq = best_local_seq
            best_score = best_local_score
            best_pos = pos
            best_cov = best_local_cov

    return best_seq, best_pos, best_score, best_cov

###############################################
# VIENNARNA FOLDING
###############################################

def vienna_fold_dna(dna_seq):
    rna_seq = dna_seq.replace("T","U")
    structure, mfe = RNA.fold(rna_seq)
    return {
        "rna_seq": rna_seq,
        "structure": structure,
        "mfe": mfe
    }

def vienna_validate_set(lset, cutoff=-5.0):
    roles = ["F3","F2","F1","B3","B2","B1"]
    out = {}
    all_pass = True
    for r in roles:
        seq = lset.get(r)
        if seq is None:
            out[r] = {"mfe": None, "structure": None, "pass": False}
            all_pass = False
            continue
        res = vienna_fold_dna(seq)
        ok = res["mfe"] > cutoff
        out[r] = {
            "mfe": res["mfe"],
            "structure": res["structure"],
            "pass": ok
        }
        if not ok:
            all_pass = False
    lset["vienna"] = out
    lset["vienna_pass"] = all_pass
    return lset

###############################################
# PERTURBATION GENERATION
###############################################

def mutate_base(b):
    return [x for x in "ACGT" if x != b]

def perturb_primer(primer):
    variants = []
    for i, b in enumerate(primer):
        for mb in mutate_base(b):
            variants.append(primer[:i] + mb + primer[i+1:])
    return variants

###############################################
# PERTURBATION ROBUSTNESS
###############################################

def perturbation_analysis(lset, aln, names, cost_max=4.0):
    roles = ["F3","F2","F1","B3","B2","B1"]
    results = []

    for r in roles:
        p = lset.get(r)
        pos = lset.get(f"{r}_pos")
        if p is None or pos is None:
            continue

        variants = perturb_primer(p)

        for var in variants:
            matched = []
            failed = []

            for i, seq in enumerate(aln):
                seg = seq[pos:pos+len(var)]
                if len(seg) != len(var):
                    failed.append(names[i])
                    continue
                cost = mismatch_cost(var, seg)
                if cost <= cost_max:
                    matched.append(names[i])
                else:
                    failed.append(names[i])

            ta = primer3.thermoanalysis.ThermoAnalysis()
            tm = ta.calc_tm(var)
            hp = ta.calc_hairpin(var)
            hd = ta.calc_homodimer(var)

            v = vienna_fold_dna(var)

            results.append({
                "role": r,
                "original": p,
                "variant": var,
                "pos": pos,
                "species_matched": len(matched),
                "species_failed": ",".join(failed),
                "tm": tm,
                "hairpin_tm": hp.tm if hp.structure_found else 0.0,
                "dimer_tm": hd.tm if hd.structure_found else 0.0,
                "mfe": v["mfe"],
                "structure": v["structure"]
            })

    return results

###############################################
# ROBUSTNESS SCORE
###############################################

def compute_robustness_score(pert_results):
    scores = []
    for r in pert_results:
        score = (
            1.0 * r["species_matched"] +
            0.5 * (r["tm"] - r["hairpin_tm"] - r["dimer_tm"]) +
            0.5 * (-r["mfe"])
        )
        scores.append(score)
    return {
        "mean": sum(scores)/len(scores) if scores else None,
        "min": min(scores) if scores else None,
        "max": max(scores) if scores else None
    }

###############################################
# BUILD FULL PRIMER SET (INCLUDING FIP/BIP)
###############################################

def build_full_lset(primers, optimised, ref):
    lset = {}
    for r in ["F3","F2","F1","B3","B2","B1"]:
        seq = optimised[r]["seq"] if r in optimised else primers[r]
        pos = optimised[r]["pos"] if r in optimised else ref.find(primers[r])
        lset[r] = seq
        lset[f"{r}_pos"] = pos

    lset["FIP"] = rc(lset["F1"]) + lset["F2"]
    lset["BIP"] = rc(lset["B1"]) + lset["B2"]

    return lset

###############################################
# MAIN
###############################################

def main():
    names, raw = load_alignment("DIO3_pinnipedia.fas")
    aln = clean_alignment(raw)
    ref = aln[0]

    # USER: specify all 6 primers
    F3 = "TCTGCTTGGATTTCCGTG"
    F2 = "ATGTCACAATAAAGCCATC"
    F1 = "AAGACCTTTGCAGACACT"
    B3 = "GTAGCTTGCATTTGAACT"
    B2 = "CCCCATATTAGTTTTGAATG"
    B1 = "ACTATCATTTTACAGTAGTCACACTTGG"

    primers = {
        "F3": F3,
        "F2": F2,
        "F1": F1,
        "B3": B3,
        "B2": B2,
        "B1": B1
    }

    optimised = {}

    for role, seq in primers.items():
        print(f"\n=== Optimising {role} ===")
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
            print(f"{role}: no candidate met binding + thermo constraints.")
            continue

        tm_final, _ = thermo_score(opt_seq)
        print(f"Original {role}: {seq}")
        print(f"Optimised {role}: {opt_seq}")
        print(f"{role} position (cleaned ref): {opt_pos}")
        print(f"{role} final Tm: {tm_final:.2f}")
        print(f"{role} score: {opt_score:.2f}")
        print(f"{role} species matched: {opt_cov}")

        optimised[role] = {
            "seq": opt_seq,
            "pos": opt_pos,
            "tm": tm_final,
            "score": opt_score,
            "species_matched": opt_cov
        }

    if "F1" in optimised and "F2" in optimised:
        FIP = rc(optimised["F1"]["seq"]) + optimised["F2"]["seq"]
        print("\n=== FIP (from optimised F1/F2) ===")
        print("FIP:", FIP)

    if "B1" in optimised and "B2" in optimised:
        BIP = rc(optimised["B1"]["seq"]) + optimised["B2"]["seq"]
        print("\n=== BIP (from optimised B1/B2) ===")
        print("BIP:", BIP)

    print("\n=== SUMMARY: Optimised primers ===")
    for r, info in optimised.items():
        print(f"{r}: {info['seq']} (pos={info['pos']}, Tm={info['tm']:.2f}, species={info['species_matched']})")

    # Build full primer set for Vienna + perturbation
    lset = build_full_lset(primers, optimised, ref)

    # ViennaRNA folding
    lset = vienna_validate_set(lset)

    # Perturbation analysis
    pert = perturbation_analysis(lset, aln, names)
    robust = compute_robustness_score(pert)

    print("\n=== ViennaRNA Folding ===")
    for r in ["F3","F2","F1","B3","B2","B1"]:
        v = lset["vienna"][r]
        print(f"{r}: mfe={v['mfe'] if v['mfe'] is not None else 'NA'}, pass={v['pass']}, structure={v['structure']}")

    print("\n=== Perturbation Robustness ===")
    print("Mean robustness:", robust["mean"])
    print("Min robustness:", robust["min"])
    print("Max robustness:", robust["max"])

if __name__ == "__main__":
    main()

