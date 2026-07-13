#!/usr/bin/env python

import math
import csv
from Bio import SeqIO
import primer3
import RNA
from collections import Counter

###############################################
# ALIGNMENT
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

def compute_conservation(aln, ref):
    N = len(aln)
    L = len(ref)
    cons = []
    for j in range(L):
        r = ref[j]
        if r not in "ACGT":
            cons.append(0.0)
            continue
        matches = sum(1 for seq in aln if seq[j] == r)
        cons.append(matches / N)
    return cons

###############################################
# BASIC UTILITIES
###############################################

def gc(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)

def bad_run(seq, max_run=6):
    run = 1
    for a,b in zip(seq, seq[1:]):
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
# MISMATCH COST FUNCTION
###############################################

def mismatch_cost(primer, seg):
    cost = 0.0
    L = len(primer)
    for i,(a,b) in enumerate(zip(primer,seg), start=1):
        if a == b:
            continue
        if b == "-":
            c_type = 3
        elif (a,b) in [("A","G"),("G","A"),("C","T"),("T","C")]:
            c_type = 1
        else:
            c_type = 2
        w_pos = 2.0 if i > L-5 else 1.0
        cost += w_pos * c_type
    return cost

###############################################
# PRIMER FILTERS
###############################################

def primer_ok(p):
    if not (18 <= len(p) <= 35): return False
    if any(b not in "ACGT" for b in p): return False
    if not (0.35 <= gc(p) <= 0.70): return False
    if bad_run(p): return False
    return True

def scan_primers(ref):
    primers=[]
    pos=0
    while pos < len(ref):
        found=False
        for k in range(18,36):
            if pos+k > len(ref): break
            p = ref[pos:pos+k]
            if primer_ok(p):
                primers.append((pos,p))
                pos += max(k,15)
                found=True
                break
        if not found:
            pos += 1
    return primers

###############################################
# SPECIES COVERAGE (COST-BASED)
###############################################

def primer_species_coverage(primer,pos,aln,names,cost_max=4.0):
    ok=[]
    for i,seq in enumerate(aln):
        seg = seq[pos:pos+len(primer)]
        if len(seg)!=len(primer): continue
        if mismatch_cost(primer,seg) <= cost_max:
            ok.append(names[i])
    return ok

def filter_primers_by_species(primers,aln,names,cost_max=4.0,min_species=5):
    return [(pos,p) for pos,p in primers
            if len(primer_species_coverage(p,pos,aln,names,cost_max))>=min_species]

def species_coverage(lset,aln,names,cost_max=4.0):
    ok=[]
    roles=["F3","F2","F1","B1","B2","B3"]
    for i,seq in enumerate(aln):
        good=True
        for r in roles:
            p=lset[r]
            pos=lset[f"{r}_pos"]
            seg=seq[pos:pos+len(p)]
            if len(seg)!=len(p): good=False; break
            if mismatch_cost(p,seg) > cost_max:
                good=False; break
        if good: ok.append(names[i])
    return ok

###############################################
# CLASH MAPPING
###############################################

def primers_clash(pos1,len1,pos2,len2):
    s1,e1=pos1,pos1+len1
    s2,e2=pos2,pos2+len2
    if e1<=s2: return (s2-e1)<15
    if e2<=s1: return (s1-e2)<15
    return True

def build_clash_map(primers):
    clash={}
    for i,(p1_pos,p1_seq) in enumerate(primers):
        clash[i]=set()
        for j,(p2_pos,p2_seq) in enumerate(primers):
            if i==j: continue
            if primers_clash(p1_pos,len(p1_seq),p2_pos,len(p2_seq)):
                clash[i].add(j)
    return clash

###############################################
# AMPLICON HEURISTIC
###############################################

def amplicon_heuristic(seq):
    L=len(seq)
    g=gc(seq)
    length_pen=abs(L-200)*0.5 if not (120<=L<=280) else 0
    gc_pen=0
    if g<0.40: gc_pen=(0.40-g)*100
    if g>0.65: gc_pen=(g-0.65)*100
    homo_pen=30 if bad_run(seq) else 0
    return {"length":L,"gc":g,"penalty":length_pen+gc_pen+homo_pen}

###############################################
# DETECT LAMP SETS
###############################################

def detect_lamp_sets(ref,primers,aln,names,clash):
    best=[]
    best_cov=0
    L=len(ref)

    def cov(s): return len(species_coverage(s,aln,names))

    for win_start in range(0,L,100):
        win_end=win_start+500
        window=[(i,pos,p) for i,(pos,p) in enumerate(primers)
                if win_start<=pos<=win_end]
        if len(window)<6: continue

        for f3_i,f3_pos,f3 in window:
            for f2_i,f2_pos,f2 in window:
                if f2_pos<=f3_pos: continue
                if not (40<=f2_pos-f3_pos<=300): continue
                if f2_i in clash[f3_i]: continue

                for f1_i,f1_pos,f1 in window:
                    if f1_pos<=f2_pos: continue
                    if not (20<=f1_pos-f2_pos<=200): continue
                    if f1_i in clash[f2_i] or f1_i in clash[f3_i]: continue

                    F1c=rc(f1)

                    for b3_i,b3_pos,b3 in window:
                        if b3_pos<=f1_pos: continue
                        if not (100<=b3_pos-f1_pos<=600): continue
                        if b3_i in clash[f1_i] or b3_i in clash[f2_i] or b3_i in clash[f3_i]: continue

                        for b2_i,b2_pos,b2 in window:
                            if b2_pos<=b3_pos: continue
                            if not (20<=b2_pos-b3_pos<=200): continue
                            if b2_i in clash[b3_i]: continue

                            for b1_i,b1_pos,b1 in window:
                                if b1_pos<=b2_pos: continue
                                if not (20<=b1_pos-b2_pos<=200): continue
                                if b1_i in clash[b2_i] or b1_i in clash[b3_i]: continue

                                B1c=rc(b1)
                                FIP=F1c+f2
                                BIP=B1c+b2

                                s={
                                    "F3":f3,"F3_pos":f3_pos,
                                    "F2":f2,"F2_pos":f2_pos,
                                    "F1":f1,"F1_pos":f1_pos,"F1c":F1c,
                                    "B3":b3,"B3_pos":b3_pos,
                                    "B2":b2,"B2_pos":b2_pos,
                                    "B1":b1,"B1_pos":b1_pos,"B1c":B1c,
                                    "FIP":FIP,"BIP":BIP,
                                    "amplicon_len":(b2_pos+len(b2))-f2_pos
                                }

                                c=cov(s)
                                if c<best_cov: continue
                                if c>best_cov:
                                    best_cov=c
                                    best=[]
                                best.append(s)
                                if len(best)>50:
                                    best=best[-50:]
    return best

###############################################
# PRIMER3 THERMODYNAMICS
###############################################

def evaluate_primer3(s):
    ta=primer3.thermoanalysis.ThermoAnalysis()
    out={}
    for r in ["F3","F2","F1","B3","B2","B1"]:
        p=s[r]
        tm=ta.calc_tm(p)
        hp=ta.calc_hairpin(p)
        hd=ta.calc_homodimer(p)
        out[r]={
            "tm":tm,
            "hairpin_tm":hp.tm if hp.structure_found else 0.0,
            "dimer_tm":hd.tm if hd.structure_found else 0.0
        }
    return out

###############################################
# VIENNARNA STRUCTURE
###############################################

def vienna_fold_dna(dna_seq):
    rna_seq=dna_seq.replace("T","U")
    structure,mfe=RNA.fold(rna_seq)
    return {"rna_seq":rna_seq,"structure":structure,"mfe":mfe}

def validate_primer_vienna(dna_seq,cutoff=-5.0):
    res=vienna_fold_dna(dna_seq)
    return {
        "primer":dna_seq,
        "rna_seq":res["rna_seq"],
        "structure":res["structure"],
        "mfe":res["mfe"],
        "pass":res["mfe"]>cutoff
    }

def validate_lamp_set_vienna(s,cutoff=-5.0):
    roles=["F3","F2","F1","B3","B2","B1"]
    results={}
    all_pass=True
    for r in roles:
        res=validate_primer_vienna(s[r],cutoff=cutoff)
        results[r]=res
        if not res["pass"]:
            all_pass=False
    s["vienna_validation"]=results
    s["vienna_pass"]=all_pass
    return s

###############################################
# SCORING MODELS
###############################################

def primer_basic_score(p):
    score=0
    if 0.40<=gc(p)<=0.65: score+=20
    if 18<=len(p)<=30: score+=20
    if not bad_run(p): score+=10
    return score

def primer_thermo_score(ev):
    tm=ev["tm"]
    hairpin_tm=ev["hairpin_tm"]
    dimer_tm=ev["dimer_tm"]
    score=0
    if 60<=tm<=68: score+=20
    elif 55<=tm<=72: score+=10
    if hairpin_tm>40: score-=10
    if dimer_tm>40: score-=10
    return score

def primer_struct_score(vres):
    return 10 if vres["mfe"]>-5.0 else 0

def primer_conservation_score(pos,length,cons_array):
    c=sum(cons_array[pos+k] for k in range(length))/length
    return 20*c

def amplicon_score(amp_info):
    L=amp_info["length"]
    g=amp_info["gc"]
    pen=amp_info["penalty"]
    score=0
    if 140<=L<=240: score+=40
    elif 100<=L<=280: score+=20
    if 0.40<=g<=0.65: score+=20
    score -= 0.5*pen
    return score

def lamp_score(s,cons_array,N_species):
    roles=["F3","F2","F1","B1","B2","B3"]
    primer_score=0
    for r in roles:
        p=s[r]
        pos=s[f"{r}_pos"]
        basic=primer_basic_score(p)
        thermo=primer_thermo_score(s["primer3_eval"][r])
        struct=primer_struct_score(s["vienna_validation"][r])
        cons=primer_conservation_score(pos,len(p),cons_array)
        primer_score += basic + thermo + struct + cons
    cov_score = 100.0*(s["species_matched"]/N_species)
    amp_score_val = s["amp_score"]
    return primer_score + cov_score + amp_score_val

def rtlamp_score(s):
    score=s["final_score"]
    L=s["amplicon_len"]
    if 140<=L<=240: score+=80
    elif 100<=L<=280: score+=40
    score -= s["amp_penalty"]*0.5
    return score

###############################################
# DEGENERATE BASE SUGGESTIONS
###############################################

DEGEN_MAP = {
    frozenset(["A","G"]): "R",
    frozenset(["C","T"]): "Y",
    frozenset(["G","C"]): "S",
    frozenset(["A","T"]): "W",
    frozenset(["G","T"]): "K",
    frozenset(["A","C"]): "M",
    frozenset(["A","C","G"]): "V",
    frozenset(["A","C","T"]): "H",
    frozenset(["A","G","T"]): "D",
    frozenset(["C","G","T"]): "B",
    frozenset(["A","C","G","T"]): "N",
}

def suggest_degenerate(primer, idx, bases):
    base_set = frozenset(bases)
    if base_set in DEGEN_MAP:
        deg = DEGEN_MAP[base_set]
        return primer[:idx] + deg + primer[idx+1:]
    return primer

###############################################
# FAILURE MODE ANALYSIS
###############################################

def analyze_failure_modes(lset, aln, names, cost_max=4.0):
    failures = []
    roles = ["F3","F2","F1","B1","B2","B3"]

    for i, seq in enumerate(aln):
        species_name = names[i]

        for r in roles:
            p = lset[r]
            pos = lset[f"{r}_pos"]
            seg = seq[pos:pos+len(p)]

            if len(seg) != len(p):
                continue

            cost = mismatch_cost(p, seg)
            if cost <= cost_max:
                continue

            mismatch_positions = [
                idx for idx,(a,b) in enumerate(zip(p,seg))
                if a != b and b in "ACGT"
            ]

            for mp in mismatch_positions:
                genomic_pos = pos + mp

                bases = []
                for sseq in aln:
                    if genomic_pos < len(sseq):
                        b = sseq[genomic_pos]
                        if b in "ACGT":
                            bases.append(b)

                if not bases:
                    continue

                alt_bases = sorted(set(bases))
                suggested = suggest_degenerate(p, mp, alt_bases)

                failures.append({
                    "species": species_name,
                    "role": r,
                    "primer_pos": pos,
                    "primer_index": mp,
                    "original_primer": p,
                    "ref_base": p[mp],
                    "alt_bases": "".join(alt_bases),
                    "suggested_primer": suggested,
                    "mismatch_cost": cost
                })

    return failures

###############################################
# PERTURBATION ANALYSIS
###############################################

def mutate_base(b):
    """Return all possible single-base mutations of a nucleotide."""
    return [x for x in "ACGT" if x != b]

def perturb_primer(primer):
    """Generate all single-base perturbations of a primer."""
    variants = []
    for i, b in enumerate(primer):
        for mb in mutate_base(b):
            variants.append(primer[:i] + mb + primer[i+1:])
    return variants

def perturbation_analysis(lset, aln, names, cost_max=4.0):
    roles = ["F3","F2","F1","B1","B2","B3"]
    results = []

    for r in roles:
        p = lset[r]
        pos = lset[f"{r}_pos"]

        variants = perturb_primer(p)

        for var in variants:
            # Species coverage under perturbation
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

            # Thermodynamics
            ta = primer3.thermoanalysis.ThermoAnalysis()
            tm = ta.calc_tm(var)
            hp = ta.calc_hairpin(var)
            hd = ta.calc_homodimer(var)

            # ViennaRNA
            v = vienna_fold_dna(var)

            results.append({
                "role": r,
                "original_primer": p,
                "perturbed_primer": var,
                "primer_pos": pos,
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
# ROBUSTNESS HEATMAP (Z-SCORES, COMBINED)
###############################################

import numpy as np

def build_robustness_matrix(lamp_sets, aln, names):
    roles = ["F3","F2","F1","B1","B2","B3"]
    matrix = []
    col_labels = []

    for idx, s in enumerate(lamp_sets):
        row_scores = []

        for r in roles:
            p = s[r]
            pos = s[f"{r}_pos"]
            variants = perturb_primer(p)

            for var in variants:
                # Species coverage
                matched = 0
                for i, seq in enumerate(aln):
                    seg = seq[pos:pos+len(var)]
                    if len(seg) != len(var):
                        continue
                    cost = mismatch_cost(var, seg)
                    if cost <= 4.0:
                        matched += 1

                # Thermodynamics
                ta = primer3.thermoanalysis.ThermoAnalysis()
                tm = ta.calc_tm(var)
                hp = ta.calc_hairpin(var)
                hd = ta.calc_homodimer(var)
                hairpin_tm = hp.tm if hp.structure_found else 0.0
                dimer_tm = hd.tm if hd.structure_found else 0.0

                # ViennaRNA
                v = vienna_fold_dna(var)
                mfe = v["mfe"]

                # Robustness score
                R = (1.0 * matched) + (0.5 * (tm - hairpin_tm - dimer_tm)) + (0.5 * (-mfe))

                row_scores.append(R)
                col_labels.append(f"{r}_{var}")

        matrix.append(row_scores)

    # Pad rows to equal length
    max_len = max(len(row) for row in matrix)
    mat_fixed = np.array([row + [np.nan]*(max_len - len(row)) for row in matrix])

    return mat_fixed, col_labels

def compute_perturbation_robustness_scores(mat, scored, threshold=10.0):
    import numpy as np

    scores = []
    for i in range(mat.shape[0]):
        row = mat[i, :]
        row = row[~np.isnan(row)]

        if len(row) == 0:
            scores.append({
                "index": i,
                "mean": np.nan,
                "min": np.nan,
                "frac_above": np.nan,
                "robust_score": np.nan
            })
            continue

        mean_val = np.mean(row)
        min_val  = np.min(row)
        frac_above = np.sum(row > threshold) / len(row)

        robust_score = (1.0 * mean_val) + (0.5 * min_val) + (0.5 * frac_above)

        scores.append({
            "index": i,
            "mean": mean_val,
            "min": min_val,
            "frac_above": frac_above,
            "robust_score": robust_score
        })

    # Attach to sets
    for s, sc in zip(scored, scores):
        s["perturb_mean"]   = sc["mean"]
        s["perturb_min"]    = sc["min"]
        s["perturb_frac"]   = sc["frac_above"]
        s["perturb_score"]  = sc["robust_score"]

    return scores

def debug_no_matches(scored, mat, labels, roles=["F3","F2","F1","B1","B2","B3"]):
    import numpy as np

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

def main():
    
    ###############################################
    # LOAD + CLEAN ALIGNMENT
    ###############################################
    names, raw = load_alignment("plasminogen_receptor_aligned.fasta")
    aln = clean_alignment(raw)
    ref = aln[0]

    ###############################################
    # CONSERVATION + PRIMER SCANNING
    ###############################################
    cons_array = compute_conservation(aln, ref)

    primers = scan_primers(ref)
    primers = filter_primers_by_species(primers, aln, names)

    clash = build_clash_map(primers)

    ###############################################
    # DETECT LAMP SETS
    ###############################################
    lamp_sets = detect_lamp_sets(ref, primers, aln, names, clash)

    scored = []
    N_species = len(names)

    ###############################################
    # SCORE EACH LAMP SET (BASE SCORES)
    ###############################################
    for s in lamp_sets:

        species = species_coverage(s, aln, names)
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

    ###############################################
    # PERTURBATION ROBUSTNESS (BEFORE VALIDATION)
    ###############################################
    mat, labels = build_robustness_matrix(scored, aln, names)
    compute_perturbation_robustness_scores(mat, scored)

    ###############################################
    # COMBINED SCORE
    ###############################################
    for s in scored:
        s["combined_score"] = (
            0.7 * s["rtlamp_score"] +
            0.3 * s["perturb_score"]
        )

    ###############################################
    # FILTER VALIDATED SETS
    ###############################################
    validated = [s for s in scored if s["vienna_pass"]]

    ###############################################
    # DEBUG: NO VALIDATED SETS
    ###############################################
    if len(validated) == 0:
        debug_no_matches(scored, mat, labels)
        print("No validated sets found. Exiting.")
        return

    ###############################################
    # SORT VALIDATED BY COMBINED SCORE
    ###############################################
    validated.sort(key=lambda x: x["combined_score"], reverse=True)

    ###############################################
    # MAIN CSV OUTPUT
    ###############################################
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

    with open("HSP70_rtlamp_validated_top20.csv","w",newline="") as f:
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

    ###############################################
    # FAILURE MODES CSV
    ###############################################
    fail_fields = [
        "species","role","primer_pos","primer_index",
        "original_primer","ref_base","alt_bases",
        "suggested_primer","mismatch_cost"
    ]

    with open("plas_rtlamp_failure_modes.csv","w",newline="") as f2:
        writer2 = csv.DictWriter(f2, fieldnames=fail_fields)
        writer2.writeheader()

        for s in validated[:20]:
            failures = analyze_failure_modes(s, aln, names)
            for fail in failures:
                writer2.writerow(fail)

    ###############################################
    # PERTURBATION CSV
    ###############################################
    pert_fields = [
        "set_index","role","original_primer","perturbed_primer",
        "primer_pos","species_matched","species_failed",
        "tm","hairpin_tm","dimer_tm","mfe","structure",
        "perturb_score","perturb_mean","perturb_min","perturb_frac"
    ]

    with open("plas_rtlamp_perturbation.csv","w",newline="") as f4:
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

    print("Pipeline complete.")
