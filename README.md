# wildlamp: a tool for identification of primer binding sites among conserved gene families of primers. 

wildlamp/
|
|-- alignment/
|    |-- conservation.py # compute_conservation
|    |-- load_seq.py # load_alignment and clean_alignment
|
|-- data/
|    |--FILL IN LATER
|
|-- failures/
|    |-- debug.py # debug_no_matches
|    |-- failure_mode.py # analyze_failure_modes
|    |-- suggest.py # suggest_degenerate
|
|-- perturbation/
|    |-- analysis.py # perturbation_analysis
|    |-- mutate.py  # mutate_base and perturb_primer
|    |-- robustness.py  # built_robustness_matrix, compute_perturbation_robustness_scores
|
|-- pipeline/
|    |-- init_primers.py    # init_primers
|
|-- primers/
|    |-- amplicons.py   # amplicon_heuristic and amplicon_score
|    |-- clash.py   # primer_clash, build_clash_map
|    |-- detect_lamp.py # detect_lamp_sets
|    |-- lamp_score.py  # primer_basic_score, primer_thermo_score, primer_struct_score
|    |                  # primer_conservation_score, lamp_score, rtlamp_score
|    |-- mismatch.py    # mismatch_cost
|    |-- scan.py    # primer_ok, scan_primers
|    |-- species.py # primer_species_coverage, filter_primers_by_species
|
|-- scripts/ # contains the scripts used in the study
|    |-- acc2spp.py # extract_accession, accession_to_gene_id, get_species_from_accession, 
|    |              # extract_accessions_from_csv, write_species_lookup
|    |-- fasta-stats.py # no commands - just basic code
|    |-- final_dilution-LoD_3.py    # rc, build_FIP, build_BIP, load_gene, best_match_mismatches, 
|    |                              # compute_hairpin, compute_homodimer, compute_heterodimer,
|    |                              # has_3prime_dimer, compute_primer_features, 
|    |                              # fraction_sequestered, fraction_free, binding_efficiency,
|    |                              # dimer_sink_factor, lamp_tm_score, 
|    |                              # lamp_primer_biophysical_score, rt_lamp_rate, simulate_rt_qlamp
|    |                              # color_hnb_from_copies, color_pr_from_copies, delta_rgb, 
|    |                              # simulate_quant_curve, fit_standard_curve, 
|    |                              # copies_from_Tt_colour, plot_multi_gene_curve, plot_delta_rgb,
|    |                              # generate_all_thermo_tables, write_thermo_csv, 
|    |                              # generate_standard_curve_summary, main
|    |-- final-validation.py    # load_gene, gc_content, fold_rna, best_match, binding_accessibility
|    |                          # thermo, analyse_gene, write_gene_tables, write_gc_table, 
|    |                          # write_long_all, main
|    |-- LAMP_Assay.py  # load_alignment, clean_alignment, compute_conservation, gc, bad_run, 
|    |                  # mismatch_cost, primers_ok, scan_primers, primer_species_coverage, 
|    |                  # filter_primers_by_species, species_coverage, primers_clash, 
|    |                  # build_clash_map, amplicon_heuristic, detect_lamp_sets, evaluate_primer3,
|    |                  # vienna_fold_dna, validate_primer_vienna, validate_lamp_set_vienna, 
|    |                  # primer_basic_score, primer_thermo_score, primer_struct_score, 
|    |                  # primer_conservation_score, amplicon_score, lamp_score, rtlamp_score,
|    |                  # suggest_degenerate, analyze_failure_modes, mutate_base, perturb_primer,
|    |                  # perturbation_analysis, build_robustness_matrix, 
|    |                  # compute_perturbation_robustness_scores, debug_no_matches, main
|    |-- LAMP-hits.py   # blast_sequences, biologically_meaningful
|    |-- pertubation_final.py   # best_mismatch_score, 
|    |-- primer_optimise.py # load_alignment, clean_alignment, gc, bad_run, rc, mismatch_cost, 
|                           # primer_species_coverage, thermo_score, primer_ok, optimise_primer, 
|                           # vienna_fold_dna, vienna_validate_set, mutate_base, perturb_primer,
|                           # perturbation_analysis, compute_robustness_score, build_full_lset, main
|
|-- thermodynamics/
|    |-- primer3_eval.py    # evaluate_primer3
|    |-- vienna_eval.py     # vienna_fold_dna, validate_primer_vienna, validate_lamp_set_vienna
|
|-- utils/
|    |-- seq.py # gc, bad_run, rc
|
|-- __init__.py 
|-- LICENSE # open access code - must cite
|-- README.md