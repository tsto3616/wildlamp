# wildlamp: a tool for identification of primer binding sites among conserved gene families of primers. 

This package enables the construction of the LAMP primers, their optimisation and checking for quality/ cross-reactivity. The package can be downloaded as: 

```
pip install wildlamp
```

The main command for the construction of LAMP primers is the command:

```
wl.init_primers(alignment_path, output="fill in for your specified output.csv")
```

Where the user specifies the file path for the aligned set of sequences and the rest occurs by computer. 

Following on from that the user can retrieve the species from the aligned sequence accessions provided they are from NCBI, using the command: 

```
wl.species_fetch(csv_file)
```

The primers can then be checked for hits against the genome of a species of their choice using either the csv2hits or the primer2hits command. 

```
wl.csv2hits(csv_file, taxid, output="fill in your preferred output prefix.csv")

# or use this command if you want to specify your own primers independent of the csv
panel: [
    {"FIP": "ATCG...", "BIP": "GGTA..."},
    {"FIP": "...", "BIP": "..."},
    ...
]

wl.primers2hits(panel, taxid, output="preferred specified csv output")
```

The final primers can be generated through the command: 

```
wl.final_primers(alignment_path, primer_csv_path, output="preferred csv output")
```

The final perturbation checks can be performed with the command below - with both lightweight and heavy perturbation analysis. Heavyweight analysis is primer-thermodynamic stress test and the lightweight is a biological/ transcript specificity/ sensitivity analysis... or more simply a multitransript mismatch scan. The first transcript in the fasta file is assumed to be the reference for the heavyweight analysis. Note only use unaligned fasta files for this section. The command can be called by: 

```
wl.final_perturbation(unaligned_fasta_path, 
    primer_csv_path_from_earlier, 
    output="perturbation_output.csv path")
```

Final validation check for the primer sets:

```
wl.final_valid(primer_csv_path, unaligned_fasta_path, output="output csv file")
```

Once the primers are settled, run the csv file through this code to generate the standard curve. But the csv must be edited first to have a column "Gene" and "fasta":
``` 
wl.stand_curves(primer_gene_csv_path, output="standard_curve.csv, but fill in yourself")
```

If one wants to run the standard curve quantification in the absence of experimental data the following code will do that, but first caveates, must fill in the user R, G, B values, the gene name in the csv, the curve_csv from the step above, the output and the dilution factor. The dilution factor is 1:X where X is the user specified value: 
```
wl.estimate_copies(user_R, user_G, user_B,
                    gene_name, curve_csv,
                    output="final_csv.csv",
                    dilution_factor=1)
```

File directory is below:

```
wildlamp/
|
|-- alignment/
|    |-- conservation.py # compute_conservation
|    |-- load_seq.py # load_alignment and clean_alignment
|
|-- data/
|    |-- csv_results/
|    |   |-- FILL IN LATER
|    |
|    |-- study_fasta/
|    |   |-- DIO1_N_fur.fasta # DIO1 reference for followup
|    |   |-- DIO1_pinniped_aligned.fasta    # DIO1 aligned
|    |   |-- DIO1_pinniped.fasta    # unaligned
|    |   |-- DIO2_N_fur.fasta   # DIO2 reference
|    |   |-- DIO2_pinniped_aligned.fasta    # DIO2 aligned
|    |   |-- DIO2_pinniped.fasta    # DIO2 unaligned
|    |   |-- DIO3_N_fur.fasta   # DIO3 reference
|    |   |-- DIO3_pinnipedia.fas    # DIO3 ALIGNED
|    |   |-- DIO3_pinnipeds.fasta   # DIO3 UNALIGNED 
|    |   |-- HSP70_all_aligned.fasta    # HSP70 aligned
|    |   |-- HSP70_all.fasta   # HSP70 unaligned
|    |   |-- HSP70_N_fur.fasta  # HSP70 reference
|    |   |-- plasminogen_receptor_aligned.fasta # plasminogen receptor aligned
|    |   |-- plasminogen_receptor.fasta # plasminogen receptor unaligned (note no ref - not studied)
|    |   |-- slc25a43_aligned.fasta # slc25a43 aligned
|    |   |-- slc25a43_N_fur.fasta   # slc25a43 reference
|    |   |-- slc25a43.fasta # unaligned
|
|-- failures/
|    |-- debug.py # debug_no_matches
|    |-- failure_mode.py # analyze_failure_modes
|    |-- suggest.py # suggest_degenerate
|
|-- perturbation/
|    |-- analysis.py # perturbation_analysis
|    |-- best_mismatch_score.py # best_mismatch_score, find_best_matches
|    |-- mismatch_counter.py    # best_match
|    |-- mutate.py  # mutate_base and perturb_primer
|    |-- robustness.py  # built_robustness_matrix, compute_perturbation_robustness_scores
|
|-- pipeline/
|    |-- conc.py    # estimate_copies
|    |-- csv2hits.py    # csv2hits
|    |-- final_perturb.py   # fianl_perturbation
|    |-- final_primers.py   # final_primers
|    |-- init_primers.py    # init_primers
|    |-- primer2hits.py     # primers2hits
|    |-- species_fetch.py   # species_fetch
|    |-- stand_curve.py # stand_curve
|    |-- validate.py    # final_validate
|
|-- primer_hits/
|    |-- crossreactive.py   # blast_sequence, biologically_meaningful
|
|-- primers/
|    |-- amplicons.py   # amplicon_heuristic and amplicon_score
|    |-- analyse_primer_rows.py # analyse_row
|    |-- build_full_lset.py # build_full_lset
|    |-- clash.py   # primer_clash, build_clash_map
|    |-- detect_lamp.py # detect_lamp_sets
|    |-- lamp_score.py  # primer_basic_score, primer_thermo_score, primer_struct_score
|    |                  # primer_conservation_score, lamp_score, rtlamp_score
|    |-- mismatch.py    # mismatch_cost
|    |-- optimise_primer.py # optimise_primer
|    |-- primer_from_csv.py # load_primers_from_csv
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
|-- simulations/
|    |-- curve_fit.py   # fit_standard_curve, copies_from_Tt
|    |-- dye.py # color_pr_from_copies, delta_rgb, parameters: PR_BASELINE
|    |-- model_rtlamp.py    # best_match_mismatches, compute_hairpin, compute_homodimer, 
|    |                      # compute_heterodimer, has_3prime_dimer, compute_primer_features, 
|    |                      # fraction_sequestered, fraction_free, binding_efficiency, 
|    |                      # dimer_sink_factor, lamp_tm_score, lamp_primer_biophysical_score,
|    |                      # rt_lamp_rate, simulate_rt_qlamp
|    |-- sim_curve.py   # simulate_quant_curve_pr, parameters: QUANT_CONC
|
|-- species/
|    |-- accession.py   # extract_accession, extract_accessions_from_csv
|    |-- ncbi_lookup.py # accession_to_gene_id, get_species_from_accession
|    |-- table.py   # write_species_lookup
|
|-- thermodynamics/
|    |-- fold_rna.py    # fold_rna
|    |-- primer_access.py   # binding_accessibility
|    |-- primer3_eval.py    # evaluate_primer3
|    |-- scoring.py # thermo_score
|    |-- simple_thermo.py   # thermo
|    |-- vienna_eval.py     # vienna_fold_dna, validate_primer_vienna, validate_lamp_set_vienna
|
|-- utils/
|    |-- build_primers.py   # build_BIP, build_FIP
|    |-- gc_content.py # gc_content
|    |-- gene_primer.py # load_gene_primer_csv
|    |-- load_seq.py    # load_primers_csv, load_fasta
|    |-- seq.py # gc, bad_run, rc
|    |-- valid_csv_write.py # write_outputs
|
|-- __init__.py 
|-- LICENSE # open access code - must cite
|-- README.md
```