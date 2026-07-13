import csv

def write_outputs(rows, filename="final_valid_output.csv"):
    fieldnames = [
        "Primer","Sequence","GC%","Mismatches",
        "Binding_position","Accessibility_%",
        "Hairpin_dG","Homodimer_dG","Transcript_MFE"
    ]

    with open(filename, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"Output written to {filename}")
