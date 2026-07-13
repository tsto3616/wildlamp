def gc_content(seq):
    gc = seq.count("G") + seq.count("C")
    return (gc / len(seq)) * 100
