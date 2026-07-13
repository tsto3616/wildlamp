def binding_accessibility(structure, start, end):
    region = structure[start:end]
    unpaired = region.count('.')
    return (unpaired / len(region)) * 100
