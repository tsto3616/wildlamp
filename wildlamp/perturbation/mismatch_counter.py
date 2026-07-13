def best_match(primer, template):
    Lp = len(primer)
    Lt = len(template)
    best_mm = Lp
    best_pos = -1

    for i in range(Lt - Lp + 1):
        window = template[i:i+Lp]
        mm = sum(1 for a, b in zip(primer, window) if a != b)
        if mm < best_mm:
            best_mm = mm
            best_pos = i
            if mm == 0:
                break

    return best_mm, best_pos
