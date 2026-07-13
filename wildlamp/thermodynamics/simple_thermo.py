import primer3.bindings as p3

def thermo(primer):
    hp = p3.calc_hairpin(primer).dg
    hd = p3.calc_homodimer(primer).dg
    return hp, hd
