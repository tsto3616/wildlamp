from ..utils.seq import rc

def build_FIP(F1: str, F2: str) -> str:
    return rc(F1) + F2

def build_BIP(B1: str, B2: str) -> str:
    return rc(B1) + B2