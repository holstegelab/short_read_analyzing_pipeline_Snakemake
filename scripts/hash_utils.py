FNV_OFFSET = 0xcbf29ce484222325
FNV_PRIME = 0x100000001b3

def fnv1a64(data: str) -> int:
    h = FNV_OFFSET
    for ch in data:
        h ^= ord(ch)
        h = (h * FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
    return h

def normalize_qual(s: str) -> str:
    # dragmap sometimes replaces '#' with '!' in quality
    return s.replace('#', '!')
