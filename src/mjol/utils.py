def load_attributes(s: str, kv_sep: str = '=') -> dict:
    return {
        k.strip(): v.strip()
        for x in s.strip().split(';') if x
        for k, v in [x.strip().split(kv_sep, 1)]
    }