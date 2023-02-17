

def SI_to_AU(r=0, v=0):
    if r is not None:
        r = r/149597871
    if v is not None:
        v = v/149597871 * 3600 * 24

    return r, v


def AU_to_SI(r=None, v=None):
    if r is not None:
        r = r*149597871

    if v is not None:
        v = v * 149597871 / (3600 * 24)

    return r, v