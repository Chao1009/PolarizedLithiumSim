"""Backend selector: TOY parameterizations vs LHAPDF grids (via parton)."""


def get_backends(pdf="toy", nuclear=None, r_func=None):
    """Return dict with 'base' (F2 source), 'g1' (g1 model), 'tag'.

    `nuclear` is an `beams.Ion`; give one and the dict gains a 'nuclear'
    key holding the whole-nucleus F2A object that goes with `pdf` --
    `structure.NuclearF2` over the free-nucleon backend on 'toy',
    `structure.NuclearF2FromGrid` (the EPPS21nlo_CT18Anlo_Li6 nuclear set)
    on 'grid'.  `r_func(x, q2) -> R = sigma_L/sigma_T` is threaded to both
    R consumers this function builds, the g1 model and that nuclear F2.

    Both default to None, and at None the returned dict is exactly the
    three keys and the three objects it has always had -- so every one of
    the nine existing call sites is bit-for-bit unchanged.
    """
    if pdf == "grid":
        from .polarized import PartonG1
        from .structure import PartonF2
        base = PartonF2()  # CT18NLO
        out = {"base": base, "g1": PartonG1(base=base, r_func=r_func),
               "tag": "grid"}
        if nuclear is not None:
            from .structure import NuclearF2FromGrid
            out["nuclear"] = NuclearF2FromGrid(nuclear, r_func=r_func)
        return out
    from .polarized import ToyG1
    from .structure import ToyF2
    base = ToyF2()
    out = {"base": base, "g1": ToyG1(base=base, r_func=r_func), "tag": "toy"}
    if nuclear is not None:
        from .structure import NuclearF2
        out["nuclear"] = NuclearF2(nuclear, base=base, r_func=r_func)
    return out
