"""
fluidtools.similitude — dimensional analysis (Buckingham Pi) and model-test
scaling.

Industry analog: wind-tunnel and tow-tank test planning. Unlike the other
tools, no commercial package does this *for* engineers — it's the analysis
they do on paper before spending tunnel time. This module automates the
linear algebra so you can focus on choosing variables, which is the actual
engineering.

The math: each variable's dimensions form a column of exponents over the
base dimensions (M, L, T, K). Dimensionless groups are vectors in the NULL
SPACE of that matrix. The number of groups = (variables) - (rank), which is
the Buckingham Pi theorem.
"""
from fractions import Fraction
import numpy as np

BASE_DIMS = ("M", "L", "T", "K")

# Dimension signatures of common physical quantities.
DIMS = {
    "length":          {"L": 1},
    "area":            {"L": 2},
    "volume":          {"L": 3},
    "time":            {"T": 1},
    "velocity":        {"L": 1, "T": -1},
    "acceleration":    {"L": 1, "T": -2},
    "gravity":         {"L": 1, "T": -2},
    "angular_velocity": {"T": -1},
    "frequency":       {"T": -1},
    "flow_rate":       {"L": 3, "T": -1},
    "mass":            {"M": 1},
    "mass_flow":       {"M": 1, "T": -1},
    "density":         {"M": 1, "L": -3},
    "force":           {"M": 1, "L": 1, "T": -2},
    "pressure":        {"M": 1, "L": -1, "T": -2},
    "stress":          {"M": 1, "L": -1, "T": -2},
    "energy":          {"M": 1, "L": 2, "T": -2},
    "power":           {"M": 1, "L": 2, "T": -3},
    "viscosity":       {"M": 1, "L": -1, "T": -1},   # dynamic, mu
    "kinematic_viscosity": {"L": 2, "T": -1},
    "surface_tension": {"M": 1, "T": -2},
    "temperature":     {"K": 1},
    "specific_heat":   {"L": 2, "T": -2, "K": -1},
    "conductivity":    {"M": 1, "L": 1, "T": -3, "K": -1},
    "dimensionless":   {},
}


def _dims_of(spec):
    """Accept a DIMS key string or an exponent dict."""
    if isinstance(spec, str):
        try:
            return DIMS[spec]
        except KeyError:
            raise KeyError(f"Unknown quantity '{spec}'. "
                           f"Known: {sorted(DIMS)} — or pass an exponent "
                           "dict like {'M':1,'L':-1,'T':-2}.")
    return dict(spec)


def _rref(M):
    """Reduced row echelon form over exact Fractions. Returns (R, pivots)."""
    M = [row[:] for row in M]
    rows, cols = len(M), len(M[0]) if M else 0
    pivots, r = [], 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if pivot is None:
            continue
        M[r], M[pivot] = M[pivot], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return M, pivots


def pi_groups(variables, repeating=None, verbose=True):
    """
    Find the dimensionless groups for a variable set.

    Parameters
    ----------
    variables : dict of {name: quantity}, where quantity is a DIMS key
                ('velocity', 'density', ...) or an exponent dict.
    repeating : optional list of variable names to prefer as the repeating
                set (each remaining variable then appears in exactly one
                group, textbook-style). They must be dimensionally
                independent or a ValueError is raised.

    Returns
    -------
    list of dicts {name: integer exponent}, one per group. Also prints the
    groups in readable form when verbose.

    Example
    -------
    >>> pi_groups({'F': 'force', 'V': 'velocity', 'D': 'length',
    ...            'rho': 'density', 'mu': 'viscosity'},
    ...           repeating=['rho', 'V', 'D'])
    Pi_1 = F / (rho V^2 D^2)     <- a drag coefficient
    Pi_2 = mu / (rho V D)        <- 1/Reynolds
    """
    names = list(variables)
    if repeating:
        missing = [r for r in repeating if r not in names]
        if missing:
            raise ValueError(f"repeating variables not in set: {missing}")
        names = list(repeating) + [n for n in names if n not in repeating]

    dim_maps = {n: _dims_of(variables[n]) for n in names}
    used_dims = [d for d in BASE_DIMS if any(m.get(d) for m in dim_maps.values())]

    # dimension matrix: rows = base dims, cols = variables
    A = [[Fraction(dim_maps[n].get(d, 0)) for n in names] for d in used_dims]
    R, pivots = _rref(A)
    rank = len(pivots)
    n_groups = len(names) - rank

    if repeating and any(p >= len(repeating) for p in pivots[:min(rank, len(repeating))]) \
            or (repeating and len(repeating) != rank):
        # repeating set must supply exactly the pivot columns
        if repeating and (len(repeating) != rank or set(pivots) != set(range(len(repeating)))):
            raise ValueError(
                f"The repeating set {repeating} is not a valid dimensionally "
                f"independent basis (rank = {rank}). Choose {rank} variables "
                "that between them contain all the dimensions present.")

    if n_groups == 0:
        if verbose:
            print("No dimensionless groups: the variables are dimensionally "
                  "independent (did you forget the dependent variable?).")
        return []

    free_cols = [c for c in range(len(names)) if c not in pivots]
    groups = []
    for fc in free_cols:
        # null-space basis vector: free var exponent = 1
        vec = [Fraction(0)] * len(names)
        vec[fc] = Fraction(1)
        for r_i, pc in enumerate(pivots):
            vec[pc] = -R[r_i][fc]
        # clear denominators -> integer exponents
        den = np.lcm.reduce([f.denominator for f in vec]) if vec else 1
        ints = [int(f * den) for f in vec]
        g = np.gcd.reduce([abs(i) for i in ints if i != 0])
        ints = [i // g for i in ints]
        groups.append({n: e for n, e in zip(names, ints) if e != 0})

    if verbose:
        for i, grp in enumerate(groups, 1):
            print(f"Pi_{i} = {format_group(grp)}")
    return groups


def format_group(group):
    """Render {name: exp} as a readable 'num / (den)' string."""
    num = [f"{n}" + (f"^{e}" if e != 1 else "") for n, e in group.items() if e > 0]
    den = [f"{n}" + (f"^{-e}" if e != -1 else "") for n, e in group.items() if e < 0]
    s = " ".join(num) if num else "1"
    if den:
        s += " / (" + " ".join(den) + ")"
    return s


def check_dimensionless(group, variables):
    """True if the exponent set is dimensionless — use it to verify by hand."""
    total = {}
    for n, e in group.items():
        for d, p in _dims_of(variables[n]).items():
            total[d] = total.get(d, 0) + p * e
    return all(v == 0 for v in total.values())


# ----------------------------------------------------------------------------
# Model-test scaling
# ----------------------------------------------------------------------------

def model_scale(law, Lr, rho_r=1.0, mu_r=1.0, g_r=1.0, sigma_r=1.0, a_r=1.0):
    """
    Scale ratios (model/prototype) for a similarity law.

    Parameters
    ----------
    law : 'reynolds' | 'froude' | 'mach' | 'weber'
    Lr  : length ratio, L_model / L_prototype (e.g. 1/20)
    rho_r, mu_r, g_r, sigma_r, a_r : fluid-property / gravity / sonic-speed
        ratios (model/prototype), default 1 (same fluid, same planet).

    Returns dict of ratios: V (velocity), t (time), f (frequency),
    Q (flow rate), F (force), P (power). Force uses dynamic-pressure scaling
    F_r = rho_r V_r^2 L_r^2, valid when the matched law governs.
    """
    law = law.lower()
    if law == "reynolds":
        Vr = (mu_r / rho_r) / Lr
    elif law == "froude":
        Vr = np.sqrt(g_r * Lr)
    elif law == "mach":
        Vr = a_r
    elif law == "weber":
        Vr = np.sqrt(sigma_r / (rho_r * Lr))
    else:
        raise ValueError("law must be reynolds, froude, mach, or weber")
    Fr = rho_r * Vr**2 * Lr**2
    return {"V": Vr, "t": Lr / Vr, "f": Vr / Lr, "Q": Vr * Lr**2,
            "F": Fr, "P": Fr * Vr}


def scale_report(law, Lr, **kw):
    """Human-readable version of model_scale for lab planning."""
    r = model_scale(law, Lr, **kw)
    print(f"{law.title()} similarity at length ratio Lr = {Lr:g}:")
    print(f"  velocity ratio  V_m/V_p = {r['V']:.4g}")
    print(f"  flow-rate ratio Q_m/Q_p = {r['Q']:.4g}")
    print(f"  force ratio     F_m/F_p = {r['F']:.4g}"
          f"   (prototype force = model force x {1/r['F']:.4g})")
    print(f"  power ratio     P_m/P_p = {r['P']:.4g}")
    print(f"  time ratio      t_m/t_p = {r['t']:.4g}")
    return r


# ----------------------------------------------------------------------------
# Named numbers & regime identification — "which physics is in charge?"
# ----------------------------------------------------------------------------
G = 9.81


def reynolds(rho, V, L, mu):
    """Re = rho V L / mu (inertia vs viscosity). L is the context length:
    pipe diameter, plate length, sphere diameter, ..."""
    return rho * V * L / mu


def froude(V, L):
    """Fr = V / sqrt(g L) (inertia vs gravity: channels, waves, ships)."""
    return V / np.sqrt(G * L)


def mach(V, c=343.0):
    """Ma = V / c (compressibility; c = 343 m/s for air at 20 C)."""
    return V / c


def weber(rho, V, L, sigma):
    """We = rho V^2 L / sigma (inertia vs surface tension: drops, sprays)."""
    return rho * V**2 * L / sigma


def regime(number, value, context="pipe"):
    """
    The physical regime a dimensionless number implies — i.e., which model
    you are allowed to use. Returns a short verdict string.

    number : 'reynolds' | 'froude' | 'mach' | 'weber'
    context (Reynolds only): 'pipe' | 'plate' | 'sphere'
    """
    n = number.lower()
    if n == "reynolds":
        if context == "pipe":
            if value < 2300:  return "laminar pipe flow (Poiseuille applies)"
            if value < 4000:  return "transitional (avoid designing here)"
            return "turbulent pipe flow (Colebrook/Moody territory)"
        if context == "plate":
            return ("laminar boundary layer over the whole plate"
                    if value < 5e5 else
                    "transition occurs on the plate: mixed laminar/turbulent")
        if context == "sphere":
            if value < 1:      return "creeping (Stokes) flow: Cd = 24/Re"
            if value < 1000:   return "intermediate: use a correlation"
            if value < 3.5e5:  return "Newton plateau: Cd ~ 0.44"
            return "post-drag-crisis: Cd ~ 0.1"
        raise ValueError("context must be pipe, plate, or sphere")
    if n == "froude":
        if value < 1:  return ("subcritical: gravity waves outrun the flow, "
                               "so downstream conditions reach upstream")
        if value > 1:  return ("supercritical: the flow outruns its own "
                               "waves; upstream cannot hear downstream")
        return "critical (Fr = 1): the choking condition"
    if n == "mach":
        if value < 0.3: return ("incompressible treatment is fine "
                                "(density error under ~5%)")
        if value < 0.8: return "subsonic but compressible: density matters"
        if value < 1.2: return "transonic: expect shocks, expect trouble"
        return "supersonic: shock waves are part of the flow"
    if n == "weber":
        return ("surface tension dominates: drops hold together"
                if value < 10 else
                "inertia dominates: drops deform and break up")
    raise ValueError("number must be reynolds, froude, mach, or weber")
