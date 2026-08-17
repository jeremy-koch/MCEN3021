"""
pipeflow.py — A teaching toolbox for steady, incompressible pipe flow analysis.

Designed for undergraduate fluid mechanics. Students assemble piping systems
from components (pipes, fittings, pumps, turbines); the toolbox handles the
energy equation, friction factors, and numerical solving.

Conventions
-----------
* SI units throughout: meters, seconds, kilograms, Pascals, Watts.
* The energy equation between inlet (1) and outlet (2) of a system:

    p1/(rho*g) + a1 V1^2/(2g) + z1 + h_pump = p2/(rho*g) + a2 V2^2/(2g) + z2 + h_turbine + h_L

* h_L = sum of major losses (Darcy-Weisbach) and minor losses (K*V^2/2g).
* Friction factor: laminar f = 64/Re; turbulent f from the Colebrook equation
  (solved numerically). Swamee-Jain is available for comparison.

Author: course toolbox, v1.0
"""

import numpy as np
from dataclasses import dataclass, field
from scipy.optimize import brentq

G = 9.81  # m/s^2


# ----------------------------------------------------------------------------
# Fluids
# ----------------------------------------------------------------------------

@dataclass
class Fluid:
    """
    A fluid with density rho [kg/m^3], dynamic viscosity mu [Pa*s], and
    (optionally) vapor pressure pv [Pa, absolute] for cavitation analysis.
    """
    rho: float
    mu: float
    name: str = "fluid"
    pv: float = None   # vapor pressure [Pa abs]; None = unknown

    @property
    def nu(self):
        """Kinematic viscosity [m^2/s]."""
        return self.mu / self.rho

    def __repr__(self):
        pv = f", pv={self.pv:.0f} Pa" if self.pv is not None else ""
        return f"Fluid({self.name}: rho={self.rho} kg/m3, mu={self.mu:.3e} Pa*s{pv})"


def water(T_celsius=20.0):
    """
    Liquid water properties at 1 atm, valid ~0-100 C (correlation fits).
    Includes vapor pressure (Antoine equation) so cavitation checks work
    automatically -- and so that TEMPERATURE becomes a design variable.
    """
    T = T_celsius
    if not (0 <= T <= 100):
        raise ValueError("water() correlation valid for 0-100 C")
    rho = 999.84 + 0.0673 * T - 0.00894 * T**2 + 8.78e-5 * T**3 - 6.62e-7 * T**4
    # Vogel-type viscosity fit
    mu = 2.414e-5 * 10 ** (247.8 / (T + 133.15))
    # Antoine equation (1-100 C), pressure in mmHg -> Pa
    pv = 133.322 * 10 ** (8.07131 - 1730.63 / (233.426 + T))
    return Fluid(rho=rho, mu=mu, name=f"water @ {T:g} C", pv=pv)


# Common fluids at ~20 C for convenience
AIR_20C = Fluid(rho=1.204, mu=1.825e-5, name="air @ 20 C")
GASOLINE = Fluid(rho=737.0, mu=2.92e-4, name="gasoline @ 20 C")
SAE30_OIL = Fluid(rho=891.0, mu=0.29, name="SAE 30 oil @ 20 C")


# ----------------------------------------------------------------------------
# Friction factors
# ----------------------------------------------------------------------------

def friction_factor(Re, rel_roughness=0.0, method="colebrook"):
    """
    Darcy friction factor.

    Parameters
    ----------
    Re : Reynolds number (V*D/nu), must be > 0.
    rel_roughness : epsilon/D (dimensionless).
    method : 'colebrook' (implicit, solved numerically) or 'swamee-jain'
             (explicit approximation). Laminar flow (Re < 2300) always uses
             f = 64/Re regardless of method.

    Notes
    -----
    In the transition region 2300 < Re < 4000 the flow is unpredictable;
    we apply the turbulent correlation there (a common, conservative choice)
    but students should treat such results with suspicion.
    """
    if Re <= 0:
        raise ValueError("Reynolds number must be positive.")
    if Re < 2300.0:
        return 64.0 / Re
    if method == "swamee-jain":
        return 0.25 / (np.log10(rel_roughness / 3.7 + 5.74 / Re**0.9)) ** 2

    # Colebrook: 1/sqrt(f) = -2 log10( eps/D/3.7 + 2.51/(Re sqrt(f)) )
    def residual(f):
        return 1.0 / np.sqrt(f) + 2.0 * np.log10(
            rel_roughness / 3.7 + 2.51 / (Re * np.sqrt(f))
        )

    return brentq(residual, 1e-4, 0.2, xtol=1e-12)


# ----------------------------------------------------------------------------
# Roughness and minor-loss databases
# ----------------------------------------------------------------------------

# Absolute roughness epsilon [m] for common pipe materials (typical new values)
ROUGHNESS = {
    "drawn_tubing":     1.5e-6,
    "commercial_steel": 4.5e-5,
    "galvanized_iron":  1.5e-4,
    "cast_iron":        2.6e-4,
    "concrete":         1.0e-3,   # midpoint of 0.3-3 mm range
    "riveted_steel":    3.0e-3,
    "pvc":              1.5e-6,   # smooth plastic
    "copper":           1.5e-6,
    "smooth":           0.0,
}

# Loss coefficients K for fittings and components (typical textbook values).
# h_minor = K * V^2 / (2g), with V the velocity in the attached pipe.
FITTINGS = {
    "elbow_90_standard":     0.9,
    "elbow_90_long_radius":  0.6,
    "elbow_45":              0.4,
    "tee_line_flow":         0.2,   # flow straight through the run
    "tee_branch_flow":       1.0,   # flow out the branch
    "return_bend_180":       1.5,
    "gate_valve_open":       0.15,
    "gate_valve_half":       2.1,
    "globe_valve_open":      10.0,
    "angle_valve_open":      2.0,
    "ball_valve_open":       0.05,
    "check_valve_swing":     2.0,
    "entrance_sharp":        0.5,
    "entrance_reentrant":    0.8,
    "entrance_rounded":      0.03,
    "exit":                  1.0,   # into a reservoir: all KE lost
    "strainer":              2.0,
    "flow_meter_orifice":    2.5,
}


# ----------------------------------------------------------------------------
# Components
# ----------------------------------------------------------------------------

class Pipe:
    """
    A straight pipe segment, optionally carrying minor-loss fittings.

    Parameters
    ----------
    L : length [m]
    D : inside diameter [m]
    roughness : absolute roughness epsilon [m], OR a material-name string
                from the ROUGHNESS database (e.g. 'commercial_steel').
    name : optional label used in reports.

    Examples
    --------
    >>> p = Pipe(L=30, D=0.05, roughness='commercial_steel', name='suction line')
    >>> p.add_fitting('elbow_90_standard', n=2)
    >>> p.add_fitting(K=0.3)   # custom coefficient
    """

    def __init__(self, L, D, roughness=0.0, name="pipe"):
        if L <= 0 or D <= 0:
            raise ValueError("Pipe length and diameter must be positive.")
        self.L = float(L)
        self.D = float(D)
        if isinstance(roughness, str):
            try:
                roughness = ROUGHNESS[roughness]
            except KeyError:
                raise KeyError(
                    f"Unknown material '{roughness}'. "
                    f"Options: {sorted(ROUGHNESS)}"
                )
        self.eps = float(roughness)
        self.name = name
        self.fittings = []          # list of (label, K, count)
        self.f_method = "colebrook"

    # -- geometry ------------------------------------------------------------
    @property
    def area(self):
        return np.pi * self.D**2 / 4.0

    @property
    def rel_roughness(self):
        return self.eps / self.D

    @property
    def K_total(self):
        return sum(K * n for (_, K, n) in self.fittings)

    # -- fittings ------------------------------------------------------------
    def add_fitting(self, kind=None, n=1, K=None):
        """Attach n fittings by database name, or a custom K value."""
        if K is not None:
            self.fittings.append((kind or "custom", float(K), n))
        elif kind is not None:
            try:
                self.fittings.append((kind, FITTINGS[kind], n))
            except KeyError:
                raise KeyError(
                    f"Unknown fitting '{kind}'. Options: {sorted(FITTINGS)}"
                )
        else:
            raise ValueError("Provide a fitting name or a K value.")
        return self  # allow chaining

    # -- hydraulics ----------------------------------------------------------
    def velocity(self, Q):
        return Q / self.area

    def reynolds(self, Q, fluid):
        return fluid.rho * self.velocity(Q) * self.D / fluid.mu

    def friction_factor(self, Q, fluid):
        return friction_factor(self.reynolds(Q, fluid),
                               self.rel_roughness, self.f_method)

    def head_loss(self, Q, fluid):
        """Total head loss [m] = major (Darcy-Weisbach) + minor (sum K)."""
        if Q == 0:
            return 0.0
        V = self.velocity(abs(Q))
        f = self.friction_factor(abs(Q), fluid)
        h_major = f * (self.L / self.D) * V**2 / (2 * G)
        h_minor = self.K_total * V**2 / (2 * G)
        return h_major + h_minor

    def report(self, Q, fluid):
        """Print a breakdown of losses at flow rate Q."""
        V = self.velocity(Q)
        Re = self.reynolds(Q, fluid)
        f = self.friction_factor(Q, fluid)
        hv = V**2 / (2 * G)
        h_major = f * (self.L / self.D) * hv
        h_minor = self.K_total * hv
        regime = "laminar" if Re < 2300 else (
            "transitional (!)" if Re < 4000 else "turbulent")
        print(f"--- {self.name} ---")
        print(f"  V = {V:8.3f} m/s   Re = {Re:10.3e}  ({regime})")
        print(f"  f = {f:8.4f}      eps/D = {self.rel_roughness:.2e}")
        print(f"  major loss: {h_major:8.3f} m  (f L/D = {f*self.L/self.D:.2f})")
        print(f"  minor loss: {h_minor:8.3f} m  (sum K  = {self.K_total:.2f})")
        for label, K, n in self.fittings:
            print(f"      {n} x {label} (K={K})")
        print(f"  TOTAL:      {h_major + h_minor:8.3f} m")
        return h_major + h_minor

    def __repr__(self):
        return (f"Pipe({self.name}: L={self.L} m, D={self.D} m, "
                f"eps={self.eps:.2e} m, sumK={self.K_total:.2f})")


class Parallel:
    """
    Two or more pipe branches in parallel.

    All branches share the same head loss; the total flow divides among them.
    Given total Q, the toolbox finds the head loss h such that the branch
    flows sum to Q.

    Example
    -------
    >>> branch_a = Pipe(L=100, D=0.05, roughness='commercial_steel')
    >>> branch_b = Pipe(L=60,  D=0.04, roughness='commercial_steel')
    >>> par = Parallel([branch_a, branch_b], name='split section')
    """

    def __init__(self, branches, name="parallel"):
        if len(branches) < 2:
            raise ValueError("Parallel needs at least two branches.")
        self.branches = list(branches)
        self.name = name

    def _branch_flow(self, h, branch, fluid):
        """Flow through one branch given head loss h across it."""
        if h <= 0:
            return 0.0
        f_res = lambda Q: branch.head_loss(Q, fluid) - h
        # bracket: zero flow (zero loss) up to a flow whose loss exceeds h
        Q_hi = 1.0
        while branch.head_loss(Q_hi, fluid) < h:
            Q_hi *= 2.0
            if Q_hi > 1e6:
                raise RuntimeError("Could not bracket branch flow.")
        return brentq(f_res, 0.0, Q_hi, xtol=1e-12)

    def head_loss(self, Q, fluid):
        """Head loss across the parallel section for total flow Q."""
        if Q == 0:
            return 0.0
        Q = abs(Q)

        def residual(h):
            return sum(self._branch_flow(h, b, fluid) for b in self.branches) - Q

        # bracket h between ~0 and the loss if ALL flow went down each branch
        h_hi = max(b.head_loss(Q, fluid) for b in self.branches)
        return brentq(residual, 1e-12, h_hi * 1.0000001, xtol=1e-10)

    def flow_split(self, Q, fluid):
        """Return list of branch flow rates for total flow Q."""
        h = self.head_loss(Q, fluid)
        return [self._branch_flow(h, b, fluid) for b in self.branches]

    def report(self, Q, fluid):
        h = self.head_loss(Q, fluid)
        print(f"--- {self.name} (parallel, {len(self.branches)} branches) ---")
        print(f"  head loss across section: {h:.3f} m")
        for b, Qb in zip(self.branches, self.flow_split(Q, fluid)):
            print(f"  {b.name}: Q = {Qb*1000:8.3f} L/s "
                  f"({100*Qb/Q:5.1f} % of total)")
        return h

    def __repr__(self):
        return f"Parallel({self.name}: {len(self.branches)} branches)"


class Pump:
    """
    A pump adds head h_p(Q) to the flow.

    Three ways to define one:
      Pump(head=25)                          constant-head idealization
      Pump.from_coeffs(h0=40, a=1200)        h = h0 - a*Q^2 (Q in m^3/s)
      Pump.from_data(Q_pts, h_pts)           quadratic fit through test data

    efficiency : pump efficiency (0-1], used for power calculations.
    npsh_r : NPSH-required [m] -- a scalar, or a function of Q. Comes from
             the manufacturer datasheet. Needed for cavitation checks.
    z : pump centerline elevation [m]. Defaults to the system inlet
        elevation z1 if not given. Critical for suction analysis!
    """

    def __init__(self, head=None, curve=None, efficiency=0.75, name="pump",
                 npsh_r=None, z=None):
        self._head = head
        self._curve = curve
        self.eta = efficiency
        self.name = name
        self._npshr = npsh_r
        self.z = z

    @classmethod
    def from_coeffs(cls, h0, a, efficiency=0.75, name="pump",
                    npsh_r=None, z=None):
        return cls(curve=lambda Q: h0 - a * Q**2, efficiency=efficiency,
                   name=name, npsh_r=npsh_r, z=z)

    @classmethod
    def from_data(cls, Q_pts, h_pts, efficiency=0.75, name="pump",
                  npsh_r=None, z=None):
        """Fit h = c0 + c1*Q + c2*Q^2 through manufacturer test points."""
        c = np.polyfit(np.asarray(Q_pts, float), np.asarray(h_pts, float), 2)
        return cls(curve=lambda Q: np.polyval(c, Q), efficiency=efficiency,
                   name=name, npsh_r=npsh_r, z=z)

    def head(self, Q):
        if self._curve is not None:
            return float(max(self._curve(Q), 0.0))
        return float(self._head)

    def npsh_required(self, Q):
        """NPSH-required [m] at flow Q, from the manufacturer spec."""
        if self._npshr is None:
            return None
        if callable(self._npshr):
            return float(self._npshr(Q))
        return float(self._npshr)

    def power(self, Q, fluid):
        """Shaft power required [W] = rho g Q h / eta."""
        return fluid.rho * G * Q * self.head(Q) / self.eta

    # -- affinity laws ---------------------------------------------------
    def at_speed(self, ratio, name=None):
        """
        The same pump run at (ratio x) its rated speed, via the affinity
        laws: Q ~ N, h ~ N^2 (so P ~ N^3). NPSH-required scales ~ N^2.

        Example: pump.at_speed(0.8) is the pump on a VFD at 80% speed.
        """
        base_h = self._curve if self._curve is not None \
            else (lambda Q, h=self._head: h)
        new_curve = lambda Q: ratio**2 * base_h(Q / ratio)
        if self._npshr is None:
            new_npshr = None
        elif callable(self._npshr):
            new_npshr = lambda Q, f=self._npshr: ratio**2 * f(Q / ratio)
        else:
            new_npshr = ratio**2 * self._npshr
        return Pump(curve=new_curve, efficiency=self.eta,
                    name=name or f"{self.name} @ {100*ratio:.0f}% speed",
                    npsh_r=new_npshr, z=self.z)

    def trimmed(self, d_ratio, name=None):
        """
        The same pump with its impeller trimmed to (d_ratio x) the original
        diameter. To first order the affinity scaling is the same form as
        speed change: Q ~ d, h ~ d^2. (Real trims deviate a few percent --
        a good discussion point.) Trimming is permanent; a VFD is not.
        """
        return self.at_speed(d_ratio,
                             name=name or f"{self.name} trimmed to "
                                          f"{100*d_ratio:.0f}% impeller")

    def __repr__(self):
        kind = "curve" if self._curve is not None else f"h={self._head} m"
        return f"Pump({self.name}: {kind}, eta={self.eta})"


class PumpsInSeries:
    """
    Pumps in series: same Q through each, heads ADD. Used for high lifts.
    NPSH is governed by the FIRST pump (it sees the suction line).
    """

    def __init__(self, pumps, name="pumps in series"):
        self.pumps = list(pumps)
        self.name = name
        self.z = pumps[0].z
        self.eta = None  # heterogeneous; power computed per pump

    def head(self, Q):
        return sum(p.head(Q) for p in self.pumps)

    def npsh_required(self, Q):
        return self.pumps[0].npsh_required(Q)

    def power(self, Q, fluid):
        return sum(p.power(Q, fluid) for p in self.pumps)

    def __repr__(self):
        return f"PumpsInSeries({[p.name for p in self.pumps]})"


class PumpsInParallel:
    """
    Pumps in parallel: same head across each, flows ADD. Each pump is
    assumed to have a check valve (a pump asked for more head than its
    shutoff simply contributes zero flow rather than spinning backward).

    Key teaching point: two identical pumps in parallel do NOT double the
    system flow, because the operating point slides up the system curve.
    """

    def __init__(self, pumps, name="pumps in parallel"):
        self.pumps = list(pumps)
        self.name = name
        self.z = pumps[0].z
        self.eta = None

    def _pump_flow(self, h, pump):
        """Flow one pump delivers against head h (0 if h >= shutoff)."""
        if h >= pump.head(0.0):
            return 0.0
        Q_hi = 1e-3
        while pump.head(Q_hi) > h:
            Q_hi *= 2.0
            if Q_hi > 1e4:
                raise RuntimeError("Could not bracket pump runout flow.")
        return brentq(lambda Q: pump.head(Q) - h, 0.0, Q_hi, xtol=1e-12)

    def head(self, Q):
        """Head at which the pumps jointly deliver total flow Q."""
        if Q <= 0:
            return max(p.head(0.0) for p in self.pumps)
        h_shut = max(p.head(0.0) for p in self.pumps)

        def residual(h):
            return sum(self._pump_flow(h, p) for p in self.pumps) - Q

        if residual(0.0) < 0:   # beyond combined runout
            return 0.0
        return brentq(residual, 0.0, h_shut, xtol=1e-10)

    def flow_split(self, Q):
        h = self.head(Q)
        return [self._pump_flow(h, p) for p in self.pumps]

    def npsh_required(self, Q):
        """Worst-case NPSHr among pumps, each at its own branch flow."""
        vals = [p.npsh_required(Qi)
                for p, Qi in zip(self.pumps, self.flow_split(Q))]
        vals = [v for v in vals if v is not None]
        return max(vals) if vals else None

    def power(self, Q, fluid):
        return sum(p.power(Qi, fluid)
                   for p, Qi in zip(self.pumps, self.flow_split(Q)))

    def __repr__(self):
        return f"PumpsInParallel({[p.name for p in self.pumps]})"


PUMP_TYPES = (Pump, PumpsInSeries, PumpsInParallel)


class Turbine:
    """A turbine extracts head h_t from the flow. efficiency in (0,1]."""

    def __init__(self, head, efficiency=0.85, name="turbine"):
        self._head = head
        self.eta = efficiency
        self.name = name

    def head(self, Q):
        return float(self._head)

    def power(self, Q, fluid):
        """Shaft power delivered [W] = eta * rho g Q h."""
        return self.eta * fluid.rho * G * Q * self.head(Q)

    def __repr__(self):
        return f"Turbine({self.name}: h={self._head} m, eta={self.eta})"


# ----------------------------------------------------------------------------
# The system
# ----------------------------------------------------------------------------

class PipeSystem:
    """
    A flow path from point 1 to point 2, assembled from components in series.

    Endpoint conditions (defaults describe two large open reservoirs):
      z1, z2 : elevations [m]
      p1, p2 : gauge pressures [Pa]
      V1, V2 : velocities [m/s] (0 at a large reservoir surface)

    Energy equation solved (the head equation):
      p1/(rho g) + alpha1 V1^2/2g + z1 + h_pump(Q) =
      p2/(rho g) + alpha2 V2^2/2g + z2 + h_turbine(Q) + h_L(Q)

    alpha1, alpha2 : kinetic-energy correction factors for the endpoint
      velocity profiles. alpha = 2.0 for laminar pipe flow, ~1.05 for
      turbulent, 1.0 for a uniform stream (jets, reservoir surfaces).
      Default 1.0 — the usual engineering shortcut, harmless whenever the
      endpoint velocities are zero or the flow is turbulent; set it
      explicitly when an endpoint is a laminar pipe discharge and the
      velocity head actually matters.

    Usage
    -----
    >>> sys = PipeSystem(fluid=water(20), z1=0, z2=12)
    >>> sys.add(Pipe(L=50, D=0.08, roughness='commercial_steel'))
    >>> sys.add(Pump.from_coeffs(h0=30, a=900))
    >>> Q = sys.solve_flow()
    """

    def __init__(self, fluid, z1=0.0, z2=0.0, p1=0.0, p2=0.0, V1=0.0, V2=0.0,
                 alpha1=1.0, alpha2=1.0, patm=101325.0, name="system"):
        self.fluid = fluid
        self.z1, self.z2 = z1, z2
        self.p1, self.p2 = p1, p2         # GAUGE pressures
        self.V1, self.V2 = V1, V2
        self.alpha1, self.alpha2 = alpha1, alpha2   # KE correction factors
        self.patm = patm                  # absolute atmospheric [Pa]
        self.name = name
        self.components = []

    def add(self, *components):
        for c in components:
            self.components.append(c)
        return self

    # -- bookkeeping -----------------------------------------------------
    @property
    def pipes(self):
        return [c for c in self.components if isinstance(c, (Pipe, Parallel))]

    @property
    def pumps(self):
        return [c for c in self.components if isinstance(c, PUMP_TYPES)]

    @property
    def turbines(self):
        return [c for c in self.components if isinstance(c, Turbine)]

    # -- physics -----------------------------------------------------------
    def head_loss(self, Q):
        """Total friction + minor losses at flow Q [m]."""
        return sum(p.head_loss(Q, self.fluid) for p in self.pipes)

    def static_head(self):
        """Head the system demands at zero flow: elevation + pressure + KE terms."""
        rho_g = self.fluid.rho * G
        return ((self.p2 - self.p1) / rho_g
                + (self.alpha2 * self.V2**2 - self.alpha1 * self.V1**2) / (2 * G)
                + (self.z2 - self.z1))

    def required_head(self, Q):
        """System curve: head that must be supplied to move flow Q."""
        return self.static_head() + self.head_loss(Q) + \
            sum(t.head(Q) for t in self.turbines)

    def supplied_head(self, Q):
        """Head added by all pumps at flow Q."""
        return sum(p.head(Q) for p in self.pumps)

    def energy_residual(self, Q):
        """supplied - required; zero at the operating point."""
        return self.supplied_head(Q) - self.required_head(Q)

    def solve_flow(self, Q_max=None):
        """
        Find the steady operating flow rate [m^3/s] where supplied head
        equals required head (the pump/system curve intersection, or the
        gravity-driven balance if there is no pump).
        """
        if self.energy_residual(1e-9) <= 0:
            raise RuntimeError(
                "No flow: at Q ~ 0 the system demands more head than is "
                "supplied. Check elevations, pressures, and pump size.")
        if Q_max is None:
            Q_max = 1e-6
            while self.energy_residual(Q_max) > 0:
                Q_max *= 2.0
                if Q_max > 1e4:
                    raise RuntimeError("Could not bracket a solution; "
                                       "does the pump curve ever fall below "
                                       "the system curve?")
        return brentq(self.energy_residual, 1e-9, Q_max, xtol=1e-12)

    def solve_head_required(self, Q):
        """Head a pump would need to supply to drive flow Q (ignores pumps)."""
        return self.static_head() + self.head_loss(Q) + \
            sum(t.head(Q) for t in self.turbines)

    # -- cavitation --------------------------------------------------------
    def npsh_available(self, Q):
        """
        Net Positive Suction Head available at the (first) pump inlet [m]:

          NPSH_a = p_abs,1/(rho g) + V1^2/(2g) + (z1 - z_pump)
                   - h_L(suction side) - p_v/(rho g)

        The suction side is every pipe added BEFORE the pump. The pump
        elevation is pump.z (defaults to z1 if not set -- flooded suction
        at source level). Requires the fluid to have a vapor pressure.
        """
        if not self.pumps:
            raise RuntimeError("No pump in system; NPSH is undefined.")
        if self.fluid.pv is None:
            raise ValueError(
                f"{self.fluid.name} has no vapor pressure set. Use water(T) "
                "or Fluid(..., pv=...) [Pa absolute] for cavitation checks.")
        pump = self.pumps[0]
        z_pump = pump.z if pump.z is not None else self.z1
        rho_g = self.fluid.rho * G

        h_suction = 0.0
        for c in self.components:
            if c is pump:
                break
            if isinstance(c, (Pipe, Parallel)):
                h_suction += c.head_loss(Q, self.fluid)

        return ((self.p1 + self.patm) / rho_g + self.alpha1 * self.V1**2 / (2 * G)
                + (self.z1 - z_pump) - h_suction - self.fluid.pv / rho_g)

    def check_cavitation(self, Q=None, verbose=True):
        """
        Compare NPSH-available against the pump's NPSH-required at the
        operating point. Returns a dict with npsha, npshr, margin, ok.
        A margin below ~0.5-1 m is considered unsafe practice even if
        positive -- manufacturer curves are measured on clean test rigs.
        """
        if Q is None:
            Q = self.solve_flow()
        pump = self.pumps[0]
        npsha = self.npsh_available(Q)
        npshr = pump.npsh_required(Q)
        result = {"Q": Q, "npsha": npsha, "npshr": npshr,
                  "margin": None if npshr is None else npsha - npshr,
                  "ok": None if npshr is None else npsha > npshr}
        if verbose:
            rho_g = self.fluid.rho * G
            print(f"CAVITATION CHECK at Q = {Q*1000:.2f} L/s "
                  f"({self.fluid.name}, p_v = {self.fluid.pv/1000:.2f} kPa "
                  f"= {self.fluid.pv/rho_g:.2f} m of head)")
            print(f"  NPSH available: {npsha:6.2f} m")
            if npshr is None:
                print("  NPSH required : not specified for this pump "
                      "(pass npsh_r=... to enable the verdict)")
            else:
                print(f"  NPSH required : {npshr:6.2f} m")
                print(f"  margin        : {npsha-npshr:+6.2f} m  ->  "
                      + ("OK" if npsha - npshr > 1.0 else
                         "MARGINAL (< 1 m safety margin)" if npsha > npshr
                         else "*** CAVITATION -- design not viable ***"))
        return result

    # -- reporting -----------------------------------------------------------
    def report(self, Q=None):
        """Solve (if needed) and print a full system summary at flow Q."""
        if Q is None:
            Q = self.solve_flow()
        print("=" * 60)
        print(f"SYSTEM REPORT: {self.name}   ({self.fluid.name})")
        print("=" * 60)
        print(f"Flow rate: Q = {Q:.5f} m^3/s  =  {Q*1000:.2f} L/s  "
              f"=  {Q*15850.3:.1f} gpm")
        print(f"Static head (elev + pressure + KE): {self.static_head():.3f} m")
        print("-" * 60)
        for comp in self.components:
            if isinstance(comp, (Pipe, Parallel)):
                comp.report(Q, self.fluid)
            elif isinstance(comp, PUMP_TYPES):
                print(f"--- {comp.name} (pump) ---")
                print(f"  head added: {comp.head(Q):.3f} m")
                print(f"  shaft power: {comp.power(Q, self.fluid)/1000:.2f} kW")
                if isinstance(comp, PumpsInParallel):
                    for p, Qi in zip(comp.pumps, comp.flow_split(Q)):
                        print(f"      {p.name}: Q = {Qi*1000:.2f} L/s")
            elif isinstance(comp, Turbine):
                print(f"--- {comp.name} (turbine) ---")
                print(f"  head extracted: {comp.head(Q):.3f} m")
                print(f"  power output: {comp.power(Q, self.fluid)/1000:.2f} kW "
                      f"(eta = {comp.eta})")
        print("-" * 60)
        print(f"Total losses:  {self.head_loss(Q):.3f} m")
        print(f"Required head: {self.required_head(Q):.3f} m   "
              f"Supplied head: {self.supplied_head(Q):.3f} m")
        if (self.pumps and self.fluid.pv is not None
                and self.pumps[0].npsh_required(Q) is not None):
            print("-" * 60)
            self.check_cavitation(Q)
        print("=" * 60)
        return Q

    # -- plotting ------------------------------------------------------------
    def plot_curves(self, Q_max=None, n=200, ax=None):
        """Plot system curve (and pump curve + operating point if pumped)."""
        import matplotlib.pyplot as plt
        if Q_max is None:
            try:
                Q_op = self.solve_flow()
                Q_max = 1.6 * Q_op
            except RuntimeError:
                Q_op = None
                Q_max = 0.05
        else:
            try:
                Q_op = self.solve_flow(Q_max=Q_max * 2)
            except RuntimeError:
                Q_op = None

        Q = np.linspace(1e-6, Q_max, n)
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 5))
        ax.plot(Q * 1000, [self.required_head(q) for q in Q],
                lw=2, label="system curve (required head)")
        if self.pumps:
            ax.plot(Q * 1000, [self.supplied_head(q) for q in Q],
                    lw=2, label="pump curve (supplied head)")
        if Q_op is not None and self.pumps:
            h_op = self.required_head(Q_op)
            ax.plot(Q_op * 1000, h_op, "ko", ms=9, zorder=5)
            ax.annotate(f"  operating point\n  Q = {Q_op*1000:.2f} L/s\n"
                        f"  h = {h_op:.2f} m",
                        (Q_op * 1000, h_op), fontsize=9)
        ax.set_xlabel("Q  [L/s]")
        ax.set_ylabel("head  [m]")
        ax.set_title(self.name)
        ax.grid(alpha=0.3)
        ax.legend()
        return ax


# ----------------------------------------------------------------------------
# Unit helpers (conversions are the STUDENT'S responsibility at the boundary)
# ----------------------------------------------------------------------------

def gpm_to_m3s(gpm):   return gpm / 15850.3
def m3s_to_gpm(m3s):   return m3s * 15850.3
def Ls_to_m3s(Ls):     return Ls / 1000.0
def inch_to_m(inch):   return inch * 0.0254
def ft_to_m(ft):       return ft * 0.3048
def psi_to_Pa(psi):    return psi * 6894.76
def hp_to_W(hp):       return hp * 745.7
def W_to_hp(W):        return W / 745.7


# ----------------------------------------------------------------------
# The course pump catalog — eight machines spanning the specific-speed map.
# Used by the term project ("The Pump Station"); prices are teaching props.
# Curves: H = h0 - a*Q^2  [m, m^3/s]; NPSH_r = n0 + n1*Q^2  [m].
# ----------------------------------------------------------------------

PUMP_CATALOG_SPECS = {
    #  name              h0     a       eta   n0    n1      rpm    price
    'GN-25  "Gopher"':   (18,  2.0e5,  0.58, 1.0,  1.5e4,  1750,  1400),
    'GN-40  "Badger"':   (28,  1.1e5,  0.66, 1.2,  1.2e4,  1750,  2300),
    'CV-50  "Coyote"':   (38,  7.0e4,  0.78, 1.2,  1.0e4,  1750,  3200),
    'CV-80  "Elk"':      (52,  3.5e4,  0.74, 1.8,  6.5e3,  1750,  5600),
    'HS-120 "Bighorn"':  (85,  2.8e4,  0.70, 2.6,  3.4e3,  3500,  9800),
    'MF-200 "Moose"':    (16,  6.0e1,  0.80, 2.0,  8.0,    960,  11500),
    'AX-350 "Pelican"':  (7,   4.8,    0.76, 1.5,  1.2,    705,  14800),
    'AX-600 "Heron"':    (5,   0.85,   0.72, 1.4,  0.4,    590,  21000),
}


def pump_catalog():
    """The course pump catalog as a dict of ready-to-use Pump objects.

    Each entry also carries .price, .rpm and .name attributes. Radial
    machines top the list (high head, modest flow); the axials at the
    bottom move rivers over fences. Specific speed sorts them -- which is
    the point.
    """
    out = {}
    for name, (h0, a, eta, n0, n1, rpm, price) in PUMP_CATALOG_SPECS.items():
        p = Pump.from_coeffs(h0, a, efficiency=eta,
                             npsh_r=(lambda Q, n0=n0, n1=n1: n0 + n1 * Q**2))
        p.name, p.rpm, p.price = name, rpm, price
        out[name] = p
    return out


def catalog_report():
    """Print the pump catalog the way a supplier's one-pager would."""
    print(f"{'model':22s} {'shutoff':>8s} {'BEP-ish Q':>10s} {'eta':>5s} "
          f"{'rpm':>5s} {'price':>8s}")
    for name, (h0, a, eta, n0, n1, rpm, price) in PUMP_CATALOG_SPECS.items():
        q_bep = (h0 / (3 * a)) ** 0.5      # rough: 2/3 of shutoff head
        print(f"{name:22s} {h0:6.0f} m {q_bep*1000:8.1f} L/s {eta:5.2f} "
              f"{rpm:5d} ${price:7,d}")
