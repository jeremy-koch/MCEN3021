"""
fluidtools.extflow — flat-plate boundary layers, drag estimation and
build-up, and terminal velocity.

Industry analog: preliminary drag bookkeeping in vehicle and aircraft design
(the spirit of DATCOM): sum the drag of components, each estimated from a
reference area and a coefficient, before anyone runs CFD.

Correlations (standard smooth-plate results, Re_x transition ~ 5e5):
  laminar  : cf = 0.664/sqrt(Re_x),  CF = 1.328/sqrt(Re_L), delta = 5.0 x/sqrt(Re_x)
  turbulent: cf = 0.0592 Re_x^-1/5,  CF = 0.074 Re_L^-1/5,  delta = 0.38 x Re_x^-1/5
  mixed    : CF = 0.074 Re_L^-1/5 - 1742/Re_L
"""
import numpy as np
from scipy.optimize import brentq

G = 9.81
RE_TRANSITION = 5e5


# ----------------------------------------------------------------------------
# Flat-plate boundary layer
# ----------------------------------------------------------------------------

def delta(x, U, nu, transition=RE_TRANSITION):
    """Boundary-layer thickness [m] at distance x from the leading edge."""
    Rex = U * x / nu
    if Rex <= 0:
        return 0.0
    if Rex < transition:
        return 5.0 * x / np.sqrt(Rex)
    return 0.38 * x / Rex**0.2


def cf_local(Rex, transition=RE_TRANSITION):
    """Local skin-friction coefficient at Re_x."""
    if Rex <= 0:
        raise ValueError("Re_x must be positive")
    if Rex < transition:
        return 0.664 / np.sqrt(Rex)
    return 0.0592 / Rex**0.2


def CF_flat_plate(ReL, regime="auto", transition=RE_TRANSITION):
    """
    Average friction-drag coefficient for one side of a flat plate.

    regime: 'laminar', 'turbulent' (tripped from the leading edge), or
    'auto' (laminar below the transition Re; mixed correlation above).
    """
    if ReL <= 0:
        raise ValueError("Re_L must be positive")
    if regime == "laminar" or (regime == "auto" and ReL < transition):
        return 1.328 / np.sqrt(ReL)
    if regime == "turbulent":
        return 0.074 / ReL**0.2
    # mixed: turbulent average minus the credit for the laminar start
    A = transition * (0.074 / transition**0.2 - 1.328 / np.sqrt(transition))
    return 0.074 / ReL**0.2 - A / ReL


def plate_drag(U, L, width, rho, mu, sides=2, regime="auto"):
    """Friction drag [N] on a flat plate (both sides by default)."""
    ReL = rho * U * L / mu
    CF = CF_flat_plate(ReL, regime)
    return sides * CF * 0.5 * rho * U**2 * (L * width)


# ----------------------------------------------------------------------------
# Drag coefficients & build-up
# ----------------------------------------------------------------------------

# Representative Cd values (frontal area unless noted). Textbook-typical.
CD_SHAPES = {
    "sphere":               0.47,   # subcritical; use cd_sphere(Re) for detail
    "hemisphere_open_up":   1.42,   # parachute-like
    "hemisphere_round_up":  0.38,
    "disk_normal":          1.17,
    "square_plate_normal":  1.18,
    "cube_face_on":         1.05,
    "long_cylinder_cross":  1.20,   # subcritical
    "streamlined_body":     0.04,
    "half_car_modern":      0.30,   # sedan, frontal area basis
    "pickup_truck":         0.45,
    "semi_truck":           0.65,   # without fairings
    "semi_truck_faired":    0.50,
    "cyclist_upright":      1.10,
    "cyclist_racing_tuck":  0.88,
    "standing_person":      1.15,
    "parachutist_spread":   1.20,
    "ski_jumper":           1.30,
}


def cd_sphere(Re):
    """
    Sphere drag coefficient vs Reynolds number (smooth sphere).
    Stokes below Re=0.1; Schiller-Naumann to ~1000; Newton plateau 0.44 to
    the drag crisis (~3.5e5); ~0.1 just above it (supercritical).
    """
    if Re <= 0:
        raise ValueError("Re must be positive")
    if Re < 0.1:
        return 24.0 / Re
    if Re < 1000.0:
        return 24.0 / Re * (1.0 + 0.15 * Re**0.687)
    if Re < 3.5e5:
        return 0.44
    return 0.10


def drag_force(Cd, A, rho, V):
    """F_D = Cd * (1/2) rho V^2 A  [N]."""
    return Cd * 0.5 * rho * V**2 * A


def drag_power(Cd, A, rho, V):
    """Power to overcome drag at speed V [W]."""
    return drag_force(Cd, A, rho, V) * V


class DragBuildup:
    """
    Component drag bookkeeping: add parts (name, Cd, area), get totals.

    >>> bike = DragBuildup(rho=1.204)
    >>> bike.add('rider (tuck)', Cd=0.88, A=0.36)
    >>> bike.add('frame+wheels', Cd=1.0, A=0.05)
    >>> bike.report(V=12.0)
    """

    def __init__(self, rho=1.204, name="assembly"):
        self.rho = rho
        self.name = name
        self.parts = []

    def add(self, name, Cd, A):
        self.parts.append((name, float(Cd), float(A)))
        return self

    @property
    def CdA(self):
        """The industry currency: total drag area, sum of Cd*A [m^2]."""
        return sum(Cd * A for (_, Cd, A) in self.parts)

    def force(self, V):
        return 0.5 * self.rho * V**2 * self.CdA

    def power(self, V):
        return self.force(V) * V

    def report(self, V):
        q = 0.5 * self.rho * V**2
        print(f"Drag build-up: {self.name} at V = {V} m/s "
              f"(q = {q:.1f} Pa)")
        for name, Cd, A in self.parts:
            F = Cd * A * q
            print(f"  {name:24s} Cd={Cd:5.2f}  A={A:6.3f} m^2  "
                  f"CdA={Cd*A:6.3f}  F={F:8.1f} N "
                  f"({100*Cd*A/self.CdA:4.1f} %)")
        print(f"  TOTAL: CdA = {self.CdA:.3f} m^2, F = {self.force(V):.1f} N, "
              f"P = {self.power(V)/1000:.2f} kW")
        return self.force(V)


# ----------------------------------------------------------------------------
# Terminal velocity
# ----------------------------------------------------------------------------

def terminal_velocity(weight_net, Cd, A, rho):
    """
    Terminal speed [m/s] for constant Cd: W_net = 0.5 rho V^2 Cd A.
    weight_net = weight minus buoyancy [N].
    """
    if weight_net <= 0:
        raise ValueError("net weight must be positive (buoyancy exceeds weight?)")
    return np.sqrt(2.0 * weight_net / (rho * Cd * A))


def terminal_velocity_sphere(d, rho_p, rho_f, mu):
    """
    Terminal speed of a sphere with Re-dependent Cd (solved iteratively).
    Falls if rho_p > rho_f, rises (bubble) if lighter — speed returned either way.
    """
    A = np.pi * d**2 / 4
    Vol = np.pi * d**3 / 6
    W = abs(rho_p - rho_f) * G * Vol
    def residual(V):
        Re = rho_f * V * d / mu
        return 0.5 * rho_f * V**2 * cd_sphere(Re) * A - W
    return brentq(residual, 1e-12, 1e4, xtol=1e-12)
