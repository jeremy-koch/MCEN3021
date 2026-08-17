"""
fluidtools.turbo — turbomachinery fundamentals: specific speed and machine
selection, Euler turbomachine analysis, and ideal pump characteristics.

Industry analog: manufacturer selection software (every pump vendor has
one), which is at heart a specific-speed lookup. Pairs with
fluidtools.pipeflow, which handles what the machine does *in a system*
(operating points, affinity laws, cavitation).
"""
import numpy as np

G = 9.81


# ----------------------------------------------------------------------------
# Specific speed & machine selection
# ----------------------------------------------------------------------------

def specific_speed(omega, Q, H):
    """
    Dimensionless pump specific speed: N_s = omega sqrt(Q) / (g H)^(3/4),
    omega in rad/s, Q in m^3/s, H in m. Evaluated at the BEST-EFFICIENCY
    point — that's the datasheet convention.
    """
    return omega * np.sqrt(Q) / (G * H)**0.75


def specific_speed_US(rpm, gpm, H_ft):
    """US customary N_s = rpm sqrt(gpm) / H_ft^0.75 (the catalog number)."""
    return rpm * np.sqrt(gpm) / H_ft**0.75


def turbine_specific_speed(omega, P, H, rho=998.2):
    """Turbine form: N_st = omega sqrt(P/rho) / (g H)^(5/4), P in W."""
    return omega * np.sqrt(P / rho) / (G * H)**1.25


def classify_pump(Ns):
    """
    Impeller type suggested by dimensionless specific speed (approximate
    textbook boundaries — real selection charts have soft edges):
      Ns < ~1.5 : radial (centrifugal)  — high head, low flow
      1.5 - 3.5 : mixed flow
      Ns > ~3.5 : axial (propeller)     — low head, high flow
    """
    if Ns < 1.5:
        return "radial (centrifugal)"
    if Ns < 3.5:
        return "mixed flow"
    return "axial (propeller)"


def select_machine(omega, Q, H, verbose=True):
    """One-call selection helper: computes Ns and names the machine type."""
    Ns = specific_speed(omega, Q, H)
    kind = classify_pump(Ns)
    if verbose:
        print(f"N_s = {Ns:.2f}  ->  {kind}")
        print(f"(duty: Q = {Q*1000:.1f} L/s at H = {H:.1f} m, "
              f"{omega*60/(2*np.pi):.0f} rpm)")
    return Ns, kind


def rpm_to_rad(rpm):
    return rpm * 2 * np.pi / 60


# ----------------------------------------------------------------------------
# Euler turbomachine analysis
# ----------------------------------------------------------------------------

def euler_head(u2, vt2, u1=0.0, vt1=0.0):
    """
    Euler head [m]: h = (u2 vt2 - u1 vt1) / g.
    u = blade speed, vt = tangential (swirl) component of absolute velocity.
    Radial inflow with no pre-swirl -> u1 vt1 = 0 (the usual pump case).
    """
    return (u2 * vt2 - u1 * vt1) / G


class IdealPump:
    """
    Ideal (frictionless, infinite-blade) centrifugal pump from impeller
    geometry, via the Euler equation and exit velocity triangle:

        u2   = omega r2
        vr2  = Q / (2 pi r2 b2)            (radial component, from continuity)
        vt2  = u2 - vr2 / tan(beta2)       (backswept blades: beta2 < 90)
        h    = u2 vt2 / g                  (no inlet swirl)

    The ideal h(Q) is a straight line: h = u2^2/g - [u2/(g 2 pi r2 b2 tan b2)] Q.
    A real pump falls below it (slip, friction, incidence) — comparing this
    line to a measured curve is the whole story of pump losses.
    """

    def __init__(self, r2, b2, beta2_deg, rpm):
        self.r2, self.b2 = float(r2), float(b2)
        self.beta2 = np.radians(beta2_deg)
        self.omega = rpm_to_rad(rpm)

    @property
    def u2(self):
        return self.omega * self.r2

    def vr2(self, Q):
        return Q / (2 * np.pi * self.r2 * self.b2)

    def vt2(self, Q):
        return self.u2 - self.vr2(Q) / np.tan(self.beta2)

    def head(self, Q):
        """Ideal head [m] at flow Q (shutoff = u2^2/g at Q = 0)."""
        return euler_head(self.u2, self.vt2(Q))

    def power(self, Q, rho=998.2):
        """Ideal shaft power [W] = rho g Q h_ideal."""
        return rho * G * Q * self.head(Q)

    def shutoff_head(self):
        return self.u2**2 / G

    def report(self, Q):
        print(f"IdealPump: r2={self.r2*1000:.0f} mm, b2={self.b2*1000:.0f} mm, "
              f"beta2={np.degrees(self.beta2):.0f} deg, "
              f"{self.omega*60/(2*np.pi):.0f} rpm")
        print(f"  u2 = {self.u2:.2f} m/s, shutoff head = {self.shutoff_head():.1f} m")
        print(f"  at Q = {Q*1000:.1f} L/s: vr2 = {self.vr2(Q):.2f} m/s, "
              f"vt2 = {self.vt2(Q):.2f} m/s")
        print(f"  ideal head = {self.head(Q):.2f} m, "
              f"ideal power = {self.power(Q)/1000:.2f} kW")
        return self.head(Q)


def hydraulic_efficiency(h_actual, h_euler):
    """eta_h = actual head / Euler head."""
    return h_actual / h_euler
