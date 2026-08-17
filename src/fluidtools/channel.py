"""
fluidtools.channel — open-channel flow: uniform flow (Manning), critical
flow, hydraulic jumps, and gradually-varied surface profiles.

Industry analog: HEC-RAS — arguably the most-used fluids software on Earth
(every floodplain study, culvert, and drainage channel). This module is the
1% of HEC-RAS that carries 90% of the concepts.

SI Manning equation:  V = (1/n) R_h^(2/3) S^(1/2)
"""
import numpy as np
from scipy.optimize import brentq

G = 9.81

# Manning roughness n (typical values)
MANNING_N = {
    "glass":              0.010,
    "finished_concrete":  0.012,
    "unfinished_concrete": 0.014,
    "brick":              0.015,
    "clean_earth":        0.022,
    "gravel":             0.025,
    "earth_with_weeds":   0.030,
    "natural_clean":      0.030,
    "natural_stony":      0.045,
    "floodplain_trees":   0.100,
}


class Channel:
    """
    A prismatic open channel: rectangular (m=0) or trapezoidal.

    Parameters
    ----------
    b : bottom width [m]
    m : side slope, HORIZONTAL per unit vertical (0 = rectangular)
    n : Manning roughness (value or MANNING_N key)
    S : bottom slope (dimensionless, e.g. 0.001)
    """

    def __init__(self, b, m=0.0, n=0.013, S=0.001, name="channel"):
        self.b, self.m = float(b), float(m)
        self.n = MANNING_N[n] if isinstance(n, str) else float(n)
        self.S = float(S)
        self.name = name

    # -- geometry as functions of depth y --------------------------------
    def area(self, y):        return (self.b + self.m * y) * y
    def perimeter(self, y):   return self.b + 2 * y * np.sqrt(1 + self.m**2)
    def top_width(self, y):   return self.b + 2 * self.m * y
    def Rh(self, y):          return self.area(y) / self.perimeter(y)

    # -- uniform (normal) flow --------------------------------------------
    def manning_Q(self, y):
        """Discharge [m^3/s] at depth y if flow were uniform on slope S."""
        A = self.area(y)
        return A / self.n * self.Rh(y)**(2/3) * np.sqrt(self.S)

    def normal_depth(self, Q):
        """Depth [m] at which Manning gives discharge Q on this slope."""
        return brentq(lambda y: self.manning_Q(y) - Q, 1e-6, 100.0, xtol=1e-10)

    # -- critical flow -----------------------------------------------------
    def froude(self, Q, y):
        """Fr = V / sqrt(g * A/T) — the hydraulic depth form."""
        A, T = self.area(y), self.top_width(y)
        return (Q / A) / np.sqrt(G * A / T)

    def critical_depth(self, Q):
        """Depth where Fr = 1 (Q^2 T / g A^3 = 1)."""
        f = lambda y: Q**2 * self.top_width(y) / (G * self.area(y)**3) - 1.0
        return brentq(f, 1e-6, 100.0, xtol=1e-10)

    def specific_energy(self, Q, y):
        """E = y + V^2/2g [m]."""
        return y + (Q / self.area(y))**2 / (2 * G)

    def alternate_depth(self, Q, y):
        """The other depth with the same specific energy (sub<->supercritical)."""
        E = self.specific_energy(Q, y)
        yc = self.critical_depth(Q)
        f = lambda yy: self.specific_energy(Q, yy) - E
        if y > yc:      # given subcritical -> find supercritical root
            return brentq(f, 1e-6, yc, xtol=1e-10)
        return brentq(f, yc, 100.0, xtol=1e-10)

    # -- hydraulic jump ----------------------------------------------------
    def momentum_function(self, Q, y):
        """M(y) = A*y_bar + Q^2/(g A): conserved across a jump."""
        A = self.area(y)
        # centroid depth of a trapezoid section below the surface
        y_bar = (y / 6) * (3 * self.b + 2 * (self.top_width(y) - self.b)) / \
                ((self.b + self.top_width(y)) / 2) * ((self.b + self.top_width(y)) / 2) \
                if False else self._centroid_depth(y)
        return A * y_bar + Q**2 / (G * A)

    def _centroid_depth(self, y):
        """Depth of the area centroid below the free surface."""
        # rectangle part + two triangles
        A_rect, A_tri = self.b * y, self.m * y**2
        ybar_rect, ybar_tri = y / 2, y / 3
        return (A_rect * ybar_rect + A_tri * ybar_tri) / (A_rect + A_tri)

    def conjugate_depth(self, Q, y1):
        """Downstream depth of a hydraulic jump from supercritical y1."""
        yc = self.critical_depth(Q)
        if y1 >= yc:
            raise ValueError("y1 must be supercritical (y1 < critical depth) "
                             "for a jump to occur.")
        M1 = self.momentum_function(Q, y1)
        f = lambda y2: self.momentum_function(Q, y2) - M1
        return brentq(f, yc, 1000.0, xtol=1e-10)

    def jump_energy_loss(self, Q, y1):
        """Head destroyed by the jump [m] — where the whitewater energy goes."""
        y2 = self.conjugate_depth(Q, y1)
        return self.specific_energy(Q, y1) - self.specific_energy(Q, y2)

    # -- gradually varied flow ---------------------------------------------
    def friction_slope(self, Q, y):
        """S_f from Manning at the local depth."""
        return (self.n * Q / (self.area(y) * self.Rh(y)**(2/3)))**2

    def gvf_profile(self, Q, y_start, x_end, steps=400):
        """
        Integrate the water-surface profile dy/dx = (S0 - Sf)/(1 - Fr^2)
        from x=0 (depth y_start) toward x_end (may be negative to march
        upstream, the usual direction for subcritical profiles).

        Returns (x array, y array). Integration stops if the depth
        approaches critical (the equation is singular there).
        """
        x = np.linspace(0.0, x_end, steps)
        dx = x[1] - x[0]
        y = np.empty_like(x)
        y[0] = y_start
        yc = self.critical_depth(Q)
        for i in range(1, len(x)):
            Fr2 = self.froude(Q, y[i-1])**2
            if abs(1 - Fr2) < 1e-3:
                return x[:i], y[:i]          # reached critical: stop cleanly
            dydx = (self.S - self.friction_slope(Q, y[i-1])) / (1 - Fr2)
            y[i] = y[i-1] + dydx * dx
            if y[i] <= 1e-4:
                return x[:i], y[:i]
        return x, y

    def classify_slope(self, Q):
        """'mild' (yn > yc), 'steep' (yn < yc), or 'critical'."""
        yn, yc = self.normal_depth(Q), self.critical_depth(Q)
        if abs(yn - yc) / yc < 0.01:
            return "critical"
        return "mild" if yn > yc else "steep"

    def report(self, Q):
        yn, yc = self.normal_depth(Q), self.critical_depth(Q)
        print(f"Channel '{self.name}': b={self.b} m, m={self.m}, "
              f"n={self.n}, S={self.S}")
        print(f"  Q = {Q} m^3/s")
        print(f"  normal depth   yn = {yn:.3f} m  (Fr = {self.froude(Q,yn):.3f})")
        print(f"  critical depth yc = {yc:.3f} m")
        print(f"  slope class: {self.classify_slope(Q).upper()}")
        return yn, yc
