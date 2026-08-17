"""
fluidtools.hydrostat — pressure fields, forces on submerged surfaces,
buoyancy, and floating stability.

Industry analog: naval-architecture stability suites (GHS, Maxsurf) and the
gate/dam force calculations in every water-resources design office.

Conventions: SI units; depths h measured positive DOWN from the free surface;
gauge pressures unless noted.
"""
import numpy as np

G = 9.81


def pressure(h, rho=998.2, p_surface=0.0):
    """Gauge pressure [Pa] at depth h [m] below a free surface."""
    return p_surface + rho * G * h


# ----------------------------------------------------------------------------
# Plane surfaces
# ----------------------------------------------------------------------------

class PlaneSurface:
    """
    A flat surface submerged in a liquid, possibly inclined.

    The classic results implemented here:
        F_R   = rho g h_c A                    (resultant force)
        y_cp  = y_c + I_xc / (y_c A)           (center of pressure, along the
                                                surface from the free-surface
                                                line; always BELOW the centroid)

    Parameters
    ----------
    shape : 'rectangle' (width b, height L along the plate) or
            'circle' (diameter d).
    b, L, d : geometry [m]; give (b, L) for a rectangle, d for a circle.
    y_top : slant distance from the free surface to the TOP edge of the
            surface, measured along the plate [m].
    theta_deg : angle of the plate from the horizontal (90 = vertical).
    rho : liquid density [kg/m^3].
    """

    def __init__(self, shape="rectangle", b=None, L=None, d=None,
                 y_top=0.0, theta_deg=90.0, rho=998.2):
        self.shape = shape
        self.theta = np.radians(theta_deg)
        self.rho = rho
        self.y_top = float(y_top)
        if shape == "rectangle":
            if b is None or L is None:
                raise ValueError("rectangle needs width b and length L")
            self.b, self.L = float(b), float(L)
            self.A = self.b * self.L
            self.y_c = self.y_top + self.L / 2.0        # centroid, slant dist
            self.I_xc = self.b * self.L**3 / 12.0
        elif shape == "circle":
            if d is None:
                raise ValueError("circle needs diameter d")
            self.d = float(d)
            self.A = np.pi * self.d**2 / 4.0
            self.y_c = self.y_top + self.d / 2.0
            self.I_xc = np.pi * self.d**4 / 64.0
        else:
            raise ValueError("shape must be 'rectangle' or 'circle'")

    @property
    def h_c(self):
        """Vertical depth of the centroid below the surface [m]."""
        return self.y_c * np.sin(self.theta)

    def force(self):
        """Resultant hydrostatic force magnitude [N], acting normal to plate."""
        return self.rho * G * self.h_c * self.A

    def center_of_pressure(self):
        """Slant distance from the free-surface line to the c.p. [m]."""
        return self.y_c + self.I_xc / (self.y_c * self.A)

    def cp_below_centroid(self):
        """How far below the centroid the resultant acts, along the plate [m]."""
        return self.center_of_pressure() - self.y_c

    def report(self):
        print(f"PlaneSurface ({self.shape}), theta = {np.degrees(self.theta):.0f} deg")
        print(f"  area A          = {self.A:.4f} m^2")
        print(f"  centroid depth  = {self.h_c:.3f} m")
        print(f"  resultant F     = {self.force()/1000:.2f} kN")
        print(f"  center of pressure: {self.center_of_pressure():.4f} m along "
              f"plate ({self.cp_below_centroid()*1000:.1f} mm below centroid)")
        return self.force()


def gate_hinge_reaction(surface: PlaneSurface, hinge_at="top"):
    """
    Force a stop/latch must resist for a gate hinged along one edge [N].
    Moment balance about the hinge; hinge_at in {'top','bottom'}.
    Returns (F_stop, F_resultant).
    """
    F = surface.force()
    y_cp_rel_top = surface.center_of_pressure() - surface.y_top
    Ldim = surface.L if surface.shape == "rectangle" else surface.d
    if hinge_at == "top":
        F_stop = F * y_cp_rel_top / Ldim            # stop at bottom edge
    elif hinge_at == "bottom":
        F_stop = F * (Ldim - y_cp_rel_top) / Ldim   # stop at top edge
    else:
        raise ValueError("hinge_at must be 'top' or 'bottom'")
    return F_stop, F


# ----------------------------------------------------------------------------
# Curved surfaces (via horizontal/vertical decomposition)
# ----------------------------------------------------------------------------

def curved_surface_force(h_proj_top, height, width, volume_above,
                         rho=998.2):
    """
    Force components on a curved surface using the projection method.

    F_H = force on the VERTICAL PROJECTION of the surface (a plane problem):
          rho g h_c A_proj, with h_c the centroid depth of the projection.
    F_V = weight of the liquid (real or imaginary) VERTICALLY ABOVE the
          surface, up to the free surface: rho g * volume_above.

    Parameters
    ----------
    h_proj_top : depth of the top of the vertical projection [m]
    height, width : dimensions of the vertical projection [m]
    volume_above : liquid volume vertically above the surface [m^3]
                   (for a surface holding liquid below/behind it, this is the
                   'imaginary' displaced volume and F_V acts upward)

    Returns dict with F_H, F_V, magnitude, and angle from horizontal [deg].
    """
    A_proj = height * width
    h_c = h_proj_top + height / 2.0
    F_H = rho * G * h_c * A_proj
    F_V = rho * G * volume_above
    return {"F_H": F_H, "F_V": F_V,
            "magnitude": float(np.hypot(F_H, F_V)),
            "angle_deg": float(np.degrees(np.arctan2(F_V, F_H)))}


# ----------------------------------------------------------------------------
# Buoyancy & floating stability
# ----------------------------------------------------------------------------

def buoyant_force(volume_displaced, rho=998.2):
    """Archimedes: F_B = rho g V_displaced [N]."""
    return rho * G * volume_displaced


class Barge:
    """
    A box-section floating body: the canonical stability teaching case and a
    fair first model of a real barge.

    Parameters
    ----------
    length, beam : plan dimensions [m] (beam = width, the stability-critical one)
    depth : hull depth [m] (bottom to deck)
    weight : total weight [N]
    KG : height of the center of gravity above the keel (bottom) [m]
    rho : water density [kg/m^3]

    Key results
    -----------
    draft T   : from W = rho g (L * B * T)
    KB = T/2  : center of buoyancy (box hull)
    BM = I/V  : metacentric radius, I = L B^3 / 12 (waterplane inertia)
    GM = KB + BM - KG : metacentric height. GM > 0 -> stable.
    Righting moment (small angles): W * GM * sin(phi).
    """

    def __init__(self, length, beam, depth, weight, KG, rho=998.2):
        self.L, self.B, self.D = float(length), float(beam), float(depth)
        self.W, self.KG, self.rho = float(weight), float(KG), float(rho)
        if self.draft > self.D:
            raise ValueError(
                f"Barge sinks: required draft {self.draft:.2f} m exceeds "
                f"hull depth {self.D:.2f} m.")

    @property
    def draft(self):
        return self.W / (self.rho * G * self.L * self.B)

    @property
    def KB(self):
        return self.draft / 2.0

    @property
    def BM(self):
        I = self.L * self.B**3 / 12.0
        V = self.L * self.B * self.draft
        return I / V

    @property
    def GM(self):
        return self.KB + self.BM - self.KG

    def is_stable(self):
        return self.GM > 0

    def righting_moment(self, heel_deg):
        """Small-angle righting (+) or capsizing (-) moment [N*m]."""
        return self.W * self.GM * np.sin(np.radians(heel_deg))

    def max_KG(self):
        """The KG at which the barge becomes neutrally stable [m]."""
        return self.KB + self.BM

    def report(self):
        print(f"Barge {self.L} x {self.B} x {self.D} m, W = {self.W/1000:.1f} kN")
        print(f"  draft T = {self.draft:.3f} m   freeboard = {self.D-self.draft:.3f} m")
        print(f"  KB = {self.KB:.3f} m   BM = {self.BM:.3f} m   KG = {self.KG:.3f} m")
        verdict = "STABLE" if self.is_stable() else "*** UNSTABLE - capsizes ***"
        print(f"  GM = {self.GM:+.3f} m   ->  {verdict}")
        print(f"  (KG could rise to {self.max_KG():.3f} m before instability)")
        return self.GM


# Common liquids for statics problems (density, kg/m^3, ~20 C)
LIQUIDS = {
    "water": 998.2, "seawater": 1025.0, "mercury": 13550.0,
    "gasoline": 737.0, "sae30_oil": 891.0, "glycerin": 1260.0,
    "ethanol": 789.0,
}
