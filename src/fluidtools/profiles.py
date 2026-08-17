"""
fluidtools.profiles — viscous flow profiles as objects to READ, not
equations to grind.

The engineering skill this module serves: given a velocity distribution
(from theory, CFD, or PIV data), read off the mechanics —
    slope      -> shear stress      (tau = mu du/dy)
    curvature  -> pressure gradient (mu d2u/dy2 = dp/dx, unidirectional)
    zero slope -> zero shear        (free surfaces, symmetry, velocity max)

Industry analog: post-processing CFD and PIV results, which is how most
practicing engineers actually encounter the Navier-Stokes equations.

Geometry: channel of height h, y measured from the STATIONARY lower wall;
upper wall moves at U. dpdx > 0 is an ADVERSE pressure gradient (pushing
against the wall motion).
"""
import numpy as np

G = 9.81


# ------------------------------------------------ Couette-Poiseuille channel
def couette_poiseuille(y, h, U=0.0, dpdx=0.0, mu=1.0e-3):
    """Velocity u(y) [m/s]: u = U y/h + (dpdx/2mu)(y^2 - h y)."""
    y = np.asarray(y, dtype=float)
    return U * y / h + dpdx / (2 * mu) * (y**2 - h * y)


def shear_couette_poiseuille(y, h, U=0.0, dpdx=0.0, mu=1.0e-3):
    """Shear stress tau(y) [Pa] = mu du/dy = mu U/h + (dpdx/2)(2y - h)."""
    y = np.asarray(y, dtype=float)
    return mu * U / h + dpdx / 2 * (2 * y - h)


def flow_rate_couette_poiseuille(h, U=0.0, dpdx=0.0, mu=1.0e-3):
    """Flow per unit width [m^2/s]: q = U h/2 - dpdx h^3/(12 mu)."""
    return U * h / 2 - dpdx * h**3 / (12 * mu)


def backflow_threshold(h, U, mu):
    """
    The adverse pressure gradient [Pa/m] at which reversed flow first
    appears at the stationary wall: dpdx_crit = 2 mu U / h^2.
    (Set du/dy = 0 at y = 0.) Above this, part of the channel flows
    backward even though the upper wall still drags fluid forward.
    """
    return 2 * mu * U / h**2


def has_backflow(h, U, dpdx, mu):
    """True if the profile contains negative velocity anywhere."""
    if U <= 0:
        return dpdx > 0 or (U < 0)
    return dpdx > backflow_threshold(h, U, mu)


# ------------------------------------------------ Poiseuille pipe flow
def pipe_poiseuille(r, R, dpdx, mu):
    """Hagen-Poiseuille pipe profile u(r) = (-dpdx/4mu)(R^2 - r^2)."""
    r = np.asarray(r, dtype=float)
    return -dpdx / (4 * mu) * (R**2 - r**2)


def pipe_poiseuille_Q(R, dpdx, mu):
    """Pipe flow rate [m^3/s]: Q = -dpdx * pi R^4 / (8 mu)."""
    return -dpdx * np.pi * R**4 / (8 * mu)


# ------------------------------------------------ gravity-driven film
def film_on_wall(y, h, mu, rho=998.2, angle_deg=90.0):
    """
    Liquid film of thickness h draining down a wall inclined angle_deg from
    horizontal; y from the wall. u = (rho g sin(a)/mu)(h y - y^2/2).
    Zero shear (zero slope) at the free surface y = h — read it off the plot!
    """
    y = np.asarray(y, dtype=float)
    a = np.radians(angle_deg)
    return rho * G * np.sin(a) / mu * (h * y - y**2 / 2)


def film_flow_rate(h, mu, rho=998.2, angle_deg=90.0):
    """Film flow per unit width [m^2/s]: q = rho g sin(a) h^3 / (3 mu)."""
    return rho * G * np.sin(np.radians(angle_deg)) * h**3 / (3 * mu)


# ------------------------------------------------ reading measured profiles
def wall_shear_from_data(y, u, mu, deg=2):
    """
    Wall shear stress [Pa] from measured near-wall velocity data (PIV,
    hot-wire traverse, CFD samples): fit a polynomial u(y) of degree `deg`
    and evaluate mu du/dy at y = 0.

    Use points close to the wall (within the region where a low-order
    polynomial is honest). Returns (tau_wall, fitted poly1d).
    """
    y = np.asarray(y, dtype=float)
    u = np.asarray(u, dtype=float)
    p = np.polyfit(y, u, deg)
    dudy0 = p[-2]                      # derivative of the fit at y = 0
    return mu * dudy0, np.poly1d(p)
