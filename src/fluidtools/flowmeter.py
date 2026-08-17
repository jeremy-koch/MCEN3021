"""
fluidtools.flowmeter — differential-pressure flow meters (orifice, venturi,
flow nozzle) and pitot tubes.

Industry analog: ISO-5167-style DP meter sizing, which every process and
instrumentation engineer touches. This is Bernoulli's most employable
application: the meter equation is Bernoulli between the pipe and the
throat, patched with a discharge coefficient for what Bernoulli ignores.

    Q = Cd * A_throat * sqrt( 2 dP / rho ) / sqrt(1 - beta^4)

with beta = d_throat / D_pipe. The sqrt(1 - beta^4) factor is the
"velocity of approach" correction — Bernoulli's V1 term, kept honest.
"""
import numpy as np
from scipy.optimize import brentq

# Typical discharge coefficients (turbulent, well-installed)
CD_DEFAULT = {"orifice": 0.61, "venturi": 0.98, "nozzle": 0.96}


def _check(meter):
    if meter not in CD_DEFAULT:
        raise ValueError(f"meter must be one of {sorted(CD_DEFAULT)}")


def meter_flow(meter, D, d, dP, rho, Cd=None):
    """
    Flow rate [m^3/s] through a DP meter from a measured pressure drop.

    D : pipe inside diameter [m];  d : throat/bore diameter [m]
    dP : pressure difference between upstream tap and throat [Pa]
    """
    _check(meter)
    Cd = CD_DEFAULT[meter] if Cd is None else Cd
    beta = d / D
    if not (0 < beta < 1):
        raise ValueError("need 0 < d < D")
    A_t = np.pi * d**2 / 4
    return Cd * A_t * np.sqrt(2 * dP / rho) / np.sqrt(1 - beta**4)


def meter_dP(meter, D, d, Q, rho, Cd=None):
    """Pressure drop [Pa] a DP meter reads at flow Q — inverse of meter_flow."""
    _check(meter)
    Cd = CD_DEFAULT[meter] if Cd is None else Cd
    beta = d / D
    A_t = np.pi * d**2 / 4
    return rho / 2 * (Q * np.sqrt(1 - beta**4) / (Cd * A_t))**2


def size_meter(meter, D, Q, dP_target, rho, Cd=None):
    """
    Bore diameter d [m] that produces dP_target at design flow Q.
    The instrument engineer's daily problem: match the meter to the
    transmitter's range.
    """
    _check(meter)
    f = lambda d: meter_dP(meter, D, d, Q, rho, Cd) - dP_target
    return brentq(f, 1e-6, D * 0.999, xtol=1e-10)


def permanent_loss(meter, dP, beta):
    """
    UNRECOVERED pressure loss [Pa] — what the meter costs the system forever
    (feed this into a pipeflow analysis as an equivalent loss!).

    Orifice: most of dP is lost (~ (1 - beta^1.9) fraction, a standard
    approximation). Venturi: the diffuser recovers most; ~10-15%% lost
    (we use 0.12). Nozzle: intermediate (~0.45).
    """
    _check(meter)
    if meter == "orifice":
        return dP * (1 - beta**1.9)
    if meter == "venturi":
        return dP * 0.12
    return dP * 0.45


def pitot_velocity(dP, rho):
    """Point velocity [m/s] from a pitot-static reading: V = sqrt(2 dP/rho)."""
    return np.sqrt(2 * dP / rho)


def meter_report(meter, D, d, Q, rho, Cd=None):
    """Full sizing summary for a design flow."""
    beta = d / D
    dP = meter_dP(meter, D, d, Q, rho, Cd)
    loss = permanent_loss(meter, dP, beta)
    V_pipe = Q / (np.pi * D**2 / 4)
    print(f"{meter.upper()} meter: D = {D*1000:.0f} mm, d = {d*1000:.1f} mm "
          f"(beta = {beta:.3f})")
    print(f"  design flow: {Q*1000:.2f} L/s (pipe velocity {V_pipe:.2f} m/s)")
    print(f"  reading at design flow: dP = {dP/1000:.2f} kPa")
    print(f"  permanent loss: {loss/1000:.2f} kPa "
          f"({100*loss/dP:.0f}% of the reading)")
    return {"beta": beta, "dP": dP, "permanent_loss": loss}
