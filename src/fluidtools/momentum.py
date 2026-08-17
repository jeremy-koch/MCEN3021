"""
fluidtools.momentum — control-volume momentum analysis: jet forces on vanes,
anchoring forces on bends and nozzles, and static thrust.

Industry analog: pipe-support and anchor-block load calculations, and
first-cut propulsion sizing.

Sign convention: +x along the incoming flow, +y left of it. All results are
forces ON the vane/fitting FROM the fluid (report both action and reaction
in your write-ups and say which you mean!).
"""
import numpy as np

G = 9.81


def jet_on_flat_plate(rho, V, A, plate_angle_deg=90.0, plate_speed=0.0):
    """
    Force on a flat plate struck by a liquid jet [N].

    plate_angle_deg : angle between jet and plate (90 = normal impact).
    plate_speed : plate velocity in the jet direction (moving-vane problem);
                  the analysis uses the RELATIVE velocity Vr = V - U.

    Returns dict: Fx, Fy (on the plate), and the relative velocity used.
    For normal impact of a stationary plate this reduces to F = rho A V^2.
    """
    Vr = V - plate_speed
    if Vr <= 0:
        return {"Fx": 0.0, "Fy": 0.0, "V_rel": Vr}
    mdot = rho * A * Vr
    th = np.radians(plate_angle_deg)
    # Jet splits along the plate; only the normal momentum component is
    # destroyed. Normal direction force = mdot * Vr * sin(theta).
    Fn = mdot * Vr * np.sin(th)
    return {"Fx": Fn * np.sin(th), "Fy": Fn * np.cos(th), "V_rel": Vr}


def jet_on_curved_vane(rho, V, A, turn_angle_deg, vane_speed=0.0,
                       friction_factor=1.0, series=False):
    """
    Force on a curved vane that turns a jet by turn_angle_deg
    (180 = full reversal, the Pelton-bucket idealization).

    vane_speed : vane velocity in the jet direction.
    friction_factor : exit-to-inlet relative speed ratio (1 = frictionless).
    series : False -> a SINGLE vane running away from the jet; it intercepts
             only mdot = rho A (V - U), and max power occurs at U = V/3.
             True -> a WHEEL of buckets (Pelton turbine); some bucket always
             takes the full jet, mdot = rho A V, and max power is at U = V/2.
             The distinction is a classic exam ambush — know which you have!

    Returns dict: Fx, Fy on the vane, power extracted (Fx * vane_speed),
    and the mass flow used.
    """
    U = vane_speed
    Vr = V - U
    if Vr <= 0:
        return {"Fx": 0.0, "Fy": 0.0, "power": 0.0, "mdot": 0.0}
    mdot = rho * A * (V if series else Vr)
    th = np.radians(turn_angle_deg)
    k = friction_factor
    Fx = mdot * Vr * (1.0 - k * np.cos(th))
    Fy = -mdot * Vr * k * np.sin(th)
    return {"Fx": Fx, "Fy": Fy, "power": Fx * U, "mdot": mdot}


def pelton_optimum(V, series=True):
    """Vane speed for max power: V/2 for a bucket WHEEL (series=True,
    the Pelton turbine result), V/3 for a SINGLE runaway vane."""
    return V / 2.0 if series else V / 3.0


def bend_anchor_force(rho, Q, D1, D2, p1_gauge, p2_gauge, turn_angle_deg,
                      weight=0.0):
    """
    Anchoring force required to hold a reducing pipe bend [N].

    Momentum + pressure balance on a control volume around the bend:
      sum F = mdot (V2_vec - V1_vec)
    with gauge pressures acting on the inlet/outlet areas.

    Inlet flow along +x; outlet turned by turn_angle_deg toward +y.
    weight : total weight of bend + contained fluid (acts -y) [N].

    Returns dict: Rx, Ry = components the ANCHOR must supply, plus magnitude.
    """
    A1, A2 = np.pi * D1**2 / 4, np.pi * D2**2 / 4
    V1, V2 = Q / A1, Q / A2
    th = np.radians(turn_angle_deg)
    mdot = rho * Q
    # x: p1 A1 - p2 A2 cos(th) + Rx = mdot (V2 cos th - V1)
    Rx = mdot * (V2 * np.cos(th) - V1) - p1_gauge * A1 + p2_gauge * A2 * np.cos(th)
    # y: -p2 A2 sin(th) + Ry - weight = mdot (V2 sin th - 0)
    Ry = mdot * V2 * np.sin(th) + p2_gauge * A2 * np.sin(th) + weight
    return {"Rx": Rx, "Ry": Ry, "magnitude": float(np.hypot(Rx, Ry)),
            "V1": V1, "V2": V2}


def nozzle_thrust(rho, Q, D_exit, p_exit_gauge=0.0):
    """
    Static thrust of a nozzle discharging to atmosphere [N]:
    T = mdot * V_e + (p_e - p_atm) A_e. Positive = pushes back on the mount.
    """
    A_e = np.pi * D_exit**2 / 4
    V_e = Q / A_e
    return rho * Q * V_e + p_exit_gauge * A_e


def rocket_thrust(mdot, V_exhaust, p_exit_gauge=0.0, A_exit=0.0):
    """T = mdot V_e + (p_e - p_atm) A_e — same physics, propulsion notation."""
    return mdot * V_exhaust + p_exit_gauge * A_exit


def sluice_gate_force(rho, q, y1, y2):
    """
    Horizontal force per unit width on a sluice gate [N/m], from the
    momentum equation between upstream depth y1 and downstream depth y2
    with discharge per unit width q [m^2/s]. Hydrostatic pressure at 1, 2.
    """
    V1, V2 = q / y1, q / y2
    F = 0.5 * rho * G * (y1**2 - y2**2) - rho * q * (V2 - V1)
    return F
