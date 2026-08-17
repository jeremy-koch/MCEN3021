"""Validation suite for fluidtools' new modules (hydrostat, momentum,
similitude, extflow, flowmeter, channel, turbo). Textbook-value checks."""
import numpy as np
import pytest
from fluidtools import hydrostat as hs
from fluidtools import momentum as mo
from fluidtools import similitude as si
from fluidtools import extflow as ef
from fluidtools import flowmeter as fm
from fluidtools import channel as ch
from fluidtools import turbo as tb

G = 9.81


# ------------------------------------------------------------- hydrostat
def test_vertical_gate_classic():
    # 3 m tall x 2 m wide vertical gate, top at surface, water:
    # F = rho g (H/2)(bH); cp at 2/3 H from surface.
    g = hs.PlaneSurface("rectangle", b=2, L=3, y_top=0, theta_deg=90)
    assert abs(g.force() - 998.2 * G * 1.5 * 6) < 1e-6
    assert abs(g.center_of_pressure() - 2.0) < 1e-12

def test_cp_always_below_centroid():
    for ytop in [0, 1, 5, 20]:
        g = hs.PlaneSurface("rectangle", b=1, L=2, y_top=ytop)
        assert g.cp_below_centroid() > 0
    # ...and it approaches the centroid as the surface goes deep
    shallow = hs.PlaneSurface("rectangle", b=1, L=2, y_top=0).cp_below_centroid()
    deep = hs.PlaneSurface("rectangle", b=1, L=2, y_top=50).cp_below_centroid()
    assert deep < shallow / 10

def test_inclined_gate_uses_vertical_depth():
    # 45-degree plate: centroid depth = slant * sin(45)
    g = hs.PlaneSurface("rectangle", b=1, L=2, y_top=1, theta_deg=45)
    assert abs(g.h_c - 2 * np.sin(np.radians(45))) < 1e-12

def test_hinge_reaction_splits_force():
    g = hs.PlaneSurface("rectangle", b=2, L=3, y_top=0)
    F_stop, F = hs.gate_hinge_reaction(g, hinge_at="top")
    F_stop2, _ = hs.gate_hinge_reaction(g, hinge_at="bottom")
    assert abs((F_stop + F_stop2) - F) < 1e-6      # the two stops carry it all

def test_curved_surface_pythagoras():
    r = hs.curved_surface_force(h_proj_top=0, height=2, width=3,
                                volume_above=4.0)
    assert abs(r["magnitude"] - np.hypot(r["F_H"], r["F_V"])) < 1e-9
    assert abs(r["F_V"] - 998.2 * G * 4.0) < 1e-6

def test_barge_BM_formula():
    # Box barge: BM = B^2/(12 T). W chosen for T = 2 m draft.
    L, B, T = 10, 6, 2
    W = 998.2 * G * L * B * T
    b = hs.Barge(length=L, beam=B, depth=4, weight=W, KG=1.0)
    assert abs(b.draft - T) < 1e-9
    assert abs(b.BM - B**2 / (12 * T)) < 1e-9
    assert abs(b.GM - (T/2 + 1.5 - 1.0)) < 1e-9
    assert b.is_stable()

def test_barge_capsizes_with_high_KG():
    L, B, T = 10, 3, 2
    W = 998.2 * G * L * B * T
    b = hs.Barge(length=L, beam=B, depth=4, weight=W, KG=2.5)
    assert not b.is_stable()          # BM = 9/24 = 0.375; KB=1 -> max KG 1.375
    assert b.righting_moment(5) < 0

def test_barge_sinks_raises():
    with pytest.raises(ValueError):
        hs.Barge(10, 6, depth=1.0, weight=998.2*G*10*6*2, KG=0.5)


# ------------------------------------------------------------- momentum
def test_jet_normal_plate():
    r = mo.jet_on_flat_plate(rho=1000, V=10, A=0.001, plate_angle_deg=90)
    assert abs(r["Fx"] - 1000 * 0.001 * 100) < 1e-9        # rho A V^2

def test_pelton_cup_doubles():
    flat = mo.jet_on_curved_vane(1000, 10, 0.001, turn_angle_deg=90)
    cup = mo.jet_on_curved_vane(1000, 10, 0.001, turn_angle_deg=180)
    assert abs(cup["Fx"] - 2 * flat["Fx"]) < 1e-9

def test_pelton_optimum_speeds():
    # Wheel of buckets peaks at U = V/2; a single vane at U = V/3.
    V, rho, A = 20.0, 1000.0, 0.002
    Us = np.linspace(0.1, V - 0.1, 2000)
    P_wheel = [mo.jet_on_curved_vane(rho, V, A, 180, u, series=True)["power"] for u in Us]
    P_single = [mo.jet_on_curved_vane(rho, V, A, 180, u, series=False)["power"] for u in Us]
    assert abs(Us[np.argmax(P_wheel)] - V/2) < 0.05
    assert abs(Us[np.argmax(P_single)] - V/3) < 0.05

def test_bend_90_momentum_only():
    # Same D, zero gauge pressures, 90-degree bend: |R| = sqrt(2) rho Q V
    rho, Q, D = 1000, 0.02, 0.1
    A = np.pi * D**2 / 4
    V = Q / A
    r = mo.bend_anchor_force(rho, Q, D, D, 0, 0, 90)
    assert abs(r["magnitude"] - np.sqrt(2) * rho * Q * V) < 1e-6

def test_nozzle_thrust():
    # T = rho Q V_e for atmospheric exit
    rho, Q, De = 1000, 0.01, 0.02
    Ve = Q / (np.pi * De**2 / 4)
    assert abs(mo.nozzle_thrust(rho, Q, De) - rho * Q * Ve) < 1e-9


# ------------------------------------------------------------- similitude
def test_pi_count_drag_problem():
    grps = si.pi_groups({'F': 'force', 'V': 'velocity', 'D': 'length',
                         'rho': 'density', 'mu': 'viscosity'},
                        repeating=['rho', 'V', 'D'], verbose=False)
    assert len(grps) == 2               # 5 vars - 3 dims

def test_pi_groups_are_dimensionless():
    varset = {'dP': 'pressure', 'V': 'velocity', 'D': 'length', 'L': 'length',
              'rho': 'density', 'mu': 'viscosity'}
    for grp in si.pi_groups(varset, verbose=False):
        assert si.check_dimensionless(grp, varset)

def test_reynolds_emerges():
    grps = si.pi_groups({'F': 'force', 'V': 'velocity', 'D': 'length',
                         'rho': 'density', 'mu': 'viscosity'},
                        repeating=['rho', 'V', 'D'], verbose=False)
    # one group must involve mu with rho,V,D exponents proportional to
    # (1,1,1) against mu^-1 (i.e. Reynolds or its inverse)
    mu_grp = next(g for g in grps if 'mu' in g)
    s = -mu_grp['mu']
    assert mu_grp.get('rho', 0) == s and mu_grp.get('V', 0) == s \
        and mu_grp.get('D', 0) == s

def test_bad_repeating_set_raises():
    with pytest.raises(ValueError):
        si.pi_groups({'F': 'force', 'V': 'velocity', 'D': 'length',
                      'L2': 'length', 'rho': 'density'},
                     repeating=['D', 'L2', 'V'], verbose=False)  # D,L2 dependent

def test_froude_scaling():
    r = si.model_scale('froude', Lr=1/25)
    assert abs(r['V'] - 1/5) < 1e-12
    assert abs(r['Q'] - (1/5)*(1/25)**2) < 1e-15
    assert abs(r['F'] - (1/25)**3) < 1e-15      # same fluid, same g

def test_reynolds_scaling_same_fluid():
    r = si.model_scale('reynolds', Lr=1/10)
    assert abs(r['V'] - 10) < 1e-12             # model must run 10x faster
    assert abs(r['F'] - 1.0) < 1e-12            # forces match at Re equality


# ------------------------------------------------------------- extflow
def test_blasius_numbers():
    assert abs(ef.CF_flat_plate(1e5, regime='laminar') - 1.328/np.sqrt(1e5)) < 1e-12
    assert abs(ef.cf_local(1e4) - 0.664/100) < 1e-12

def test_mixed_CF_between_pure_regimes():
    ReL = 1e7
    cf_mixed = ef.CF_flat_plate(ReL, 'auto')
    assert ef.CF_flat_plate(ReL, 'laminar') < cf_mixed < ef.CF_flat_plate(ReL, 'turbulent')

def test_delta_grows_and_jumps_at_transition():
    U, nu = 10.0, 1.5e-5
    x_tr = ef.RE_TRANSITION * nu / U
    assert ef.delta(0.9*x_tr, U, nu) < ef.delta(1.5*x_tr, U, nu)

def test_stokes_terminal_velocity():
    # 50-micron sand in water: Stokes V = (rho_p-rho_f) g d^2 / (18 mu)
    d, rho_p, rho_f, mu = 30e-6, 2650, 998.2, 1.002e-3
    V_stokes = (rho_p - rho_f) * G * d**2 / (18 * mu)
    V = ef.terminal_velocity_sphere(d, rho_p, rho_f, mu)
    Re = rho_f * V * d / mu
    assert Re < 0.1                       # confirm Stokes regime applies
    assert abs(V - V_stokes) / V_stokes < 0.01

def test_drag_buildup_sums():
    b = ef.DragBuildup(rho=1.204)
    b.add('a', 1.0, 0.5).add('b', 0.5, 1.0)
    assert abs(b.CdA - 1.0) < 1e-12
    assert abs(b.force(10) - 0.5*1.204*100*1.0) < 1e-9

def test_cd_sphere_regimes():
    assert abs(ef.cd_sphere(0.01) - 2400) < 1e-9          # Stokes 24/Re
    assert abs(ef.cd_sphere(1e4) - 0.44) < 1e-12          # Newton plateau
    assert ef.cd_sphere(1e6) < 0.2                        # post-crisis


# ------------------------------------------------------------- flowmeter
def test_meter_roundtrip():
    Q = fm.meter_flow('orifice', D=0.1, d=0.05, dP=10e3, rho=998.2)
    dP = fm.meter_dP('orifice', D=0.1, d=0.05, Q=Q, rho=998.2)
    assert abs(dP - 10e3) < 1e-6

def test_orifice_hand_value():
    # beta=0.5: Q = 0.61 * A_t * sqrt(2 dP/rho) / sqrt(1-b^4)
    A_t = np.pi * 0.05**2 / 4
    Q_hand = 0.61 * A_t * np.sqrt(2*10e3/998.2) / np.sqrt(1 - 0.5**4)
    assert abs(fm.meter_flow('orifice', 0.1, 0.05, 10e3, 998.2) - Q_hand) < 1e-12

def test_size_meter_inverts():
    d = fm.size_meter('venturi', D=0.15, Q=0.03, dP_target=25e3, rho=998.2)
    assert abs(fm.meter_dP('venturi', 0.15, d, 0.03, 998.2) - 25e3) < 1e-3

def test_venturi_recovers_more_than_orifice():
    lo = fm.permanent_loss('orifice', 10e3, beta=0.5)
    lv = fm.permanent_loss('venturi', 10e3, beta=0.5)
    assert lv < lo / 3

def test_pitot():
    assert abs(fm.pitot_velocity(500, 1.204) - np.sqrt(2*500/1.204)) < 1e-12


# ------------------------------------------------------------- channel
def test_manning_roundtrip():
    c = ch.Channel(b=3.0, n=0.013, S=0.001)
    yn = c.normal_depth(5.0)
    assert abs(c.manning_Q(yn) - 5.0) < 1e-8

def test_critical_depth_rectangular_formula():
    # rect: yc = (q^2/g)^(1/3)
    c = ch.Channel(b=3.0)
    q = 5.0 / 3.0
    assert abs(c.critical_depth(5.0) - (q**2/G)**(1/3)) < 1e-9

def test_froude_at_critical_is_one():
    c = ch.Channel(b=2.0, m=1.5)
    yc = c.critical_depth(8.0)
    assert abs(c.froude(8.0, yc) - 1.0) < 1e-8

def test_jump_conjugate_matches_belanger():
    # rectangular jump: y2/y1 = 0.5 (sqrt(1+8 Fr1^2) - 1)
    c = ch.Channel(b=4.0)
    y1 = 0.3
    Q = 4.0 * y1 * 6.0                    # V1 = 6 m/s
    Fr1 = c.froude(Q, y1)
    y2 = c.conjugate_depth(Q, y1)
    assert abs(y2/y1 - 0.5*(np.sqrt(1 + 8*Fr1**2) - 1)) < 1e-6

def test_jump_loses_energy():
    c = ch.Channel(b=4.0)
    assert c.jump_energy_loss(Q=7.2, y1=0.3) > 0

def test_jump_from_subcritical_raises():
    c = ch.Channel(b=4.0)
    with pytest.raises(ValueError):
        c.conjugate_depth(Q=1.0, y1=2.0)

def test_gvf_m1_backwater_rises_upstream():
    # Mild slope, depth above normal at a dam: M1 curve, deepening downstream
    c = ch.Channel(b=5.0, n=0.015, S=0.0005)
    Q = 8.0
    yn = c.normal_depth(Q)
    x, y = c.gvf_profile(Q, y_start=1.5*yn, x_end=-3000)  # march upstream
    assert y[-1] < y[0]                    # approaches yn going upstream
    assert y[-1] > yn * 0.99

def test_slope_classification():
    steep = ch.Channel(b=3.0, n=0.012, S=0.05)
    mild = ch.Channel(b=3.0, n=0.012, S=0.0005)
    assert steep.classify_slope(5.0) == "steep"
    assert mild.classify_slope(5.0) == "mild"


# ------------------------------------------------------------- turbo
def test_specific_speed_and_classification():
    # High head, low flow -> radial. 3000 rpm, 10 L/s, 80 m:
    Ns = tb.specific_speed(tb.rpm_to_rad(3000), 0.010, 80.0)
    assert tb.classify_pump(Ns) == "radial (centrifugal)"
    # Low head, big flow -> axial. 300 rpm, 3 m^3/s, 2 m:
    Ns2 = tb.specific_speed(tb.rpm_to_rad(300), 3.0, 2.0)
    assert tb.classify_pump(Ns2) == "axial (propeller)"

def test_euler_shutoff_head():
    p = tb.IdealPump(r2=0.10, b2=0.02, beta2_deg=25, rpm=1750)
    assert abs(p.head(0.0) - p.u2**2 / G) < 1e-9

def test_ideal_head_is_linear_decreasing_for_backswept():
    p = tb.IdealPump(r2=0.10, b2=0.02, beta2_deg=25, rpm=1750)
    Qs = np.linspace(0, 0.05, 5)
    hs_ = [p.head(q) for q in Qs]
    d2 = np.diff(hs_, 2)
    assert np.all(np.abs(d2) < 1e-9)      # linear
    assert hs_[0] > hs_[-1]               # decreasing

def test_us_specific_speed():
    # a familiar catalog magnitude: 1750 rpm, 500 gpm, 100 ft -> ~1237
    assert abs(tb.specific_speed_US(1750, 500, 100) - 1237.4) < 1.0


# ---------------------------------------------------- regime identification
def test_named_numbers():
    assert abs(si.reynolds(998.2, 1.0, 0.05, 1.0e-3) - 49910) < 1
    assert abs(si.froude(3.0, 0.918) - 1.0) < 1e-3
    assert abs(si.mach(343.0) - 1.0) < 1e-12
    assert abs(si.weber(998.2, 2.0, 0.002, 0.0728) - 998.2*4*0.002/0.0728) < 1e-9

def test_regime_verdicts():
    assert "laminar" in si.regime('reynolds', 1500, 'pipe')
    assert "turbulent" in si.regime('reynolds', 1e5, 'pipe')
    assert "Stokes" in si.regime('reynolds', 0.5, 'sphere')
    assert "subcritical" in si.regime('froude', 0.6)
    assert "supercritical" in si.regime('froude', 2.0)
    assert "incompressible" in si.regime('mach', 0.1)
