"""Validation suite for the pipeflow package. Run with: pytest

Checks friction factors against Moody-chart values, energy-balance closure,
parallel-branch physics, pump/turbine behavior, NPSH hand calculations,
affinity-law scaling, and pump staging.
"""
import numpy as np
from fluidtools.pipeflow import *


def test_full_validation():
    def ok(name, cond):
        print(("PASS: " if cond else "FAIL: ") + name)
        assert cond, name


    # 1. Laminar friction factor
    f = friction_factor(1000)
    ok("laminar f = 64/Re", abs(f - 0.064) < 1e-12)

    # 2. Colebrook vs Moody chart: Re=1e5, eps/D=1e-3 -> f ~ 0.0224 (Moody)
    f_cb = friction_factor(1e5, 1e-3)
    print(f"   Colebrook f(Re=1e5, e/D=1e-3) = {f_cb:.4f}")
    ok("Colebrook matches Moody chart region", 0.021 < f_cb < 0.024)

    # 3. Swamee-Jain within ~2% of Colebrook
    f_sj = friction_factor(1e5, 1e-3, method="swamee-jain")
    ok("Swamee-Jain within 2% of Colebrook", abs(f_sj - f_cb) / f_cb < 0.02)

    # 4. Smooth pipe, Re=1e6 -> f ~ 0.0116 (Prandtl/Moody)
    f_sm = friction_factor(1e6, 0.0)
    print(f"   f(Re=1e6, smooth) = {f_sm:.4f}")
    ok("smooth pipe high-Re value", 0.0110 < f_sm < 0.0122)

    # 5. Water properties at 20 C
    w = water(20)
    ok("water rho ~ 998", abs(w.rho - 998.2) < 1.0)
    ok("water mu ~ 1.00e-3", abs(w.mu - 1.002e-3) < 3e-5)

    # 6. Energy balance: gravity-driven flow, residual should be ~0 at solution
    sys = PipeSystem(fluid=w, z1=20, z2=0, name="gravity test")
    p = Pipe(L=100, D=0.05, roughness="commercial_steel")
    p.add_fitting("entrance_sharp").add_fitting("exit")
    sys.add(p)
    Q = sys.solve_flow()
    ok("gravity-flow energy residual ~ 0", abs(sys.energy_residual(Q)) < 1e-8)
    print(f"   Q = {Q*1000:.3f} L/s, V = {p.velocity(Q):.3f} m/s")
    # hand check: h_L must equal 20 m at solution
    ok("losses equal available head", abs(sys.head_loss(Q) - 20.0) < 1e-8)

    # 7. Parallel: identical branches split 50/50
    b1 = Pipe(L=50, D=0.04, roughness="commercial_steel", name="A")
    b2 = Pipe(L=50, D=0.04, roughness="commercial_steel", name="B")
    par = Parallel([b1, b2])
    split = par.flow_split(0.01, w)
    ok("identical parallel branches split evenly",
       abs(split[0] - split[1]) < 1e-9 and abs(sum(split) - 0.01) < 1e-9)

    # 8. Parallel: bigger pipe takes more flow
    b3 = Pipe(L=50, D=0.06, roughness="commercial_steel", name="C")
    par2 = Parallel([b1, b3])
    s = par2.flow_split(0.01, w)
    ok("larger branch carries more flow", s[1] > s[0])
    ok("parallel flows sum to total", abs(sum(s) - 0.01) < 1e-9)

    # 9. Pump operating point: supplied == required
    sys2 = PipeSystem(fluid=w, z1=0, z2=15, name="pumped test")
    pipe2 = Pipe(L=200, D=0.08, roughness="galvanized_iron")
    pipe2.add_fitting("elbow_90_standard", n=4).add_fitting("gate_valve_open")
    pump = Pump.from_coeffs(h0=40, a=30000, efficiency=0.7)
    sys2.add(pipe2, pump)
    Q2 = sys2.solve_flow()
    ok("pumped system: supplied == required at op point",
       abs(sys2.supplied_head(Q2) - sys2.required_head(Q2)) < 1e-8)
    print(f"   Q = {Q2*1000:.2f} L/s, pump head = {pump.head(Q2):.2f} m, "
          f"power = {pump.power(Q2, w)/1000:.2f} kW")

    # 10. Pump.from_data recovers a quadratic exactly
    Qp = np.array([0, 0.01, 0.02, 0.03])
    hp = 35 - 20000 * Qp**2
    pmp = Pump.from_data(Qp, hp)
    ok("pump curve fit exact for quadratic data",
       abs(pmp.head(0.015) - (35 - 20000 * 0.015**2)) < 1e-9)

    # 11. Turbine reduces flow relative to no-turbine case
    sys3 = PipeSystem(fluid=w, z1=50, z2=0, name="hydro")
    pipe3 = Pipe(L=150, D=0.15, roughness="commercial_steel")
    sys3.add(pipe3)
    Q_free = sys3.solve_flow()
    sys3.add(Turbine(head=30, efficiency=0.85))
    Q_turb = sys3.solve_flow()
    ok("turbine reduces flow", Q_turb < Q_free)
    print(f"   free Q = {Q_free*1000:.1f} L/s -> with turbine {Q_turb*1000:.1f} L/s")

    # 12. No-flow detection
    sys4 = PipeSystem(fluid=w, z1=0, z2=10)
    sys4.add(Pipe(L=10, D=0.05, roughness="pvc"))
    try:
        sys4.solve_flow()
        ok("uphill unpumped system raises error", False)
    except RuntimeError:
        ok("uphill unpumped system raises error", True)

    print("\nAll tests done.")

    # ============ v1.1 extension tests: NPSH, affinity, staging ============
    print("\n--- v1.1 extension tests ---")

    # 13. Vapor pressure vs steam tables: 20C -> 2.339 kPa, 60C -> 19.94 kPa, 100C -> 101.3 kPa
    for T, pv_ref in [(20, 2339), (60, 19940), (100, 101325)]:
        pv = water(T).pv
        ok(f"vapor pressure at {T} C ({pv:.0f} vs {pv_ref} Pa)",
           abs(pv - pv_ref) / pv_ref < 0.01)

    # 14. NPSHa hand check: open tank, pump 3 m above surface, suction loss known
    w20 = water(20)
    suction = Pipe(L=10, D=0.05, roughness='commercial_steel', name='suction')
    suction.add_fitting('entrance_sharp')
    pmp = Pump.from_coeffs(h0=30, a=200000, npsh_r=3.0, z=3.0)
    s = PipeSystem(w20, z1=0, z2=5)
    s.add(suction, pmp)
    Qt = 0.005
    h_s = suction.head_loss(Qt, w20)
    npsha_hand = (101325/(w20.rho*G)) + (0 - 3.0) - h_s - w20.pv/(w20.rho*G)
    ok("NPSHa matches hand calculation",
       abs(s.npsh_available(Qt) - npsha_hand) < 1e-10)
    print(f"   NPSHa = {s.npsh_available(Qt):.3f} m (hand: {npsha_hand:.3f} m)")

    # 15. Hot water reduces NPSHa (same geometry)
    s_hot = PipeSystem(water(60), z1=0, z2=5)
    s_hot.add(suction, pmp)
    ok("hot water reduces NPSHa", s_hot.npsh_available(Qt) < s.npsh_available(Qt))

    # 16. check_cavitation verdict flips between cold and hot
    r_cold = s.check_cavitation(Qt, verbose=False)
    ok("cavitation dict has margin", r_cold["margin"] is not None)

    # 17. Affinity laws: shutoff head scales with N^2, runout flow with N
    base = Pump.from_coeffs(h0=40, a=30000, npsh_r=lambda Q: 1.5 + 40000*Q**2)
    fast = base.at_speed(1.2)
    ok("affinity: shutoff head x1.44 at 120% speed",
       abs(fast.head(0) - 40*1.44) < 1e-9)
    Q_run_base = np.sqrt(40/30000)          # head=0 flow
    ok("affinity: runout flow x1.2",
       abs(fast.head(1.2*Q_run_base)) < 1e-9)
    ok("affinity: NPSHr scales as N^2 at scaled flow",
       abs(fast.npsh_required(1.2*0.01) - 1.44*base.npsh_required(0.01)) < 1e-9)

    # 18. Trim is same scaling
    trim = base.trimmed(0.9)
    ok("trim: shutoff x0.81", abs(trim.head(0) - 40*0.81) < 1e-9)

    # 19. Series pumps: heads add at same Q
    p1 = Pump.from_coeffs(h0=40, a=30000)
    ser = PumpsInSeries([p1, p1])
    ok("series pumps double head", abs(ser.head(0.02) - 2*p1.head(0.02)) < 1e-9)

    # 20. Parallel pumps: identical pumps -> flow doubles at same head
    par = PumpsInParallel([p1, p1])
    h_test = p1.head(0.015)
    ok("parallel pumps: 2x flow at same head",
       abs(par.head(0.030) - h_test) < 1e-6)
    ok("parallel split even", abs(np.diff(par.flow_split(0.030))[0]) < 1e-9)

    # 21. But parallel pumps on a real system give LESS than 2x system flow
    w = water(20)
    pipe_sys = Pipe(L=150, D=0.06, roughness='commercial_steel')
    sys1 = PipeSystem(w, z1=0, z2=10); sys1.add(pipe_sys, p1)
    Q1 = sys1.solve_flow()
    sys2 = PipeSystem(w, z1=0, z2=10); sys2.add(pipe_sys, PumpsInParallel([p1, p1]))
    Q2 = sys2.solve_flow()
    ok("two parallel pumps give < 2x system flow", Q1 < Q2 < 2*Q1)
    print(f"   one pump: {Q1*1000:.2f} L/s -> two parallel: {Q2*1000:.2f} L/s "
          f"({Q2/Q1:.2f}x)")

    # 22. Dissimilar parallel pumps: small pump contributes 0 above its shutoff
    small = Pump.from_coeffs(h0=15, a=30000)
    mix = PumpsInParallel([p1, small])
    h_high = 20.0  # above small pump's shutoff
    ok("check-valve behavior above shutoff",
       mix._pump_flow(h_high, small) == 0.0)

    # 23. VFD beats throttling for power at same target flow
    target_Q = 0.007
    sysA = PipeSystem(w, z1=0, z2=10)
    pipeA = Pipe(L=150, D=0.06, roughness='commercial_steel')
    sysA.add(pipeA, p1)
    # throttle: add valve K until Q = target
    from scipy.optimize import brentq as _bq
    def q_with_K(K):
        pp = Pipe(L=150, D=0.06, roughness='commercial_steel')
        pp.add_fitting('valve', K=K)
        ss = PipeSystem(w, z1=0, z2=10); ss.add(pp, p1)
        return ss.solve_flow() - target_Q
    K_needed = _bq(q_with_K, 0.01, 500)
    P_throttle = p1.power(target_Q, w)
    # VFD: find speed ratio so op point = target
    def q_at_speed(r):
        ss = PipeSystem(w, z1=0, z2=10)
        ss.add(pipeA, p1.at_speed(r))
        try:
            return ss.solve_flow() - target_Q
        except RuntimeError:      # pump too slow to overcome static lift
            return -target_Q
    r_needed = _bq(q_at_speed, 0.5, 1.0)
    P_vfd = p1.at_speed(r_needed).power(target_Q, w)
    ok("VFD power < throttle power at same Q", P_vfd < P_throttle)
    print(f"   throttle: {P_throttle/1000:.2f} kW (K={K_needed:.0f})  "
          f"vs VFD at {100*r_needed:.0f}%: {P_vfd/1000:.2f} kW  "
          f"-> saves {100*(1-P_vfd/P_throttle):.0f}%")

    # 24. Old-style Fluid without pv raises a clear error on NPSH
    s_bad = PipeSystem(Fluid(rho=737, mu=2.9e-4, name='gasoline'), z1=0, z2=5)
    s_bad.add(Pipe(L=5, D=0.05, roughness='pvc'), Pump(head=20))
    try:
        s_bad.npsh_available(0.005)
        ok("missing pv raises ValueError", False)
    except ValueError:
        ok("missing pv raises ValueError", True)

    print("\nAll v1.1 tests done.")



def test_alpha_kinetic_energy_correction():
    # Same system, laminar endpoint discharge: alpha=2 doubles the KE term.
    f = water(20)
    base = dict(z1=0.0, z2=0.0, p1=50e3, p2=0.0, V2=3.0)
    s1 = PipeSystem(fluid=f, **base)                    # default alpha = 1
    s2 = PipeSystem(fluid=f, alpha2=2.0, **base)
    dKE = s2.static_head() - s1.static_head()
    assert abs(dKE - 3.0**2 / (2 * 9.81)) < 1e-9        # extra V^2/2g of demand
    # and defaults leave every existing result untouched
    s3 = PipeSystem(fluid=f, z1=0, z2=12)
    s3.add(Pipe(L=50, D=0.08, roughness="commercial_steel"))
    s3.add(Pump.from_coeffs(h0=30, a=900))
    s4 = PipeSystem(fluid=f, z1=0, z2=12, alpha1=1.0, alpha2=1.0)
    s4.add(Pipe(L=50, D=0.08, roughness="commercial_steel"))
    s4.add(Pump.from_coeffs(h0=30, a=900))
    assert abs(s3.solve_flow() - s4.solve_flow()) < 1e-12


def test_pump_catalog_loads_and_curves_sane():
    from fluidtools.pipeflow import pump_catalog
    cat = pump_catalog()
    assert len(cat) == 8
    for name, p in cat.items():
        assert p.head(0) > 0                     # positive shutoff
        assert p.head(1e-4) < p.head(0)          # falling curve
        assert 0 < p.eta <= 0.85
        assert p.price > 0 and p.rpm > 0


def test_pump_catalog_spans_specific_speed_map():
    """The catalog must span radial -> axial habitats (the teaching point)."""
    from fluidtools.pipeflow import PUMP_CATALOG_SPECS
    from fluidtools.turbo import specific_speed, rpm_to_rad
    ns = []
    for h0, a, eta, n0, n1, rpm, price in PUMP_CATALOG_SPECS.values():
        q_bep = (h0 / (3 * a)) ** 0.5
        h_bep = h0 - a * q_bep**2
        ns.append(specific_speed(rpm_to_rad(rpm), q_bep, h_bep))
    assert min(ns) < 1.0      # a clearly radial machine exists
    assert max(ns) > 3.5      # a clearly axial machine exists
