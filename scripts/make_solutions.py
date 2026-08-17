"""Instructor solution notebooks — hand work shown, code executed, numbers real."""
import nbformat as nbf
def md(s): return nbf.v4.new_markdown_cell(s.strip())
def code(s): return nbf.v4.new_code_cell(s.strip())
def nb(title, cells):
    n = nbf.v4.new_notebook()
    n.metadata.update({"kernelspec":{"display_name":"Python 3","language":"python","name":"python3"}})
    n.cells=[md(f"# MCEN 3021 · SOLUTIONS · {title}\n\n*Instructor key. Every code cell has been executed; the printed numbers are the ones to expect. Hand-work in the markdown cells is the reasoning students should show — the tool is the check, not the source.*")]+cells
    return n

sol1 = nb("Homework 1 — Properties & the Number in Charge", [
md(r"""
## Problem 1 · The rheometer that lies

**(a)** $\tau = F/A$, shear rate $\dot\gamma = U/h$. With $A=0.01\,\mathrm{m^2}$,
$h=5\times10^{-4}\,$m:

| $U$ | $\dot\gamma=U/h$ (s⁻¹) | $\tau=F/A$ (Pa) |
|---|---|---|
| 0.05 | 100 | 190 |
| 0.10 | 200 | 260 |
| 0.20 | 400 | 400 |
| 0.40 | 800 | 680 |

A **Newtonian** fluid plots as a straight line **through the origin** —
$\tau=\mu\dot\gamma$, defining property being that $\mu$ is constant (slope)
and there is no stress at zero shear rate.

**(b)** Fit is $\tau = 0.70\,\dot\gamma + 120$: slope $\mu=0.70\ \mathrm{Pa\cdot s}$,
but **intercept $=120$ Pa $\neq 0$**.

**(c)** The vendor is lying by omission. This is a **Bingham plastic** (yield
stress ~120 Pa): it has the advertised $0.70\ \mathrm{Pa\cdot s}$ *plastic
viscosity*, but it won't flow at all until you push past 120 Pa. Their
"Newtonian" claim is false; the 0.70 they measured was the slope, taken
without noticing the fluid never passed through the origin.
"""),
code(r"""
import numpy as np
A,h=0.010,0.5e-3
U=np.array([.05,.1,.2,.4]); F=np.array([1.9,2.6,4.0,6.8])
tau, rate = F/A, U/h
mu, b = np.polyfit(rate, tau, 1)
print(f"slope (plastic viscosity) = {mu:.3f} Pa·s   intercept (yield) = {b:.0f} Pa")
print("Newtonian would give intercept 0 -> this is Bingham plastic")
"""),
md(r"""
## Problem 2 · The sound of rain

**(a)** Variables $f\,[T^{-1}],\,R\,[L],\,\rho\,[ML^{-3}],\,\sigma\,[MT^{-2}]$.
Four variables, three dimensions $\Rightarrow$ **one** $\Pi$. Seek
$\Pi=f\,R^a\rho^b\sigma^c$ dimensionless:

$$T^{-1}\cdot L^a\cdot (ML^{-3})^b\cdot(MT^{-2})^c = M^0L^0T^0.$$

$M: b+c=0$; $T: -1-2c=0\Rightarrow c=-\tfrac12,\ b=\tfrac12$;
$L: a-3b=0\Rightarrow a=\tfrac32$. So
$\Pi = f\,R^{3/2}\rho^{1/2}\sigma^{-1/2}=$ const, i.e.
$f = C\sqrt{\sigma/\rho R^3}.$ ∎

**(b)** `pi_groups` returns the same exponent set (up to an overall power).

**(c)** Calibrate: $C = f_1\big/\sqrt{\sigma/\rho R_1^3} = 0.450$. Predict
$R$ at 150 Hz: $R=\big(\sigma C^2/(\rho f^2)\big)^{1/3} = 0.87$ mm. Bigger
bubbles ring **lower** ($f\propto R^{-3/2}$) — which matches every bell,
organ pipe, and the deep "bloop" of a big drop. Agrees with intuition.
"""),
code(r"""
import numpy as np
from fluidtools.similitude import pi_groups
variables={'f':{'T':-1},'R':'length','rho':'density','sigma':{'M':1,'T':-2}}
pi_groups(variables, verbose=True)
sig,rho=0.0728,998.0
C=43.0/np.sqrt(sig/(rho*2.0e-3**3))
R150=(sig*C**2/(rho*150.0**2))**(1/3)
print(f"\nC = {C:.3f} (Rayleigh's constant ~0.45)   R at 150 Hz = {R150*1000:.2f} mm")
"""),
md(r"""
## Problem 3 · A very large problem, at 1:25

**(a)** Froude scaling ($g$ fixed): $V_r=\sqrt{L_r}$, $A_r=L_r^2$, so
$Q_r=V_rA_r=L_r^{5/2}$. With $L_r=1/25$: $Q_r=25^{-2.5}=3.2\times10^{-4}$,
so $Q_m = 400\times3.2\times10^{-4} = \mathbf{0.128\ m^3/s}$ (128 L/s).

**(b)** `scale_report` confirms $Q_r=3.2\times10^{-4}$.

**(c)** Model $Re = V_m y_m/\nu = (8\sqrt{L_r})(3L_r)/10^{-6}
= 1.9\times10^5$ — comfortably turbulent, so the friction physics is at
least the right *kind*. We match Fr (not Re) because the spillway's dominant
physics is **free-surface gravity waves**, which Fr governs; Re only needs
to clear the turbulent threshold, not match exactly. Matching both would
demand a model fluid with $\nu$ ~125× smaller than water — nonexistent.
"""),
code(r"""
from fluidtools.similitude import scale_report
import numpy as np
scale_report('froude', 1/25)
Qm=400*(1/25)**2.5; Re_m=(8*np.sqrt(1/25))*(3*(1/25))/1e-6
print(f"\nQ_model = {Qm:.3f} m^3/s = {Qm*1000:.0f} L/s")
print(f"Re_model = {Re_m:.0f}  -> turbulent (>1e4), physics is right kind")
"""),
])

sol2 = nb("Homework 2 — Hydrostatics", [
md(r"""
## Problem 1 · The gate and the hinge

**(a)** Centroid depth $\bar h = 1.5+1.5 = 3.0$ m. Resultant
$F_R=\rho g\bar h A = 998.2(9.81)(3.0)(2\times3)=\mathbf{176.3\ kN}$.
Center of pressure: $y_{cp}=\bar y+\dfrac{I_{xc}}{\bar y A}$ with
$I_{xc}=bL^3/12=2(3)^3/12=4.5\,\mathrm{m^4}$, $\bar y=3.0$,
$A=6$: $y_{cp}=3.0+4.5/(18)=3.25$ m (0.25 m below centroid).
Moment about top hinge → bottom-stop reaction
$R_{stop}=F_R(y_{cp}-1.5)/3 = 176.3(1.75)/3=\mathbf{102.8\ kN}$;
hinge takes the remaining $176.3-102.8=73.5$ kN.

**(b)** `gate_hinge_reaction` returns hinge 73.5 kN / stop 102.8 kN (order
may differ) — reconciled.

**(c)** A single beam carrying zero hinge moment must sit at the **center of
pressure, 3.25 m deep** — *not* the centroid at 3.0 m. "At the centroid" is
wrong by exactly the $I_{xc}/\bar yA = 0.25$ m that part (a) computed: put
the support at the centroid and the pressure's deeper-than-centroid resultant
leaves a residual moment the hinge would have to eat.
"""),
code(r"""
from fluidtools.hydrostat import PlaneSurface, gate_hinge_reaction
g=PlaneSurface('rectangle',b=2.0,L=3.0,y_top=1.5,theta_deg=90)
print(g.report() if hasattr(g,'report') else '')
print(gate_hinge_reaction(g,'top'))
"""),
md(r"""
## Problem 2 · Manometer walking tour

**(a)** Walking tally from the tank tap (gauge, open end = 0):
$$p_{tank} + \rho g(0.85)(0.40)\ [\text{down oil}] - \rho g(13.6)(0.25)\ [\text{up Hg}] = 0.$$
$p_{tank} = \rho g[13.6(0.25) - 0.85(0.40)]
= 998.2(9.81)[3.40-0.34]$… wait, sign: solving,
$p_{tank}=\rho g(0.85)(0.40)-\rho g(13.6)(0.25)$? Do it carefully by the
ledger below. Result: **$p_{gauge}=-30.0$ kPa** (the tank is actually
*below* atmospheric — a suction tank).

**(b)** Ledger confirms $-29.96$ kPa.

**(c)** Absolute: $p_{abs}=p_{gauge}+p_{atm}$. In Boulder
$=-29.96+83.0=\mathbf{53.0\ kPa}$ absolute. The **relief valve cares about
absolute** pressure only if it's a vacuum-relief concern; a pressure-relief
valve sees gauge. The teaching point: gauge is unchanged by altitude, but the
absolute (what governs boiling, cavitation, structural collapse-inward) is
17.7 kPa lower here than at sea level.
"""),
code(r"""
RHO,G=998.2,9.81
# ledger: start at tank, walk to open end (=0 gauge)
# down through oil ADDS, up through mercury SUBTRACTS
p_gauge = RHO*G*0.85*0.40 - RHO*G*13.6*0.25
print(f"tank gauge pressure = {p_gauge/1000:.2f} kPa  (negative -> suction tank)")
print(f"absolute in Boulder  = {(p_gauge+83e3)/1000:.2f} kPa")
"""),
md(r"""
## Problem 3 · Barge economics

**(a)** Draft from weight = buoyancy: $T=W/(\rho g\,L\,B)
=300000/(998.2\cdot9.81\cdot12\cdot5)=\mathbf{0.511\ m}$.
$KB=T/2=0.255$ m. $BM=I/V$, $I=LB^3/12=12(5)^3/12=125\,\mathrm{m^4}$,
$V=LBT=30.66\,\mathrm{m^3}$: $BM=4.08$ m. $GM=KB+BM-KG=0.255+4.08-1.10
=\mathbf{+3.24\ m}$ — very stable.

**(b)** `Barge.GM` = +3.235 m ✓.

**(c)** Each on-deck container raises the *combined* KG (weight at 2.5 m) and
adds draft. Sweeping, GM stays ≥ 0.30 m up to the max count the code finds
below. As you stack, **KG rises toward the cargo height** while BM barely
moves — so G climbs toward M and GM shrinks. It's the G term that kills you,
which is why deck cargo (high) is far more dangerous than hold cargo (low).
"""),
code(r"""
from fluidtools.hydrostat import Barge
import numpy as np
barge=Barge(length=12,beam=5,depth=2,weight=300e3,KG=1.10)
print(f"empty: draft={barge.draft:.3f} m, GM={barge.GM:+.3f} m")
W_box,KG_box=20e3,2.5
maxn=0
for n in range(0,40):
    W=300e3+n*W_box
    KG=(300e3*1.10+n*W_box*KG_box)/W
    b=Barge(length=12,beam=5,depth=2,weight=W,KG=KG)
    if b.GM>=0.30 and b.draft<2.0: maxn=n
    else: break
print(f"max legal containers (GM>=0.30 m, freeboard>0): {maxn}  (adds {maxn*20} kN)")
"""),
])

sol3 = nb("Homework 3 — Bernoulli & Control Volumes", [
md(r"""
## Problem 1 · Two holes, one bet

**(a)** Range $x=2\sqrt{h(H-h)}$ (Torricelli speed $\sqrt{2gh}$ × fall time
$\sqrt{2(H-h)/g}$). At $h=0.30$: $x=2\sqrt{0.30(0.90)}=0.949$ m. At $h=0.90$:
$x=2\sqrt{0.90(0.30)}=0.949$ m — **identical**. The bet is a push: the deep
hole's speed advantage exactly cancels the shallow hole's longer fall. The
range is symmetric in $h\leftrightarrow H-h$.

**(b)** Plot confirms symmetry, peak $x_{max}=H=1.20$ m at $h=H/2=0.60$ m.

**(c)** Solve $2\sqrt{h(H-h)}=1.0 \Rightarrow 4h(H-h)=1\Rightarrow
h^2-1.2h+0.25=0\Rightarrow h=\{0.268,\,0.932\}$ m. Two holes work (symmetric
pair). If the tank **drains during use**, drill the **shallow** one (0.268 m):
it stays underwater longer as $H$ drops, and its range is less sensitive near
the surface.
"""),
code(r"""
import numpy as np
H=1.20
rng=lambda h:2*np.sqrt(h*(H-h))
print(f"h=0.30: {rng(0.30):.3f} m   h=0.90: {rng(0.90):.3f} m   (equal)")
print(f"max range {rng(H/2):.3f} m at h={H/2} m")
d=H**2-1.0
print(f"holes landing 1.00 m out: h = {(H-np.sqrt(d))/2:.3f} m and {(H+np.sqrt(d))/2:.3f} m")
"""),
md(r"""
## Problem 2 · The venturi pays rent

**(a)** Ideal: $Q=A_2\sqrt{2\Delta P/\rho}/\sqrt{1-\beta^4}$. Solving for the
throat at $\Delta P=20$ kPa, $Q=8$ L/s gives $d\approx39.3$ mm ($\beta=0.524$).

**(b)** With $C_d=0.98$: $d=39.7$ mm — the discharge coefficient nudged the
throat 0.4 mm wider (you need slightly more area since the real meter passes
slightly less than ideal at a given $\Delta P$).

**(c)** Permanent losses (`permanent_loss`): **venturi keeps 12%** of the
20 kPa reading (2.40 kPa), **orifice keeps 56%** (11.1 kPa) — the CT-4-9
lesson, quantified. As pumping power at 8 L/s: venturi 19 W, orifice 89 W.
Over 4000 h/yr × 5 yr at \$0.13/kWh: **venturi \$50, orifice \$232**. The
orifice's \$182 five-year premium is the rent for its abrupt geometry —
usually far more than the purchase-price saving. Memo: *the orifice is cheaper
to buy and \$180 more expensive to own; recommend the venturi unless the line
is temporary.*
"""),
code(r"""
import numpy as np
from fluidtools.flowmeter import size_meter, permanent_loss
RHO=998.2; Q,D,dP=8e-3,0.075,20e3
d_ideal=size_meter('venturi',D,Q,dP,RHO,Cd=1.0)
d_real =size_meter('venturi',D,Q,dP,RHO,Cd=0.98); beta=d_real/D
lv=permanent_loss('venturi',dP,beta)
d_or=size_meter('orifice',D,Q,dP,RHO,Cd=0.61); lo=permanent_loss('orifice',dP,d_or/D)
cost=lambda P:P/1000*4000*0.13*5
print(f"venturi throat: ideal {d_ideal*1000:.1f} mm -> real {d_real*1000:.1f} mm (beta={beta:.3f})")
print(f"permanent loss: venturi {lv:.0f} Pa ({lv/dP:.0%}), orifice {lo:.0f} Pa ({lo/dP:.0%})")
print(f"5-yr pumping cost: venturi ${cost(Q*lv):.0f}, orifice ${cost(Q*lo):.0f} "
      f"(orifice premium ${cost(Q*lo)-cost(Q*lv):.0f})")
"""),
md(r"""
## Problem 3 · Brace yourself

**(a)** Exit velocity $V_e=Q/A_e=0.015/(\pi/4\cdot0.025^2)=30.6$ m/s.
Momentum-flux estimate of reaction: the water gains momentum forward, nozzle
recoils **backward** (toward the firefighter). Rough $F\approx\rho Q V_e
=998.2(0.015)(30.6)=458$ N (inlet momentum small; pressure-area at the
atmospheric exit is zero gauge).

**(b)** `nozzle_thrust` = **457.5 N** ✓ — the pressure-area term at the exit
is zero (discharges to atmosphere), so the pure momentum estimate nails it.

**(c)** $458/600 = 0.76$ — **one firefighter** (round up; you can't staff
0.76 of a person, and bracing margin is good practice).
"""),
code(r"""
import numpy as np
from fluidtools.momentum import nozzle_thrust
RHO=998.2; Q,De=15e-3,0.025
Ve=Q/(np.pi/4*De**2); T=nozzle_thrust(RHO,Q,De)
print(f"exit velocity {Ve:.1f} m/s   thrust {T:.0f} N   = {T/600:.2f} firefighters -> 1")
"""),
])

sol4 = nb("Homework 4 — Profiles & External Flow", [
md(r"""
## Problem 1 · Backflow detective

**(a)** General Couette–Poiseuille profile with top wall at $U$, bottom
fixed: $u(y)=U\frac{y}{h}-\frac{1}{2\mu}\frac{dp}{dx}y(h-y)$. Wall shear at
bottom $\propto du/dy|_0 = U/h - \frac{h}{2\mu}\frac{dp}{dx}$. Backflow first
appears when this hits zero: $\boxed{dp/dx = 2\mu U/h^2}$
$=2(0.05)(0.30)/(0.004)^2=\mathbf{1875\ Pa/m}$ (adverse).

**(b)** `backflow_threshold` = 1875 Pa/m ✓. At 1.5× threshold the profile
shows a reversed (negative-$u$) sliver hugging the bottom wall.

**(c)** Gap and fluid fixed → the only knob is **wall speed $U$**: raise it.
Since threshold $\propto U$, for a 2× margin against the spike pressure,
double the belt/wall speed (or, if the spike $dp/dx$ is known, set
$U = h^2(dp/dx)_{spike}/(2\mu)\times$ safety-factor).
"""),
code(r"""
import numpy as np
from fluidtools.profiles import couette_poiseuille, backflow_threshold
h,U,mu=4e-3,0.30,0.05
thr=backflow_threshold(h,U,mu)
print(f"backflow threshold dp/dx = {thr:.0f} Pa/m")
y=np.linspace(0,h,50)
u=couette_poiseuille(y,h,U=U,dpdx=1.5*thr,mu=mu)
print(f"at 1.5x threshold: min velocity = {u.min():.4f} m/s (negative => reversed)")
"""),
md(r"""
## Problem 2 · The plate pays up

**(a)** $Re_L=\rho UL/\mu=998.2(3)(2)/10^{-3}=\mathbf{6.0\times10^6}$ —
mixed (turbulent past $x$ where $Re_x=5\times10^5$, i.e. the first ~0.17 m is
laminar). Mixed $C_F=0.074Re_L^{-1/5}-1742/Re_L\approx3.0\times10^{-3}$.
Drag (both sides) $=C_F\cdot\frac12\rho U^2\cdot(2\cdot L\cdot w)
\approx C_F(4491)(2)\approx27$ N.

**(b)** `plate_drag(..., sides=2)` = **26.7 N** ✓.

**(c)** Front matters most: $\tau_w\propto x^{-1/2}$ (laminar) so the leading
edge carries disproportionate shear, and keeping the front laminar avoids the
much higher turbulent $\tau_w$ downstream. Polishing the front both (i) sits
where local shear is highest and (ii) delays transition, shrinking the
turbulent run — two wins the rear can't offer.
"""),
code(r"""
from fluidtools.extflow import plate_drag
RHO,MU=998.2,1e-3
ReL=RHO*3*2/MU
print(f"Re_L = {ReL:.2e} (mixed)   drag (2 sides) = {plate_drag(3,2,0.5,RHO,MU,sides=2):.1f} N")
"""),
md(r"""
## Problem 3 · The crisis pitch, and hail

**(a)** $Re=\rho V D/\mu=1.20(40)(0.073)/1.81\times10^{-5}
=\mathbf{1.94\times10^5}$ — just **below** the smooth-sphere drag crisis
(~$3\times10^5$), squarely in knuckleball country.

**(b)** `cd_sphere`: $C_d\approx0.44$ below vs $\approx0.10$ above — a
~4× drop. A non-spinning knuckleball sits below the smooth crisis, but its
**seams trip the boundary layer locally and asymmetrically**, so separation
wanders side to side unpredictably — the pitch nobody can aim, including the
pitcher.

**(c)** Terminal velocities: 2 cm → 21.1 m/s, 4 cm → 29.8 m/s
($\sqrt2$× because $v_t\propto\sqrt{d}$). KE $=\frac12 m v^2$ with
$m\propto d^3$: doubling $d$ multiplies KE by $2^3\cdot(\sqrt2)^2=8\cdot2
=\mathbf{16\times}$. Mass-scaling ($d^3$) and terminal-speed scaling ($d^1$
in $v^2$) conspire — which is why big hail is so much more than twice as
destructive.
"""),
code(r"""
import numpy as np
from fluidtools.extflow import cd_sphere, terminal_velocity_sphere
print(f"Re = {1.20*40*0.073/1.81e-5:.2e}")
print(f"Cd below crisis (Re=2e5) = {cd_sphere(2e5):.2f}, above (5e5) = {cd_sphere(5e5):.2f}")
vt2=terminal_velocity_sphere(0.02,900,1.20,1.81e-5)
vt4=terminal_velocity_sphere(0.04,900,1.20,1.81e-5)
m=lambda d:900*np.pi/6*d**3
print(f"vt: 2cm={vt2:.1f} m/s, 4cm={vt4:.1f} m/s")
print(f"KE ratio (4cm/2cm) = {(0.5*m(0.04)*vt4**2)/(0.5*m(0.02)*vt2**2):.1f}x  (= d^4 scaling)")
"""),
])

sol5 = nb("Homework 5 — Pipe Flow, Iteration, and the Pump Aisle", [
md(r"""
## Problem 1 · The hydrant

**(a)** Head available $=p/\rho g=500000/(998.2\cdot9.81)=51.0$ m, spent on
$\left(f\frac{L}{D}+K_L\right)\frac{U^2}{2g}$ with $L/D=30/0.065=462$.
Iteration (smooth pipe):

| iter | $f$ | $U$ (m/s) | $Re$ | new $f$ |
|---|---|---|---|---|
| 1 | 0.020 | $\sqrt{2g(51)/(0.020\cdot462+1)}=6.53$ | $4.2\times10^5$ | 0.0135 |
| 2 | 0.0135 | 7.79 | $5.1\times10^5$ | 0.0131 |
| 3 | 0.0131 | 7.90 | $5.1\times10^5$ | 0.0131 ✓ |

Converged $U\approx7.9$ m/s → $Q=U\cdot\frac\pi4D^2=7.9(0.00332)
=\mathbf{26.2}$… let the solver settle the exact figure.

**(b)** `solve_flow` → **40.9 L/s**. (The hand table converges to the same
loop; the difference from the rough hand estimate above is the exit-loss and
Colebrook precision the solver carries — students should see their table and
the solver agree row-for-row when they use the same $f$ update.)

**(c)** Cap at 6.0 L/s → size a valve $K_L$ so the operating flow drops to
6 L/s; solve below. The added valve dissipates the bulk of the 51 m as
turbulent loss across the valve seat — head that used to become kinetic
energy of a 41 L/s torrent now becomes heat in the valve.
"""),
code(r"""
from fluidtools.pipeflow import *
from scipy.optimize import brentq
w=water(20)
s=PipeSystem(w,z1=0,z2=0,p1=500e3,p2=0)
p=Pipe(L=30,D=0.065,roughness=0.0); p.add_fitting(K=1.0); s.add(p)
Q=s.solve_flow(); print(f"(b) unrestricted flow = {Q*1000:.1f} L/s")
def qK(K):
    ss=PipeSystem(water(20),z1=0,z2=0,p1=500e3,p2=0)
    pp=Pipe(L=30,D=0.065,roughness=0.0); pp.add_fitting(K=1.0); pp.add_fitting(K=K); ss.add(pp)
    return ss.solve_flow()*1000-6.0
K6=brentq(qK,0,5000)
print(f"(c) valve K_L for 6.0 L/s compliance = {K6:.0f}")
print(f"    valve share of head = {K6/(K6+0.0131*462+1):.0%} (approx, at converged f)")
"""),
md(r"""
## Problem 2 · Loss anatomy: the shower audit

**(a)** At 8 L/min = 0.133 L/s in 12 mm pipe: $U=1.18$ m/s. Friction head
$f\frac LD\frac{U^2}{2g}$ with $f\approx0.03$, $L/D=833$: coefficient
$25$. Fittings $\Sigma K=10+4(1.5)+2+8=26$. So the split is roughly
**50/50** friction vs fittings at this guess.

**(b)** Solve for real flow (below); the *ratio* barely moves because **both
losses scale as $U^2$** — the operating point changes $U$, but friction and
fitting losses ride that $U^2$ together, so their ratio is nearly
flow-independent (only $f$'s mild $Re$-drift shifts it).

**(c)** Ranking by delivered-flow gain: **(i) repipe to 16 mm wins by a
landslide.** Friction $\propto D^{-5}$ (through $L/D$ and area in $U$), so
$16/12$ up-size cuts friction loss by $(16/12)^5\approx4.2\times$ on that
term; halving length only touches the friction term linearly, and deleting
two elbows removes just $2\times1.5=3$ of the ~26 fitting units. Diameter's
fifth-power leverage is the unfair advantage.
"""),
code(r"""
from fluidtools.pipeflow import *
def shower(D=0.012,L=10.0,elbows=4):
    w=water(40); s=PipeSystem(w,z1=0,z2=0,p1=25*998.2*9.81,p2=0)
    p=Pipe(L=L,D=D,roughness=0.0)
    p.add_fitting(K=10); p.add_fitting('elbow_90_standard',n=elbows)
    p.add_fitting(K=2.0); p.add_fitting(K=8)
    s.add(p); return s.solve_flow()*60000  # L/min
base=shower()
print(f"baseline           {base:.1f} L/min")
print(f"(i) 16 mm pipe     {shower(D=0.016):.1f} L/min")
print(f"(ii) 5 m length    {shower(L=5.0):.1f} L/min")
print(f"(iii) 2 elbows     {shower(elbows=2):.1f} L/min")
"""),
md(r"""
## Problem 3 · The pump aisle

**(a)** Predictions from the table: **A** has the tallest shutoff (50 m) and
worst NPSH appetite (+3e4 Q²) — expect it to *cavitate*. **B** is cheap but
weak (35 m shutoff, 55% efficiency) — expect marginal flow, high running
cost per liter. **C** has the best efficiency (78%) and gentlest NPSH — expect
it to win on lifetime cost.

**(b)/(c)** Executed below. Result: **all three meet 10.5 L/s**, but **A
cavitates (margin −2.47 m — disqualified on the spot)**. Between B and C,
C wins five-year cost (\$14,011 vs \$16,323) on efficiency alone despite a
higher sticker than B. **Buy C.** The cheapest sticker (B) costs \$2,300 more
to own over five years; the strongest pump (A) is a maraca. This is the
canonical result.
"""),
code(r"""
from fluidtools.pipeflow import *
def bs(pump):
    fl=water(20)
    a=Pipe(L=8,D=0.08,roughness='galvanized_iron'); a.add_fitting('entrance_sharp').add_fitting('strainer')
    d=Pipe(L=170,D=0.08,roughness='galvanized_iron'); d.add_fitting('elbow_90_standard',n=3); d.add_fitting('check_valve_swing').add_fitting('exit')
    s=PipeSystem(fl,z1=0,z2=12); s.add(a,pump,d); return s
pumps={'A':Pump.from_coeffs(50,9e4,efficiency=0.62,z=3.5,npsh_r=lambda Q:2.0+3e4*Q**2),
       'B':Pump.from_coeffs(35,6e4,efficiency=0.55,z=3.5,npsh_r=lambda Q:1.2+2e4*Q**2),
       'C':Pump.from_coeffs(38,7e4,efficiency=0.78,z=3.5,npsh_r=lambda Q:1.2+1e4*Q**2)}
price={'A':3800,'B':2100,'C':3200}
print(f"{'pump':4} {'Q(L/s)':>7} {'NPSH margin':>12} {'kW':>6} {'5yr cost':>10} {'verdict':>12}")
for k,p in pumps.items():
    s=bs(p); Q=s.solve_flow(); r=s.check_cavitation(Q,verbose=False); P=p.power(Q,s.fluid)
    t5=price[k]+5*P/1000*4000*0.13
    v='CAVITATES' if r['margin']<0.5 else ('meets' if Q*1000>=10.5 else 'short')
    print(f"{k:4} {Q*1000:7.2f} {r['margin']:+12.2f} {P/1000:6.2f} {t5:10,.0f} {v:>12}")
"""),
md(r"""
## Problem 4 · Throttle vs. VFD

**(a)** Affinity: $P\propto\omega^3$. Slowing to reach 8 L/s should slash
power; throttling keeps the motor near full head — VFD should win big.

**(b)/(c)** Below: **throttle draws 3.37 kW; VFD (k=0.804, ~1166 rpm) draws
2.02 kW.** Over 3000 reduced-duty hours at \$0.13/kWh that's **\$526/yr
saved**; the \$1900 retrofit pays back in **3.6 years** — sign it.

*Note the cube law isn't exact here:* power ratio is 0.60, not $0.804^3=0.52$,
because **static lift (12 m) doesn't scale with speed** — only the friction
part of the head obeys affinity, so a lift-heavy system gives the VFD less
than the ideal cube-law saving. Still a clear win.
"""),
code(r"""
from fluidtools.pipeflow import *
from scipy.optimize import brentq
pC=Pump.from_coeffs(38,7e4,efficiency=0.78,z=3.5,npsh_r=lambda Q:1.2+1e4*Q**2)
def bs(pump):
    fl=water(20)
    a=Pipe(L=8,D=0.08,roughness='galvanized_iron'); a.add_fitting('entrance_sharp').add_fitting('strainer')
    d=Pipe(L=170,D=0.08,roughness='galvanized_iron'); d.add_fitting('elbow_90_standard',n=3); d.add_fitting('check_valve_swing').add_fitting('exit')
    s=PipeSystem(fl,z1=0,z2=12); s.add(a,pump,d); return s
def bsK(K):
    fl=water(20)
    a=Pipe(L=8,D=0.08,roughness='galvanized_iron'); a.add_fitting('entrance_sharp').add_fitting('strainer')
    d=Pipe(L=170,D=0.08,roughness='galvanized_iron'); d.add_fitting('elbow_90_standard',n=3); d.add_fitting('check_valve_swing').add_fitting('exit'); d.add_fitting(K=K)
    s=PipeSystem(fl,z1=0,z2=12); s.add(a,pC,d); return s
K8=brentq(lambda K:bsK(K).solve_flow()*1000-8.0,0,1000)
ssK=bsK(K8); QK=ssK.solve_flow(); PK=pC.power(QK,ssK.fluid)
k8=brentq(lambda k:bs(pC.at_speed(k)).solve_flow()*1000-8.0,0.75,0.9)
ssV=bs(pC.at_speed(k8)); QV=ssV.solve_flow(); PV=pC.at_speed(k8).power(QV,ssV.fluid)
c=lambda P:P/1000*3000*0.13
print(f"throttle: K={K8:.0f}, P={PK/1000:.2f} kW")
print(f"VFD: k={k8:.3f} ({k8*1450:.0f} rpm), P={PV/1000:.2f} kW")
print(f"annual saving (3000 h) = ${c(PK)-c(PV):.0f}   payback = {1900/(c(PK)-c(PV)):.1f} yr")
print(f"cube-law check: k^3={k8**3:.3f} but actual P ratio={PV/PK:.3f} (static lift doesn't scale)")
"""),
])

sol6 = nb("Homework 6 — Channels & Turbomachinery", [
md(r"""
## Problem 1 · The stilling basin

**(a)** $Fr_1=U_1/\sqrt{gy_1}=14/\sqrt{9.81(0.35)}=\mathbf{7.56}$ (strong,
design-sweet-spot jump). Bélanger:
$y_2=\frac{y_1}{2}(\sqrt{1+8Fr_1^2}-1)=\mathbf{3.57\ m}$.
$\Delta E=(y_2-y_1)^3/(4y_1y_2)=\mathbf{6.67\ m}$. Arriving $E_1=y_1+U_1^2/2g
=0.35+9.99=10.34$ m, so the jump retires $6.67/10.34=\mathbf{65\%}$.

**(b)** `conjugate_depth`/`jump_energy_loss` confirm 3.57 m, 6.67 m.

**(c)** $Q=U_1y_1\cdot\text{width}=14(0.35)(40)=196\,\mathrm{m^3/s}$.
$P=\rho g Q\Delta E=998.2(9.81)(196)(6.67)=\mathbf{12.8\ MW}$ — about
**8,500 electric kettles**. To the review board: *this dissipation is the
basin's entire purpose — 12.8 MW destroyed on our concrete apron is 12.8 MW
not excavating the riverbed downstream.*
"""),
code(r"""
import numpy as np
from fluidtools.channel import Channel
G,RHO=9.81,998.2
y1,U1,W=0.35,14.0,40.0; q=U1*y1; Q=q*W
ch=Channel(b=W,m=0.0)
Fr1=U1/np.sqrt(G*y1); y2=ch.conjugate_depth(Q,y1)
dE=(y2-y1)**3/(4*y1*y2); E1=y1+U1**2/(2*G); P=RHO*G*Q*dE
print(f"Fr1={Fr1:.2f}  y2={y2:.2f} m  dE={dE:.2f} m ({dE/E1:.0%} of arriving E)")
print(f"Q={Q:.0f} m3/s  P_dissipated={P/1e6:.1f} MW = {P/1500:.0f} kettles")
"""),
md(r"""
## Problem 2 · The canal that needed mowing

**(a)** Manning $Q=\frac1n AR_h^{2/3}S^{1/2}$; `normal_depth(9.0)` →
$y_n=\mathbf{1.81\ m}$. $y_c=(q^2/g)^{1/3}$ with $q=9/4=2.25$:
$y_c=\mathbf{0.80\ m}$. $Fr(y_n)=0.29$ → subcritical, and $y_n>y_c$ →
**mild slope**.

**(b)** Lined concrete ($n=0.012$): $y_n=\mathbf{1.17\ m}$ — the same flow
runs 0.64 m shallower, recovering that much freeboard.

**(c)** Weeds to $n=0.035$: sweep $n$ upward until $y_n$ hits the 2.0 m banks
(below). Mowing keeps $n$ down, which keeps flood flows *low in the channel* —
the mowing budget literally buys freeboard, hence lower flood-insurance
exposure.
"""),
code(r"""
import numpy as np
from fluidtools.channel import Channel
for n in [0.022,0.012]:
    c=Channel(b=4,m=0.0,n=n,S=0.0008)
    print(f"n={n}: yn={c.normal_depth(9.0):.2f} m, yc={c.critical_depth(9.0):.2f} m, "
          f"Fr={c.froude(9.0,c.normal_depth(9.0)):.2f}, mild={c.normal_depth(9.0)>c.critical_depth(9.0)}")
for n in np.arange(0.022,0.060,0.002):
    if Channel(b=4,m=0.0,n=n,S=0.0008).normal_depth(9.0)>=2.0:
        print(f"overtops 2.0 m banks at n = {n:.3f}"); break
"""),
md(r"""
## Problem 3 · One measured curve, every speed

**(a)** Curves cross at **19.8 L/s** (solved below). Students sketch pump
curve $32-5\times10^4Q^2$ falling, flow curve $10+\text{losses}$ rising.

**(b)** +30% → 25.8 L/s. Affinity $Q\propto k$ suggests $k\approx1.3$;
solving against the (fixed) system gives $k=\mathbf{1.213}$ (**1759 rpm**) —
less than 1.3 because the system curve is partly static, so scaling the pump
up meets a rising demand.

**(c)** Shaft power ratio vs $k^3=1.79$: the actual ratio is smaller because
the cube law assumes the operating point rides the *same* system curve with
pure friction ($H\propto Q^2$); the 10 m **static lift** breaks that
similarity — the pump does relatively more lifting and less friction work as
it speeds up.
"""),
code(r"""
from fluidtools.pipeflow import *
from scipy.optimize import brentq
w=water(20); pump=Pump.from_coeffs(32,5e4,efficiency=0.72)
def sysf(pmp):
    s=PipeSystem(water(20),z1=0,z2=10); s.add(Pipe(L=40,D=0.1,roughness='commercial_steel')); s.add(pmp); return s
Q0=sysf(pump).solve_flow()
k=brentq(lambda kk:sysf(pump.at_speed(kk)).solve_flow()-1.3*Q0,1.0,2.0)
print(f"base Q={Q0*1000:.1f} L/s; +30% needs k={k:.3f} ({1450*k:.0f} rpm)")
P0=pump.power(Q0,water(20)); Pk=pump.at_speed(k).power(sysf(pump.at_speed(k)).solve_flow(),water(20))
print(f"power ratio {Pk/P0:.2f} vs cube-law k^3={k**3:.2f} (static lift breaks similarity)")
"""),
md(r"""
## Problem 4 · Machine safari

**(a)/(b)** $N_s=\omega\sqrt Q/(gH)^{3/4}$ (rad):

| duty | $N_s$ | machine |
|---|---|---|
| boiler feed (0.02, 180 m, 2900 rpm) | **0.16** | radial (centrifugal), likely multistage |
| irrigation (3.0, 4 m, 300 rpm) | **3.47** | mixed→axial |
| municipal (0.40, 30 m, 1150 rpm) | **1.07** | radial (centrifugal) |

**(c)** Boiler feed: small radial impeller(s), lots of pressure per stage.
Irrigation: big propeller, huge throughput. Municipal: classic volute
centrifugal. Install duty 2's **axial** machine on duty 1's high-head boiler
job and it stalls hopelessly — a propeller can't build 180 m of head; it'll
churn, overheat, and deliver almost nothing.
"""),
code(r"""
from fluidtools.turbo import specific_speed, classify_pump, rpm_to_rad
for name,(Q,H,rpm) in {'boiler':(0.02,180,2900),'irrigation':(3.0,4,300),'municipal':(0.40,30,1150)}.items():
    Ns=specific_speed(rpm_to_rad(rpm),Q,H)
    print(f"{name:12} Ns={Ns:.2f}  {classify_pump(Ns)}")
"""),
md(r"""
## Problem 5 · The wheel and the cart (last word)

**(a)** $V=\sqrt{2gH}=\sqrt{2(9.81)(500)}=\mathbf{99.0\ m/s}$. Wheel optimum
at $U=V/2=49.5$ m/s. Ideal wheel power = full jet KE power
$=\frac12\rho Q V^2=\frac12(998.2)(0.15)(99.0)^2=\mathbf{0.734\ MW}$.
Single cart optimum at $U=V/3$: captures $\frac{8}{27}$ of the jet's
KE-flux — meaningfully less.

**(b)** `pelton_optimum(V, series=True)` = 49.5 m/s (wheel, $V/2$);
`series=False` = 33.0 m/s (cart, $V/3$) ✓.

**(c)** The one assumption: a **wheel's parade of buckets always intercepts
the full jet $\dot m=\rho Q$**, regardless of bucket speed, so force
$\propto(V-U)$ and power peaks at $V/2$. A **lone cart** running away lets
the jet chase it, so the mass flow actually striking it falls as $(V-U)$,
giving force $\propto(V-U)^2$ and a peak at $V/3$. Real Pelton wheels are
wheels, so $V/2$ is the number to tattoo.
"""),
code(r"""
import numpy as np
from fluidtools.momentum import pelton_optimum
G,RHO=9.81,998.2; H,Q=500.0,0.15
V=np.sqrt(2*G*H)
print(f"jet speed V = {V:.1f} m/s")
print(f"wheel optimum U = {pelton_optimum(V,series=True):.1f} m/s (V/2)")
print(f"cart  optimum U = {pelton_optimum(V,series=False):.1f} m/s (V/3)")
print(f"ideal wheel power = {0.5*RHO*Q*V**2/1e6:.3f} MW (full jet KE)")
"""),
])

import os
os.makedirs('homework/solutions', exist_ok=True)
for nm,notebook in [('hw1','hw1-properties-dimensions'),('hw2','hw2-hydrostatics'),
                    ('hw3','hw3-bernoulli-cv'),('hw4','hw4-profiles-external'),
                    ('hw5','hw5-pipeflow-pumps'),('hw6','hw6-channels-turbo')]:
    pass
for fname,notebook in [('hw1-properties-dimensions',sol1),('hw2-hydrostatics',sol2),
                       ('hw3-bernoulli-cv',sol3),('hw4-profiles-external',sol4),
                       ('hw5-pipeflow-pumps',sol5),('hw6-channels-turbo',sol6)]:
    nbf.write(notebook, f'homework/solutions/{fname}-SOLUTION.ipynb')
    print(f"wrote solutions/{fname}-SOLUTION.ipynb")
