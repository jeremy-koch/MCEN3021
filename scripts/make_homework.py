"""Generate the six MCEN 3021 homework notebooks (student versions)."""
import nbformat as nbf

def nb(title, cells):
    n = nbf.v4.new_notebook()
    n.metadata.update({"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}})
    n.cells = [nbf.v4.new_markdown_cell(HEADER.format(title=title))] + cells
    return n

def md(s): return nbf.v4.new_markdown_cell(s.strip())
def code(s): return nbf.v4.new_code_cell(s.strip())

HEADER = """# MCEN 3021 · {title}

**How these work (read once, applies all semester).** Every problem runs the
same arc: **(a) by hand** — paper, pencil, the head equation, your brain;
**(b) verify** — the toolbox reproduces your number, or you find out whose
fault that is; **(c) apply** — with the tool now trusted, answer the question
that was too tedious for paper. Hand work goes in the markdown cells (photos
of paper are fine on Canvas); code goes in the code cells. A number with no
units is not an answer. A number with no sentence is barely one.

*Setup:* `pip install fluidtools` (installation instructions on the course
site → Tools)."""

# ----------------------------------------------------------------- HW 1
hw1 = nb("Homework 1 — Properties & the Number in Charge (Units 1–2)", [
md("""
## Problem 1 · The rheometer that lies

A vendor ships you "premium Newtonian damping oil" and a spec sheet claiming
$\\mu = 0.70\\ \\mathrm{Pa\\cdot s}$. Your lab's flat-plate rheometer (plate
area $A = 0.010\\ \\mathrm{m^2}$, gap $h = 0.50$ mm) measures the force to
drag the plate at four speeds:

| $U$ (m/s) | 0.05 | 0.10 | 0.20 | 0.40 |
|---|---|---|---|---|
| $F$ (N) | 1.90 | 2.60 | 4.00 | 6.80 |

**(a) By hand:** convert each row to shear stress $\\tau$ and shear rate
$U/h$. What *should* a Newtonian fluid's $\\tau$ vs. $U/h$ plot look like —
and cite the defining property, not just "a line."

**(b) Verify:** plot your four points and fit a line (`np.polyfit`). Report
slope and intercept with units.

**(c) Apply:** is the vendor lying? One sentence to them, professional but
unambiguous, that names what kind of fluid this actually is and what the
0.70 they measured probably was.
"""),
code("""
import numpy as np
import matplotlib.pyplot as plt

A, h = 0.010, 0.5e-3            # m^2, m
U = np.array([0.05, 0.10, 0.20, 0.40])   # m/s
F = np.array([1.90, 2.60, 4.00, 6.80])   # N

# your work below: tau, rate, plot, polyfit
...
"""),
md("""
## Problem 2 · The sound of rain

Lecture derived, by Buckingham alone, that a bubble of radius $R$ in water
(density $\\rho$, surface tension $\\sigma$) oscillates at
$f = C\\sqrt{\\sigma / \\rho R^3}$ — every exponent forced, one constant $C$
left for experiment.

**(a) By hand:** reproduce the Buckingham argument. Variables $f, R, \\rho,
\\sigma$; count dimensions; form the single $\\Pi$ group; solve the exponents.

**(b) Verify:** hand the same four variables to `pi_groups` and check it
finds your group (exponent dicts are given below).

**(c) Apply:** a hydrophone in a pond hears a clean 43 Hz "bloop" from a
bubble measured at $R = 2.0$ mm ($\\sigma = 0.0728$ N/m, $\\rho = 998$
kg/m³). Calibrate $C$, then predict the radius of a bubble singing at
150 Hz. Bigger bubbles ring — which way? Does your formula agree with every
bell you've ever heard?
"""),
code("""
from fluidtools.similitude import pi_groups
import numpy as np

variables = {
    'f':     {'T': -1},
    'R':     'length',
    'rho':   'density',
    'sigma': {'M': 1, 'T': -2},   # N/m = kg/s^2
}
# your work below
...
"""),
md("""
## Problem 3 · A very large problem, at 1:25

A spillway must pass a design flood of $Q_p = 400\\ \\mathrm{m^3/s}$; the
prototype flow runs about 8 m/s at 3 m depth. You build a 1:25 Froude-scaled
model in the hydraulics lab.

**(a) By hand:** what model flow rate $Q_m$ do you order the lab pump for?
(Derive the $Q$ ratio from Froude scaling — velocity scales like
$\\sqrt{L_r}$, area like $L_r^2$.)

**(b) Verify:** `scale_report('froude', ...)` should agree with your ratio.

**(c) Apply — the honesty check:** Froude scaling silently abandons
Reynolds. Estimate the model's Re (velocity × depth / $\\nu$). Is it still
comfortably turbulent ($\\gtrsim 10^4$), so the friction physics is at least
the right *kind*? One sentence on why we match Fr and merely *police* Re,
rather than the other way around.
"""),
code("""
from fluidtools.similitude import model_scale, scale_report
import numpy as np

Qp, Vp, yp = 400.0, 8.0, 3.0
Lr = 1/25
nu = 1.0e-6
# your work below
...
"""),
])

# ----------------------------------------------------------------- HW 2
hw2 = nb("Homework 2 — Hydrostatics (Unit 3)", [
md("""
## Problem 1 · The gate and the hinge

A vertical rectangular gate ($b = 2.0$ m wide, $L = 3.0$ m tall) holds back
a freshwater reservoir; the gate's top edge sits 1.5 m below the surface.
It hinges along the top edge and rests on a stop at the bottom.

**(a) By hand:** resultant force on the gate, depth of the center of
pressure, and the reaction at the bottom stop (moment balance about the
hinge).

**(b) Verify:** build the `PlaneSurface`, call `gate_hinge_reaction`, and
reconcile every number with (a).

**(c) Apply:** the stop is corroding and the client wants a single support
beam placed so the *hinge carries no moment at all*. At what depth does the
beam go — and why is "at the centroid" the wrong answer by exactly the
margin your part (a) computed?
"""),
code("""
from fluidtools.hydrostat import PlaneSurface, gate_hinge_reaction

gate = PlaneSurface(shape='rectangle', b=2.0, L=3.0, y_top=1.5, theta_deg=90)
# your work below
...
"""),
md("""
## Problem 2 · Manometer walking tour

A pressurized air tank connects to an open U-tube: from the tank tap, the
line drops 0.40 m through oil (SG 0.85) to the mercury interface (SG 13.6);
the mercury leg rises 0.25 m to the open side. Take the walking tally from
lecture: *down adds $\\rho g h$, up subtracts, air contributes nothing
worth counting.*

**(a) By hand:** the tank's gauge pressure. Show the walk as a ledger, one
line per leg.

**(b) Verify:** the arithmetic cell below is your calculator, not your
brain — fill in the tally and compare.

**(c) Apply:** in Boulder (83 kPa atmosphere) the *gauge* reading is
unchanged. The *absolute* pressure is not. Compute both and state, in one
sentence, which one the tank's relief valve should care about and why.
""" ),
code("""
RHO_W, G = 998.2, 9.81
h_oil, SG_oil = 0.40, 0.85
h_hg,  SG_hg  = 0.25, 13.6
p_atm_boulder = 83e3
# your walking tally below
...
"""),
md("""
## Problem 3 · Barge economics

A rectangular barge — 12 m long, 5 m beam, 2 m depth, weight 300 kN, KG =
1.10 m — is about to take shipping containers *on deck*, each 20 kN with
its own center of gravity 2.5 m above the keel once loaded.

**(a) By hand:** the empty barge's draft, KB, BM ($= I/V$), and GM. Is it
stable, and generously so?

**(b) Verify:** the `Barge` class should reproduce your GM.

**(c) Apply:** company policy demands GM ≥ 0.30 m. Sweep the container
count and report the maximum legal load — so that the cargo, the barge,
and your professional license all stay dry. Where does the stability
actually go as you stack? (One sentence: which letter of K-B-M-G moves,
and which way.)
"""),
code("""
from fluidtools.hydrostat import Barge
import numpy as np

barge = Barge(length=12.0, beam=5.0, depth=2.0, weight=300e3, KG=1.10)
W_box, KG_box = 20e3, 2.5
# your work below
...
"""),
])

# ----------------------------------------------------------------- HW 3
hw3 = nb("Homework 3 — Bernoulli & Control Volumes (Units 4–5)", [
md("""
## Problem 1 · Two holes, one bet

A tall open tank holds water to $H = 1.20$ m above the table it sits on.
Your lab partner bets that a hole drilled at 0.30 m below the surface
out-squirts one drilled at 0.90 m below.

**(a) By hand:** Torricelli both jets; compute both landing ranges. Settle
the bet, and explain the tie in one sentence about the trade being made.

**(b) Verify:** sweep the hole depth in code, plot range vs. depth, and
confirm both the symmetry and the optimum at $H/2$.

**(c) Apply:** the client wants a jet that lands *exactly* 1.00 m out.
Report every hole depth that works. (There are two — say why, and which
you'd drill if the tank slowly drains during use.)
"""),
code("""
import numpy as np
import matplotlib.pyplot as plt

H = 1.20
g = 9.81
range_of = lambda h: 2*np.sqrt(h*(H-h))
# your work below
...
"""),
md("""
## Problem 2 · The venturi pays rent

A 75 mm water line carries 8.0 L/s. You must meter it, and the transducer
tops out at $\\Delta P = 20$ kPa.

**(a) By hand:** size the venturi throat $d$ from ideal Bernoulli + mass
(ignore $C_d$ for the first pass). 

**(b) Verify:** `size_meter` with the real $C_d = 0.98$ — how much did the
discharge coefficient move your throat diameter?

**(c) Apply — the rent:** compare the *permanent* loss of your venturi
against an orifice plate sized for the same 20 kPa reading
(`permanent_loss`). Convert each to pumping power at 8 L/s, then to dollars
over 4000 h/yr at \\$0.13/kWh. One sentence to the purchasing department,
who have noticed the orifice is cheaper to buy. For capitalism reasons,
show the five-year total.
"""),
code("""
from fluidtools.flowmeter import size_meter, meter_flow, permanent_loss
RHO = 998.2
Q, D, dP = 8.0e-3, 0.075, 20e3
# your work below
...
"""),
md("""
## Problem 3 · Brace yourself

A fire crew's 64 mm hose feeds a 25 mm nozzle at $Q = 15$ L/s, discharging
horizontally to atmosphere.

**(a) By hand:** exit velocity, then the momentum-flux estimate of the
force the crew must brace against. State the direction with a control-volume
sketch, not vibes.

**(b) Verify:** `nozzle_thrust` — and reconcile any difference with your
estimate (what did the pressure-area term contribute?).

**(c) Apply:** OSHA-of-the-imagination rates one braced firefighter at
600 N. How many firefighters is this hose? Round in the direction that
keeps everyone employed and upright.
"""),
code("""
from fluidtools.momentum import nozzle_thrust
import numpy as np

RHO = 998.2
Q, D_hose, D_exit = 15e-3, 0.064, 0.025
# your work below
...
"""),
])

# ----------------------------------------------------------------- HW 4
hw4 = nb("Homework 4 — Profiles & External Flow (Units 6–7)", [
md("""
## Problem 1 · Backflow detective

Oil ($\\mu = 0.05$ Pa·s) fills a 4.0 mm gap. The top wall slides at
$U = 0.30$ m/s; the bottom wall is fixed; a pressure gradient opposes the
wall's dragging.

**(a) By hand:** the adverse $dp/dx$ at which backflow *first* appears at
the bottom wall (the criterion is zero wall shear — derive
$dp/dx = 2\\mu U / h^2$ from the general Couette–Poiseuille profile).

**(b) Verify:** `backflow_threshold(h, U, mu)`, then plot the profile at
$1.5\\times$ threshold with `couette_poiseuille` and point at the reversed
region.

**(c) Apply:** your plant's coating line shows streaks consistent with
backflow when the die pressure spikes. The gap is fixed; the fluid is
fixed. What single knob remains, which way do you turn it, and by how much
for a 2× safety margin at the spike pressure?
"""),
code("""
from fluidtools.profiles import couette_poiseuille, backflow_threshold
import numpy as np
import matplotlib.pyplot as plt

h, U, mu = 4.0e-3, 0.30, 0.05
# your work below
...
"""),
md("""
## Problem 2 · The plate pays up

A 2.0 m × 0.5 m plate is towed edge-on through 20 °C water at 3.0 m/s,
wetted on both sides.

**(a) By hand:** $Re_L$; laminar, turbulent, or mixed? Use the mixed-regime
$C_F = 0.074\\,Re_L^{-1/5} - 1742/Re_L$ and estimate the drag.

**(b) Verify:** `plate_drag(U, L, width, rho, mu, sides=2)`.

**(c) Apply:** the client can afford to polish (and keep laminar) only the
leading portion. Using the tool, how much of the total drag comes from the
first 0.5 m as-is — and why does polishing the *front* buy more than the
same area at the back? (Two sentences: one for $\\tau_w(x)$, one for
transition.)
"""),
code("""
from fluidtools.extflow import plate_drag, cf_local, delta
RHO, MU = 998.2, 1.0e-3
U, L, w = 3.0, 2.0, 0.5
# your work below
...
"""),
md("""
## Problem 3 · The crisis pitch, and hail

**(a) By hand:** a baseball ($D = 73$ mm) is thrown at 40 m/s through air
($\\rho = 1.20$, $\\mu = 1.81\\times10^{-5}$). Compute Re. How close is it
to the smooth-sphere drag crisis ($\\sim 3\\times10^5$)?

**(b) Verify:** `cd_sphere` just below and just above the crisis — report
the drop, and one sentence on what the seams of a non-spinning knuckleball
are doing in this neighborhood.

**(c) Apply:** hailstones ($\\rho_p \\approx 900$ kg/m³) of 2 cm and 4 cm.
`terminal_velocity_sphere` both; then compare their *kinetic energies* on
arrival. Doubling the diameter multiplies the damage by what factor — and
which two scalings conspired to do it?
"""),
code("""
from fluidtools.extflow import cd_sphere, terminal_velocity_sphere
import numpy as np

RHO_AIR, MU_AIR = 1.20, 1.81e-5
# your work below
...
"""),
])

# ----------------------------------------------------------------- HW 5
hw5 = nb("Homework 5 — Pipe Flow, Iteration, and the Pump Aisle (Unit 8)", [
md("""
## Problem 1 · The hydrant (legally dubious, for capitalism reasons)

Your roommate has located a fire hydrant, a 30 m rental hose (smooth,
$D = 65$ mm), and an above-ground pool at the same elevation as the
hydrant. The hydrant supplies a steady **500 kPa gauge**. The hose
discharges freely into the pool (exit loss $K_L = 1.0$; neglect the
hydrant's own losses; 20 °C water).

**(a) By hand:** the head equation gives
$\\dfrac{p}{\\rho g} = \\left(f\\dfrac{L}{D} + K_L\\right)\\dfrac{U^2}{2g}$ —
but $f$ needs $Re$ needs $U$. Build the guess-and-check table from lecture:
start at $f = 0.020$, compute $U$, update $Re$, update $f$ (smooth-pipe
Colebrook), repeat until $f$ stops moving. **Three rows minimum, shown.**

**(b) Verify:** `PipeSystem` + `solve_flow()` in three lines. Your table's
final row and the solver should agree to the liter — they are the *same
loop*.

**(c) Apply:** the city (hypothetically) caps unmetered draw at 6.0 L/s.
Size the globe-valve loss coefficient $K_L$ that brings your rig into
compliance, and report the valve's share of the total head as a percentage.
One sentence on where those meters of head went.
"""),
code("""
from fluidtools.pipeflow import *
import numpy as np

w = water(20)
p_supply = 500e3   # Pa gauge
L, D = 30.0, 0.065
# (a) your iteration table (a loop or by hand — but show the rows)
...
"""),
code("""
# (b) the three lines
sys = PipeSystem(w, z1=0, z2=0, p1=500e3, p2=0)
pipe = Pipe(L=30.0, D=0.065, roughness=0.0)
pipe.add_fitting(K=1.0)
sys.add(pipe)
# your work below (solve, print, compare to your table)
...
"""),
md("""
## Problem 2 · Loss anatomy: the shower audit

A dorm shower run: 10 m of 12 mm copper (smooth), a globe valve
($K_L = 10$), four elbows ($K_L = 1.5$ each), a tee ($K_L = 2.0$), and the
showerhead itself ($K_L = 8$, engineered restriction). Supply head: 25 m.

**(a) By hand:** at a guessed 8 L/min, what fraction of the loss is pipe
friction vs. fittings? (You may take $f \\approx 0.03$ for the estimate.)

**(b) Verify:** solve the real flow with the toolbox and recompute the
split at the actual operating point. Did the *ratio* care about the flow
rate? Say why in one sentence (both losses scale how?).

**(c) Apply:** facilities will fund exactly one change: (i) repipe at
16 mm, (ii) halve the pipe length, or (iii) delete two elbows. Rank by
delivered flow gain, with numbers. The winner is not close — explain its
unfair advantage ($D$ appears where, at what power?).
"""),
code("""
from fluidtools.pipeflow import *
w = water(40)
# build, solve, audit — your work below
...
"""),
md("""
## Problem 3 · The pump aisle (choose wisely)

Your pump must lift 20 °C water from a sump ($z_1 = 0$) to a tank at
$z_2 = 12$ m: suction side 8 m of 80 mm galvanized pipe (sharp entrance +
strainer), discharge side 170 m of the same pipe (three 90° elbows, a swing
check valve, an exit). The duty: **at least 10.5 L/s**. Three candidates,
all mounted with their eye 3.5 m above the sump surface:

| Pump | curve $H(Q)$ (m, $Q$ in m³/s) | $\\eta$ | NPSH$_r$(Q) (m) | price |
|---|---|---|---|---|
| A | $50 - 9.0\\times10^4\\,Q^2$ | 0.62 | $2.0 + 3.0\\times10^4 Q^2$ | \\$3800 |
| B | $35 - 6.0\\times10^4\\,Q^2$ | 0.55 | $1.2 + 2.0\\times10^4 Q^2$ | \\$2100 |
| C | $38 - 7.0\\times10^4\\,Q^2$ | 0.78 | $1.2 + 1.0\\times10^4 Q^2$ | \\$3200 |

**(a) By hand:** *before computing*, write one sentence per pump predicting
its weakness from the table alone (shutoff head? efficiency? appetite?).

**(b) Verify & apply:** for each pump — operating point, duty check, NPSH
margin (`check_cavitation`), shaft power, and five-year cost
(price + 5 yr × 4000 h/yr × \\$0.13/kWh). Present a table.

**(c) The memo:** one paragraph to your boss naming the pump you'd buy.
The cheapest sticker is a trap and the biggest pump is a different trap —
your memo should spring both, with numbers. (A pump that cavitates is not
a pump; it is a very expensive maraca.)
"""),
code("""
from fluidtools.pipeflow import *

def build_system(pump, T=20):
    fl = water(T)
    suction = Pipe(L=8, D=0.08, roughness='galvanized_iron')
    suction.add_fitting('entrance_sharp').add_fitting('strainer')
    discharge = Pipe(L=170, D=0.08, roughness='galvanized_iron')
    discharge.add_fitting('elbow_90_standard', n=3)
    discharge.add_fitting('check_valve_swing').add_fitting('exit')
    s = PipeSystem(fl, z1=0, z2=12)
    s.add(suction, pump, discharge)
    return s

pumps = {
  'A': Pump.from_coeffs(50, 9.0e4, efficiency=0.62, z=3.5, npsh_r=lambda Q: 2.0 + 3.0e4*Q**2),
  'B': Pump.from_coeffs(35, 6.0e4, efficiency=0.55, z=3.5, npsh_r=lambda Q: 1.2 + 2.0e4*Q**2),
  'C': Pump.from_coeffs(38, 7.0e4, efficiency=0.78, z=3.5, npsh_r=lambda Q: 1.2 + 1.0e4*Q**2),
}
prices = {'A': 3800, 'B': 2100, 'C': 3200}
# your work below
...
"""),
md("""
## Problem 4 · Throttle vs. VFD (the utility bill)

The winning pump from Problem 3 turns out to need only **8.0 L/s** for
3000 of its 4000 annual hours (full duty the rest).

**(a) By hand:** predict which control strategy wins before computing.
Affinity says power scales like *what* with speed?

**(b) Verify & apply:** achieve 8.0 L/s two ways — (i) throttle: add
discharge $K_L$ until $Q = 8.0$; (ii) slow down: `pump.at_speed(k)` until
$Q = 8.0$. Shaft power for each; annual energy cost difference at
\\$0.13/kWh for the 3000 reduced-duty hours.

**(c)** The VFD retrofit quote is \\$1900 installed. Payback period, one
sentence, and whether you sign.
"""),
code("""
from fluidtools.pipeflow import *
import numpy as np
# reuse build_system and pump C from Problem 3
# your work below
...
"""),
])

# ----------------------------------------------------------------- HW 6
hw6 = nb("Homework 6 — Channels & Turbomachinery (Units 9–10)", [
md("""
## Problem 1 · The stilling basin

Flow leaves a spillway toe at $y_1 = 0.35$ m moving $U_1 = 14$ m/s, across
a 40 m wide apron.

**(a) By hand:** $Fr_1$, the conjugate depth $y_2$ (Bélanger), the energy
destroyed per unit weight $\\Delta E = (y_2-y_1)^3/(4y_1y_2)$, and the
*fraction* of the arriving specific energy retired.

**(b) Verify:** `Channel.conjugate_depth` and `jump_energy_loss`.

**(c) Apply:** total power dissipated in the basin
($P = \\rho g\\, Q\\, \\Delta E$, with $Q$ from the full width). Express it
in MW and in "electric kettles running simultaneously" (1.5 kW each). This
is the one chapter of the course where that number is the *product*, not
the waste — say so to the review board in one sentence.
"""),
code("""
from fluidtools.channel import Channel
import numpy as np

y1, U1, width = 0.35, 14.0, 40.0
g, rho = 9.81, 998.2
ch = Channel(b=width, m=0.0)
# your work below
...
"""),
md("""
## Problem 2 · The canal that needed mowing

An earth irrigation canal ($n = 0.022$), rectangular, $b = 4.0$ m, slope
$S = 0.0008$, carries $Q = 9.0$ m³/s.

**(a) By hand:** set up (don't fully iterate) the Manning normal-depth
equation, then let `normal_depth` finish the arithmetic. Also get $y_c$,
$Fr$, and the mild/steep verdict.

**(b) Apply:** the district lines the canal with finished concrete
($n = 0.012$). New normal depth, and the freeboard recovered.

**(c)** The unlined canal's weeds regrow toward $n = 0.035$ by late summer.
At what $n$ does the flow overtop the 2.0 m banks? One sentence connecting
the mowing budget to the flood insurance premium.
"""),
code("""
from fluidtools.channel import Channel, MANNING_N
import numpy as np

Q, b, S = 9.0, 4.0, 0.0008
# your work below
...
"""),
md("""
## Problem 3 · One measured curve, every speed

A pump measured at 1450 rpm fits $H = 32 - 5.0\\times10^4\\,Q^2$
($\\eta = 0.72$). It feeds a tank at $z_2 = 10$ m through 40 m of 100 mm
commercial-steel pipe.

**(a) By hand:** sketch (actually sketch) the two curves and mark the
operating point you expect. Then solve it with the toolbox.

**(b) Apply:** demand rises 30%. Using `at_speed`, find the speed ratio —
and the rpm — that delivers it. (Bracket by hand first: affinity says $Q$
scales like what?)

**(c)** Compare the shaft power at the two speeds to the cube-law
prediction $k^3$. It won't match exactly — one sentence on what the cube
law assumes about the *system* that your tank-lift system violates.
"""),
code("""
from fluidtools.pipeflow import *
import numpy as np

w = water(20)
pump = Pump.from_coeffs(32, 5.0e4, efficiency=0.72)
# your work below
...
"""),
md("""
## Problem 4 · Machine safari

Three duties, all real categories:

1. boiler feed: $Q = 0.02$ m³/s against $H = 180$ m at 2900 rpm
2. irrigation lift: $Q = 3.0$ m³/s against $H = 4$ m at 300 rpm
3. municipal supply: $Q = 0.40$ m³/s against $H = 30$ m at 1150 rpm

**(a) By hand:** compute each specific speed
$N_s = \\omega\\sqrt{Q}/(gH)^{3/4}$ (radians!).

**(b) Verify:** `specific_speed` + `classify_pump` (or let
`select_machine` narrate).

**(c)** One sentence each: what does the *shape* of each machine look
like, and what goes wrong if you install duty 2's machine on duty 1's
job?
"""),
code("""
from fluidtools.turbo import specific_speed, classify_pump, select_machine, rpm_to_rad
# your work below
...
"""),
md("""
## Problem 5 · The wheel and the cart (last word)

An alpine penstock delivers $H = 500$ m at $Q = 0.15$ m³/s to a single
Pelton nozzle.

**(a) By hand:** jet speed $V$, the wheel's ideal maximum power (buckets
at $V/2$), and the lone-cart maximum for comparison (cart at $V/3$,
$P = \\tfrac{8}{27}\\cdot\\tfrac12\\rho A V^3$… or just evaluate
$\\rho A (V-U)^2 U$ at its optimum).

**(b) Verify:** `pelton_optimum(V, series=True)` vs `series=False`.

**(c)** One sentence: the single assumption that separates the two optima
— and why a wheel gets to make it while a cart doesn't.
"""),
code("""
from fluidtools.momentum import pelton_optimum
import numpy as np

rho, g = 998.2, 9.81
H, Q = 500.0, 0.15
# your work below
...
"""),
])

import os
os.makedirs('homework', exist_ok=True)
for name, notebook in [('hw1-properties-dimensions', hw1), ('hw2-hydrostatics', hw2),
                       ('hw3-bernoulli-cv', hw3), ('hw4-profiles-external', hw4),
                       ('hw5-pipeflow-pumps', hw5), ('hw6-channels-turbo', hw6)]:
    nbf.write(notebook, f'homework/{name}.ipynb')
    print(f"wrote homework/{name}.ipynb ({len(notebook.cells)} cells)")
