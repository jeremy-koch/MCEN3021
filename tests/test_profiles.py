import numpy as np
from fluidtools import profiles as pr

h, U, mu = 0.01, 1.0, 1.0e-3

def test_pure_couette_is_linear():
    y = np.linspace(0, h, 7)
    u = pr.couette_poiseuille(y, h, U=U, dpdx=0.0, mu=mu)
    assert np.allclose(u, U * y / h)
    tau = pr.shear_couette_poiseuille(y, h, U=U, dpdx=0.0, mu=mu)
    assert np.allclose(tau, mu * U / h)     # uniform shear

def test_pure_poiseuille_symmetry_and_max():
    dpdx = -2000.0
    u_mid = pr.couette_poiseuille(h/2, h, U=0, dpdx=dpdx, mu=mu)
    assert abs(u_mid - (-dpdx) * h**2 / (8 * mu)) < 1e-12
    tau_mid = pr.shear_couette_poiseuille(h/2, h, U=0, dpdx=dpdx, mu=mu)
    assert abs(tau_mid) < 1e-12             # zero shear at centerline
    assert abs(pr.shear_couette_poiseuille(0, h, 0, dpdx, mu)) == abs(dpdx)*h/2

def test_backflow_onset():
    crit = pr.backflow_threshold(h, U, mu)   # = 2 mu U / h^2 = 20 Pa/m
    assert abs(crit - 20.0) < 1e-12
    y = np.linspace(1e-5, h, 400)
    assert (pr.couette_poiseuille(y, h, U, 0.99*crit, mu) > 0).all()
    assert (pr.couette_poiseuille(y, h, U, 1.5*crit, mu) < 0).any()
    assert not pr.has_backflow(h, U, 0.99*crit, mu)
    assert pr.has_backflow(h, U, 1.5*crit, mu)

def test_channel_flow_rate_matches_integral():
    dpdx = 800.0
    q = pr.flow_rate_couette_poiseuille(h, U, dpdx, mu)
    y = np.linspace(0, h, 20001)
    q_num = np.trapezoid(pr.couette_poiseuille(y, h, U, dpdx, mu), y)
    assert abs(q - q_num) / abs(q_num) < 1e-6

def test_pipe_poiseuille_Q():
    R, dpdx = 0.02, -500.0
    r = np.linspace(0, R, 20001)
    Q_num = np.trapezoid(pr.pipe_poiseuille(r, R, dpdx, mu) * 2*np.pi*r, r)
    assert abs(pr.pipe_poiseuille_Q(R, dpdx, mu) - Q_num) / Q_num < 1e-6

def test_film_free_surface_shear_zero():
    hf = 0.002
    y = np.array([hf - 1e-7, hf])
    u = pr.film_on_wall(y, hf, mu=mu)
    assert abs(u[1] - u[0]) / u[1] < 1e-4    # flat at the surface
    q_num = np.trapezoid(pr.film_on_wall(np.linspace(0, hf, 20001), hf, mu), 
                         np.linspace(0, hf, 20001))
    assert abs(pr.film_flow_rate(hf, mu) - q_num) / q_num < 1e-6

def test_wall_shear_from_data_recovers_analytic():
    dpdx = 1500.0
    y = np.linspace(0, 0.2*h, 9)             # near-wall samples
    u = pr.couette_poiseuille(y, h, U, dpdx, mu)
    tau_true = pr.shear_couette_poiseuille(0, h, U, dpdx, mu)
    tau_fit, _ = pr.wall_shear_from_data(y, u, mu, deg=2)
    assert abs(tau_fit - tau_true) / abs(tau_true) < 1e-6
