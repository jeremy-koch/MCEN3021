"""fluidtools: the MCEN 3021 Fluid Mechanics toolbox suite.

Submodules (import what the assignment needs):
    fluidtools.hydrostat   pressure fields, gate forces, buoyancy, stability
    fluidtools.momentum    control-volume momentum: jets, bends, thrust
    fluidtools.similitude  dimensional analysis (pi groups), model scaling
    fluidtools.pipeflow    pipe networks, losses, pumps, cavitation (v1 API)
    fluidtools.profiles    viscous flow profiles: read shear & pressure off u(y)
    fluidtools.extflow     boundary layers, drag build-up, terminal velocity
    fluidtools.flowmeter   orifice / venturi / nozzle / pitot sizing
    fluidtools.channel     open-channel flow: Manning, jumps, profiles
    fluidtools.turbo       specific speed, Euler head, machine selection

Convenience: `from fluidtools.pipeflow import *` keeps every notebook written
against the v1 `pipeflow` package working after a one-word import change.
"""
from . import hydrostat, momentum, similitude, pipeflow, profiles, extflow, flowmeter, channel, turbo

G = 9.81
__version__ = "2.1.0"
