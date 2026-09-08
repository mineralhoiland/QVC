"""
qvc.dynamics — PDE/ODE solvers for TEVC time evolution.

Tier 2 extension: implements the dynamical equations of QVC.
- TQGL solver (split-step Fourier method)
- DCE self-limiting ODE
- Lindblad dissipator for open-system dynamics
- Schwinger pair production rate
"""

from qvc.dynamics.tqgl_solver import TQGLSolver
