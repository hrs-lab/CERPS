"""Python reproduction of Huang, Situ & Lei (2024) J. Theor. Biol. 577:111664.

"Dynamics of cell-type transition mediated by epigenetic modifications".
Re-implementation (numba-accelerated) of the hybrid model released by the
authors at https://github.com/jinzhilei/Cell-type-transition (GPL-3.0).
"""
from .params import Params, DEFAULT
from .sim import simulate

__all__ = ["Params", "DEFAULT", "simulate"]
