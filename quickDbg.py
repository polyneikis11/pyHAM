from sympy import lambdify, LambertW, exp, simplify, log, Symbol, dsolve
from sympy import integrate, cos, Integral, pdsolve, sin, pi, erf, sqrt
from hamSymbols import *


def quick_dbg():
    # expr = -(sqrt(pi*p1/2)*exp(p2**2/(2*p1))*erf((x - p2)/sqrt(2*p1)))
    expr = t**(x/4)
    print('diff: ', expr.diff(x))
    assert False
