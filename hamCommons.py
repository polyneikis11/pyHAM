import sys
from enum import Enum
from sympy import simplify
from hamSymbols import *

pyham_bug_report_email = 'pyham@gmail.com'
flt_tolerance = 1e-20


class HamSolveTypes(Enum):
    UNKNOWN_TYPE = 0
    ALGEBRAIC_EQUATION = 1
    ODE_BVP = 2


def init_expr():
    return S(0.)


def h_sympify(expr):
    return ham_sympify(expr)
    # return sympify(expr)


def is_expr_valid(expr):
    return expr is not None


def is_obj_valid(obj):
    return obj is not None


def is_ode_expr(expr):
    ans = False
    symbols_list = expr.atoms(Symbol)
    for current_symbol in symbols_list:
        if current_symbol in ODE_F_EXPRESSIONS:
            ans = True
            break
    return ans


def abort_program(error_message=None):
    if is_obj_valid(error_message):
        print(error_message)
    sys.exit(1)


def is_algebraic_expr(expr):
    ans = False
    if not is_ode_expr(expr):
        diff_result = expr.diff(x)
        if simplify(diff_result) != S(0):
            ans = True
    return ans


def is_same_point_1d_flt(point_1, point_2):
    return abs(point_1 - point_2) < flt_tolerance


class HamErrorCode(Enum):
    HAM_OK = 0
    HAM_UNKNOWN_EQUATION_TYPE = 1
    HAM_ODE_WRONG_NUMBER_OF_BCS = 2
    HAM_ODE_WRONG_LINEAR_OPERATOR = 3
    HAM_ODE_INVALID_INITIAL_GUESS = 4
    HAM_ALGEBRAIC_INVALID_INITIAL_GUESS = 5


def ham_error_code_to_msg(error_code):
    msg = 'HAM Error: '
    if HamErrorCode.HAM_UNKNOWN_EQUATION_TYPE == error_code:
        msg += 'Unknown type of equation'
    elif HamErrorCode.HAM_ODE_WRONG_NUMBER_OF_BCS == error_code:
        msg += 'Valid boundary conditions do not match the order of equation'
    elif HamErrorCode.HAM_ODE_WRONG_LINEAR_OPERATOR == error_code:
        msg += 'Linear operator is invalid for the equation to be solved.'
    elif HamErrorCode.HAM_ODE_INVALID_INITIAL_GUESS == error_code:
        msg += 'Initial guess is invalid for the ODE to be solved'
    elif HamErrorCode.HAM_ALGEBRAIC_INVALID_INITIAL_GUESS == error_code:
        msg += 'Initial guess is invalid for the algebraic equation to be solved'
    return msg
