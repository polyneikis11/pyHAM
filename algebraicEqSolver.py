from sympy.abc import *
import hamTask
from deformationEquations import DeformEqAlgebr
from hamSymbols import x

ALGEBRAIC_SOLVER_TOLERANCE = 0.00001


def is_zero(expr):
    return abs(float(expr)) < ALGEBRAIC_SOLVER_TOLERANCE


def check_if_solves_eq(eq, expr):
    print('expr is', expr)
    expr_eval = eq.evalf(subs={x: expr})
    print('expr eval', expr_eval)
    print(expr_eval.diff(x))
    ans = is_zero(expr_eval)
    return ans


def eval_diff(eq, expr):
    return eq.evalf(subs={x: expr})


def valid_initial_guess(eq, x0):
    first_der = eq.diff(x)
    return not check_if_solves_eq(first_der, x0)


def eval_Q(var, x0):
    if var == 0:
        return x0
    else:
        return -1


def solve_alg_eq(expr, x0, alg_subproblem):
    c0 = alg_subproblem.c0
    print('will solve equation {} with initial guess {} and control param {}'.format(expr, x0, c0))
    approx_num = alg_subproblem.approximation
    solver = DeformEqAlgebr(c0, expr, x0)
    alg_subproblem.answer.value = solver.calc_solution(approx_num)
    alg_subproblem.answer.error = eval_diff(expr, alg_subproblem.answer.value)
    alg_subproblem.answer.status = hamTask.SubSolutionStatus.SOLVED
    print('solution is: ', alg_subproblem.answer)


def solve_kernel_algebraic(alg_subproblem):
    print('Would solve Now subproblem')
    print(alg_subproblem)
    print('-'*20)
    alg_subproblem.answer = hamTask.SolutionAnswer()
    expr = alg_subproblem.task_parent.eq
    x0 = alg_subproblem.initial_guess
    if check_if_solves_eq(expr, x0):
        alg_subproblem.answer.status = hamTask.SubSolutionStatus.SOLVED
        alg_subproblem.answer.value = x0
        alg_subproblem.answer.error = eval_diff(expr, x0)
    else:
        if valid_initial_guess(expr, x0):
            solve_alg_eq(expr, x0, alg_subproblem)
        else:
            alg_subproblem.answer.status = hamTask.SubSolutionStatus.BAD_INITIAL_GUESS

