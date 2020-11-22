from sympy import simplify, lambdify
import hamTask
from hamSymbols import *
from hamErrors import ham_error
from hamBase import ham_not_none
from deformationEquations import DeformEqPDE
from geometry import HamPoint
from hamPlotter import HamPlotter


class PdeErrorCalculator:
    def __init__(self):
        self.deq = None
        self.domain = None
        self.eval_points = []

    def _subs_solution(self, solution):
        the_symb = self.deq.atoms(Symbol)
        gen = (var for var in the_symb if var in PDE_U_EXPRESSIONS)
        subs_lexicon = {}
        for u_fun in gen:
            x_ord, y_ord, t_ord = decode_der_orders_pde(u_fun)
            res_expr = solution.diff(x, x_ord)
            res_expr = res_expr.diff(y, y_ord)
            res_expr = res_expr.diff(t, t_ord)
            subs_lexicon[u_fun] = res_expr
        ret_expr = self.deq.subs(subs_lexicon)
        return simplify(ret_expr)

    def _consistent_points(self):
        ret_vars = None
        x_absent = False
        y_absent = False
        t_absent = False
        falsy = False
        for pnt in self.eval_points:
            if not ham_not_none(pnt.x):
                if x_absent:
                    falsy = True
                    break
                else:
                    x_absent = True
            if not ham_not_none(pnt.y):
                if y_absent:
                    falsy = True
                    break
            if not ham_not_none(pnt.t):
                if t_absent:
                    falsy = True
                    break
        if not falsy:
            ret_vars = []
            if not x_absent:
                ret_vars.append(x)
            if not y_absent:
                ret_vars.append(y)
            if not t_absent:
                ret_vars.append(t)
        return ret_vars

    def calc(self, solution, bcs):
        error_value = -1.
        if ham_not_none(self.domain):
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__NOT_IMPLEMENTED_YET')
        elif len(self.eval_points) > 0:
            error_value = 0.
            the_vars = self._consistent_points()
            if ham_not_none(the_vars):
                expr_to_eval = self._subs_solution(solution)
                expr_eval_l = lambdify(the_vars, expr_to_eval, modules=['math', 'sympy'])
                if len(self.eval_points) > 0:
                    for pnt in self.eval_points:  # type: HamPoint
                        get_back_vars = pnt.p_list(the_vars)
                        error_value += expr_eval_l(get_back_vars[0], get_back_vars[1], get_back_vars[2])
                    error_value /= len(self.eval_points)
                else:
                    error_value = 'N/A'
            else:
                ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        else:
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        return error_value


def expr_solves_pdeq(expr, eq):
    subs_lexicon = {}
    for i in range(len(PDE_U_EXPRESSIONS)):
        curr_u_expr = PDE_U_EXPRESSIONS[i]
        x_ord, y_ord, t_ord = decode_der_orders_pde(curr_u_expr)
        temp_expr = expr.diff(x, x_ord)
        temp_expr = temp_expr.diff(y, y_ord)
        temp_expr = temp_expr.diff(t, t_ord)
        subs_lexicon[curr_u_expr] = temp_expr
    res = eq.subs(subs_lexicon)
    return S(0) == simplify(res)


def check_if_solves_pde_bvp(eq, bcs, expr):
    ans = True
    for curr_bcs in bcs:
        if not curr_bcs.satisfied_by(expr):
            ans = False
            break
    if ans:
        ans = expr_solves_pdeq(expr, eq)
    return ans


def solve_pde_bvp(pde_subproblem):
    pde_to_use = pde_subproblem  # type: hamTask.HamTaskSubSolution
    force = pde_to_use.force
    pdeq = pde_to_use.task_parent.eq
    bcs = pde_to_use.bcs
    x0 = pde_to_use.initial_guess
    lin_op = pde_to_use.lin_op
    approx = pde_to_use.approximation
    c0 = pde_to_use.c0
    error_calculator = pde_to_use.task_parent.error_calculator
    error_calculator.deq = pdeq
    solver = DeformEqPDE(pdeq, bcs, x0, lin_op, c0)
    if solver.check_bcs():
        if force or solver.valid_x0_for_bcs():
            pde_to_use.answer.value, sympy_solved = solver.calc_solution(approx)
            if sympy_solved:
                pde_to_use.answer.error = error_calculator.calc(pde_to_use.answer.value, bcs)
                print("SOLUTION: ", pde_to_use.answer.value)
                ham_plotter = HamPlotter()
                ham_plotter.do_the_plot(pde_to_use.answer.value, error_calculator.eval_points)
                print('TOTAL ERROR IS ', pde_to_use.answer.error)
                # TODO plot here etc...
                pde_to_use.answer.status = hamTask.SubSolutionStatus.SOLVED
            else:
                pde_to_use.answer.error = None
                pde_to_use.answer.status = hamTask.SubSolutionStatus.SYMPY_FAILED_FALL_BACK
        else:
            pde_to_use.answer.error = None
            pde_to_use.answer.value = None
            pde_to_use.answer.status = hamTask.SubSolutionStatus.BAD_INITIAL_GUESS
            print('BAD INITIAL GUESS')
    else:
        pde_to_use.error = None
        pde_to_use.answer.value = None
        pde_to_use.answer.status = hamTask.SubSolutionStatus.NOT_SUFFICIENT_BCS


def solve_kernel_PDE(pde_subproblem):
    print("------ PROBLEM TO SOLVE ---------")
    print(pde_subproblem)
    eq = pde_subproblem.task_parent.eq
    bcs = pde_subproblem.bcs
    x0 = pde_subproblem.initial_guess
    pde_subproblem_answer = hamTask.SolutionAnswer()
    pde_subproblem.answer = pde_subproblem_answer
    if check_if_solves_pde_bvp(eq, bcs, x0):
        pde_subproblem_answer.status = hamTask.SubSolutionStatus.SOLVED
        pde_subproblem_answer.value = x0
        # pde_subproblem_answer.error = eval_pde_bvp_ans_error()
    else:
        solve_pde_bvp(pde_subproblem)
