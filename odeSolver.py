import numpy as np
from scipy.special import lambertw
from matplotlib import pyplot as plt
from sympy import S, simplify, dsolve, Equality, oo, lambdify
from sympy import LambertW, solveset, integrate, exp
import hamTask
from hamSymbols import *
from hamErrors import ham_error
from deformationEquations import DeformEqODE


def expr_satisfies_bcs(expr, bcs_in):
    bcs = bcs_in  # type: hamTask.BoundaryCondition
    return bcs.satisfied_by_expr(expr)


def sympy_dsolve_sol(eq):
    subs_lexicon = {}
    for i in range(len(ODE_F_EXPRESSIONS)):
        curr_fder_symbol = ODE_F_EXPRESSIONS[i]
        if curr_fder_symbol in eq.free_symbols:
            subs_lexicon[curr_fder_symbol] = f(x).diff(x, i)
    eq_to_solve = eq.subs(subs_lexicon)
    ans = dsolve(eq_to_solve, f(x))
    if ans is not None and isinstance(ans, Equality):
        ret_ans = ans.rhs
    else:
        ret_ans = None
    return ret_ans


def expr_solves_deq(expr, eq):
    subs_lexicon = {}
    for i in range(len(ODE_F_EXPRESSIONS)):
        curr_fder_symbol = ODE_F_EXPRESSIONS[i]
        if curr_fder_symbol in eq.free_symbols:
            subs_lexicon[curr_fder_symbol] = expr.diff(x, i)
    res = eq.subs(subs_lexicon)
    return S(0) == simplify(res)


def check_if_solves_ode_bvp(eq, bcs, expr):
    ans = True
    for curr_bcs in bcs:
        if not expr_satisfies_bcs(expr, curr_bcs):
            ans = False
            break
    if ans:
        ans = expr_solves_deq(expr, eq)
    return ans


def eval_ode_bvp_ans_error(eq, bcs, expr):
    error = S(0)
    pass
    return error


def _get_terms(terms_list, number):
    ret_list = []
    for i in range(number + 1):
        ret_list.append(terms_list[i])
    return ret_list


def plot_to_present_a(solver, err_calc, analytical):
    # print('analytical solution: ', analytical)
    print('u0: ', solver.homotopy_terms[0])
    print('u1: ', solver.homotopy_terms[1])
    # print('u2: ', solver.homotopy_terms[2])
    # print('u3: ', solver.homotopy_terms[3])
    # print('u4: ', solver.homotopy_terms[4])
    # print('u5: ', solver.homotopy_terms[5])
    # print('u6: ', solver.homotopy_terms[6])
    # print('u7: ', solver.homotopy_terms[7])
    # print('u8: ', solver.homotopy_terms[8])
    # print('u9: ', solver.homotopy_terms[9])
    print('domain start: ', err_calc.domain_start)
    C0_num = -1.
    print('domain end  : ', err_calc.domain_end)
    the_vals = {p1: 1e3, p2: 21e7, p3: 1., C0: -0.5}
    # analytical = analytical.subs(the_vals)
    m1_sol = simplify((sum(_get_terms(solver.homotopy_terms, 1))).subs(the_vals))
    # m2_sol = simplify((sum(_get_terms(solver.homotopy_terms, 2))).subs({p1: 1., C0: C0_num}))
    # m3_sol = simplify((sum(_get_terms(solver.homotopy_terms, 3))).subs({p1: 1., C0: C0_num}))
    # m4_sol = simplify((sum(_get_terms(solver.homotopy_terms, 4))).subs({p1: 1., C0: C0_num}))
    # m5_sol = simplify((sum(_get_terms(solver.homotopy_terms, 5))).subs({p1: 1., C0: C0_num}))
    # m6_sol = simplify((sum(_get_terms(solver.homotopy_terms, 6))).subs({p1: 1., C0: C0_num}))
    # m7_sol = simplify((sum(_get_terms(solver.homotopy_terms, 7))).subs({p1: 1., C0: C0_num}))
    # m8_sol = simplify((sum(_get_terms(solver.homotopy_terms, 8))).subs({p1: 1., C0: C0_num}))
    # m9_sol = simplify((sum(_get_terms(solver.homotopy_terms, 9))).subs({p1: 1., C0: C0_num}))
    sample_pnts = []
    # anal_vals = []
    x0_vals = []
    m1_vals = []
    # m2_vals = []
    # m3_vals = []
    # m4_vals = []
    # m5_vals = []
    # m6_vals = []
    # m7_vals = []
    # m8_vals = []
    # m9_vals = []
    num_of_pnts = 1000
    next_point = err_calc.domain_start
    step = (err_calc.domain_end - err_calc.domain_start)/(num_of_pnts + 1)
    # anal_sol_l = lambdify(x, analytical, modules=['math', 'sympy'])
    x0_l = lambdify(x, solver.x0, modules=['math', 'sympy'])
    m1_l = lambdify(x, m1_sol, modules=['math', 'sympy'])
    # m2_l = lambdify(x, m2_sol, modules=['math', 'sympy'])
    # m3_l = lambdify(x, m3_sol, modules=['math', 'sympy'])
    # m4_l = lambdify(x, m4_sol, modules=['math', 'sympy'])
    # m5_l = lambdify(x, m5_sol, modules=['math', 'sympy'])
    # m6_l = lambdify(x, m6_sol, modules=['math', 'sympy'])
    # m7_l = lambdify(x, m7_sol, modules=['math', 'sympy'])
    # m8_l = lambdify(x, m8_sol, modules=['math', 'sympy'])
    # m9_l = lambdify(x, m9_sol, modules=['math', 'sympy'])
    # TODO test also the following if produces the same
    # anal_sol_l = lambdify(x, analytical, modules=['numpy', {'LambertW': lambertw}])
    for i in range(num_of_pnts):
        sample_pnts.append(next_point)
        # anal_vals.append(anal_sol_l(next_point))
        x0_vals.append(x0_l(next_point))
        m1_vals.append(m1_l(next_point))
        # m2_vals.append(m2_l(next_point))
        # m3_vals.append(m3_l(next_point))
        # m4_vals.append(m4_l(next_point))
        # m5_vals.append(m5_l(next_point))
        # m6_vals.append(m6_l(next_point))
        # m7_vals.append(m7_l(next_point))
        # m8_vals.append(m8_l(next_point))
        # m9_vals.append(m9_l(next_point))
        next_point += step
    sample_pnts = np.array(sample_pnts)
    # anal_vals = np.array(anal_vals)
    # plt.plot(sample_pnts, anal_vals, color='red', label='Analytical solution')
    plt.plot(sample_pnts, x0_vals, color='black', label='Initial guess')
    plt.plot(sample_pnts, m1_vals, color='red', label='HAM_1st_Order')
    # plt.plot(sample_pnts, m2_vals, color='rosybrown', label='HAM_2nd_Order')
    # plt.plot(sample_pnts, m3_vals, color='saddlebrown', label='HAM_3rd_Order')
    # plt.plot(sample_pnts, m4_vals, color='orange', label='HAM_4th_Order')
    # plt.plot(sample_pnts, m5_vals, color='yellow', label='HAM_5th_Order')
    # plt.plot(sample_pnts, m6_vals, color='lime', label='HAM_6th_Order')
    # plt.plot(sample_pnts, m7_vals, color='aqua', label='HAM_7th_Order')
    # plt.plot(sample_pnts, m8_vals, color='dodgerblue', label='HAM_8th_Order')
    # plt.plot(sample_pnts, m9_vals, color='navy', label='HAM_9th_Order')
    axes = plt.gca()
    plt.title('pyHAM solution to B4 BVP')
    textstr = '\n'.join([
        r'$\bf{BVP \/\ & \/\ HAM \/\ parameters}$',
        r'$EI\frac{d^4w}{dx^4} = q, \/\ w(0)=0, \/\ w\'(0)=0, w(1)=0, w\"(1)=0$',
        r'$w_0 = 0., \/\ \mathcal{L} \mathrm{ (w) = \frac{d^4w}{dx^4}}$',
        r'$C_0 = -1$'
    ])
    axes.set_xlabel(r'$x$', fontsize=15)
    axes.set_ylabel(r'$w(x)$', fontsize=15)
    axes.set_xlim(left=0., right=100.)
    # axes.set_ylim(bottom=0.)
    plt.legend(loc='upper left')
    # plt.text(0.5, 0.5*(axes.get_ylim())[1], textstr, fontsize=14, verticalalignment='top', bbox=dict(facecolor='cyan', edgecolor='black'))
    plt.show()


def solve_ode_bvp(ode_subproblem):
    ode_to_use = ode_subproblem  # type: hamTask.HamTaskSubSolution
    force = ode_to_use.task_parent.force
    odeq = ode_to_use.task_parent.eq
    bcs = ode_to_use.bcs
    x0 = ode_to_use.initial_guess
    lin_op = ode_to_use.lin_op
    approx = ode_to_use.approximation
    c0 = ode_to_use.c0
    error_calculator = ode_to_use.task_parent.error_calculator
    error_calculator.deq = ode_to_use.task_parent.eq
    solver = DeformEqODE(odeq, bcs, x0, lin_op, c0)
    if solver.check_bcs():
        if solver.valid_x0_for_bcs() or force:
            ode_to_use.answer.value, sympy_solved = solver.calc_solution(approx)
            if sympy_solved:
                ode_to_use.answer.error = error_calculator.calc(ode_to_use.answer.value, bcs)
                anal_sol = p1*x**2*(3*p3**2 - 5*p3*x + 2*x**2)/(48*p2)
                # plot_to_present_a(solver, error_calculator, anal_sol)
                ode_to_use.answer.status = hamTask.SubSolutionStatus.SOLVED
            else:
                ode_to_use.answer.error = None
                ode_to_use.answer.status = hamTask.SubSolutionStatus.SYMPY_FAILED_FALL_BACK
        else:
            ode_to_use.answer.value = None
            ode_to_use.answer.error = None
            ode_to_use.answer.status = hamTask.SubSolutionStatus.BAD_INITIAL_GUESS
            print('BAD INITIAL GUESS')
    else:
        ode_to_use.answer.error = None
        ode_to_use.answer.value = None
        ode_to_use.answer.status = hamTask.SubSolutionStatus.NOT_SUFFICIENT_BCS


def solve_kernel_ODE(ode_subproblem):
    # print("going to solve subproblem: ")
    # print(ode_subproblem)
    eq = ode_subproblem.task_parent.eq
    bcs = ode_subproblem.bcs
    x0 = ode_subproblem.initial_guess
    ode_subproblem_answer = hamTask.SolutionAnswer()
    ode_subproblem.answer = ode_subproblem_answer
    if check_if_solves_ode_bvp(eq, bcs, x0) and False:
        ode_subproblem_answer.status = hamTask.SubSolutionStatus.SOLVED
        ode_subproblem_answer.value = x0
        ode_subproblem_answer.error = eval_ode_bvp_ans_error(eq, bcs, x0, ode_subproblem.ending_times)
    else:
        solve_ode_bvp(ode_subproblem)


class OdeErrorCalculator:
    def __init__(self):
        self.points = 10
        self.step = None
        self.domain_start = 0
        self.domain_end = 1
        self.deq = None

    def _subs_solution(self, solution):
        ret_expr = S(0)
        lexicon = {}
        for i in range(len(ODE_F_EXPRESSIONS)):
            curr_fder_symbol = ODE_F_EXPRESSIONS[i]
            if curr_fder_symbol in self.deq.free_symbols:
                lexicon[curr_fder_symbol] = solution.diff(x, i)
        # print('deq: ', self.deq)
        # print('lexicon: ', lexicon)
        ret_expr += self.deq.subs(lexicon)
        return ret_expr

    def _domain_def(self, bcs):
        start_pnt = +oo
        end_pnt = -oo
        for curr_bcs in bcs:  # type: hamTask.BoundaryCondition
            if curr_bcs.x < start_pnt:
                start_pnt = curr_bcs.x
            if curr_bcs.x > end_pnt:
                end_pnt = curr_bcs.x
        if start_pnt is not +oo:
            self.domain_start = start_pnt
        if end_pnt is not -oo and end_pnt > self.domain_end:
            self.domain_end = end_pnt
        if self.domain_start == self.domain_end:
            self.domain_end = self.domain_start + 1
            # ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')


    def _calc_with_subs(self, solution, bcs):
        ans = S(0)
        point = self.domain_start
        if self.step:
            step_to_use = self.step
        else:
            step_to_use = (self.domain_end - self.domain_start)/(self.points + 1)
        expr_to_eval = simplify(self._subs_solution(solution))
        while point <= self.domain_end:
            ans += expr_to_eval.evalf(subs={x: point})
            point += step_to_use
        return ans

    def _calc_optimal_C0(self, solution, start, end):
        ret_val = C0
        anal_sol = -1/x
        if C0 in solution.free_symbols:
            integrand = self._subs_solution(solution)
            # integrand = solution - anal_sol
            integrand *= integrand
            # integrand = integrand.subs({C0: -0.11077})
            squared_residual = integrate(integrand, (x, start, end))
            # print("total error is ", squared_residual)
            # assert False
            det_eq = squared_residual.diff(C0)
            print('squared res is ')
            print(det_eq)
            ans = solveset(det_eq, C0, domain=S.Reals)
            print("C0 opt vals are: ", ans)
            ISR_l = lambdify(C0, squared_residual)
            # print('min error is ', ISR_l(-0.0633790771568749))
            # print('min error is ', ISR_l(-0.313647060539033))
            # print('min error is ', ISR_l(-0.177974008726031))
            assert False
            pass
        else:
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        return ret_val

    def calc(self, solution, bcs):
        ans = 0.
        # TODO if C0 from parse -> calc optimal C0 else just produce error
        lambdify_holds = True
        self._domain_def(bcs)
        point = self.domain_start
        # optimal_C0 = -1.
        # sol_to_use = solution.subs({C0: optimal_C0})
        sol_to_use = solution
        # sol_to_use = solution.subs({p1: 1.})
        # self.deq = self.deq.subs({p1: 1.})
        # optimal_C0 = self._calc_optimal_C0(sol_to_use, self.domain_start, self.domain_end)
        # print('optimal C0 is ', optimal_C0)
        # print('sol to use:')
        # print(sol_to_use)
        # return S(0)
        if self.step:
            step_to_use = self.step
        else:
            step_to_use = (self.domain_end - self.domain_start)/(self.points + 1)
        # print('sol to use: ', sol_to_use)
        expr_to_eval = simplify(self._subs_solution(sol_to_use))
        # print('expr eval: ', expr_to_eval)
        # anal_sol = p1*x**2*(3*p3**2 - 5*p3*x + 2*x**2)/(48*p2)
        # anal_sol = anal_sol.subs({p3: 1.})
        # expr_to_eval = simplify(sol_to_use - anal_sol)
        expr_eval_l = lambdify(x, expr_to_eval, modules=['math', 'sympy'])
        # print('expr to eval')
        while point <= self.domain_end and lambdify_holds:
            try:
                curr_error = expr_eval_l(point)
                ans += curr_error*curr_error
                point += step_to_use
            except TypeError:
                print('lambdify failed')
                assert False
                lambdify_holds = False
        if lambdify_holds:
            ans = simplify(ans)
        else:
            ans = self._calc_with_subs(sol_to_use, bcs)
        if x in ans.free_symbols:
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        # print('error ')
        # print(ans)
        # ans /= self.points
        # opt = ans.diff(C0)
        # opt_ans = solveset(opt, C0, domain=S.Reals)
        # print('opt C0 is ', opt_ans)
        # l_opt = lambdify(C0, ans)
        # print('opt error is ', l_opt(-0.5))
        # print('opt error2 is ', l_opt(-0.319094219626965))
        # print('opt error3 is ', l_opt(-0.15812531498891))
        return ans

