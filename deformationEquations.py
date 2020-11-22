from sympy import simplify, solve, sympify, dsolve, Eq, linsolve, nonlinsolve
from sympy import FiniteSet, EmptySet, exp, Derivative, Integral, pdsolve

import hamTask
from hamErrors import ham_error
from hamSymbols import *


def zero_val_float_list(num):
    a_list = []
    for i in range(num):
        a_list.append(0.)
    return a_list


def construct_zero_deform_eq(c0):
    deform_eq = (1 - q)*(f(h(q)) - fh0) - c0*q*f(h(q))
    return deform_eq


def try_clear_expr(expr):
    ret_expr = expr
    uneval_integrals = expr.atoms(Integral)
    if uneval_integrals and len(uneval_integrals) > 0:
        times = 0
        unevals = len(uneval_integrals)
        while unevals > 0 and times < 10:
            ret_expr = ret_expr.doit()
            uneval_integrals = ret_expr.atoms(Integral)
            unevals = len(uneval_integrals)
            times += 1
        if times >= 10:
            ret_expr = expr
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
    return ret_expr


class DeformEqAlgebr:
    def __init__(self, c0, expr_in, x0):
        self.c0 = c0
        self.x0 = x0
        self.expr = expr_in
        self.homotopy_terms = []
        self.homotopy_terms.append(x0)
        self.deform_eqs = []
        self.deform_eqs.append(construct_zero_deform_eq(c0))
        self.f_vals = []
        self.f_vals.append(self._eval_fun(self.x0))

    def _eval_fun(self, val):
        return self.expr.evalf(subs={x: val})

    def _build_f_vals(self, n):
        prev_len = len(self.f_vals)
        times = n - len(self.f_vals) + 1
        for i in range(times):
            der_num = i + prev_len
            ret_val = (self.expr.diff(x, der_num)).evalf(subs={x: self.x0})
            self.f_vals.append(ret_val)

    def _build_temp_lexicon(self, n):
        ret_dict = {}
        for i in range(n + 1):
            ret_dict[FDER_EVAL_LIST[i]] = self.f_vals[i]
        return ret_dict



    def _calc_terms(self, n):
        prev_num = len(self.deform_eqs)
        if prev_num != len(self.homotopy_terms):
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        if n > prev_num - 1:
            self._build_f_vals(n)
            temp_lex = self._build_temp_lexicon(n)
            d_times = n - prev_num + 1
            for i in range(d_times):
                self._populate(temp_lex)

    def _enrich_dict_with_x_vals(self, dict_to_enrich):
        for i in range(len(self.homotopy_terms)):
            dict_to_enrich[XTERMS_LIST[i]] = self.homotopy_terms[i]

    def _populate(self, fvals_lexicon):
        last_deform_eq = self.deform_eqs[-1]
        print('last deform eq is ', last_deform_eq)
        ldeq_diffed = simplify(last_deform_eq.diff(q))
        print('diffed: ', ldeq_diffed)
        self.deform_eqs.append(ldeq_diffed)
        ldeq_subs_homotopies = ldeq_diffed.xreplace(ALG_LEXICON)
        print('homotopies after xreplace: ', ldeq_subs_homotopies)
        ldeq_subs_homotopies = ldeq_subs_homotopies.xreplace(ALG_LEXICON_SECOND)
        print("after secondary xrepl: ", ldeq_subs_homotopies)
        ldeq_subs_homotopies = ldeq_subs_homotopies.xreplace(ALG_LEXICON_THIRD)
        print("after third xrepl: ", ldeq_subs_homotopies)
        var_to_solve_for = XTERMS_LIST[len(self.deform_eqs) - 1]
        print('var to solve for: ', var_to_solve_for)
        next_term = solve(ldeq_subs_homotopies, var_to_solve_for)
        print("next term is ", next_term)
        temp_lex = fvals_lexicon.copy()
        self._enrich_dict_with_x_vals(temp_lex)
        if next_term and len(next_term) == 1:
            next_term_calced = next_term[0].subs(temp_lex)
        else:
            ham_error('HAM_ERROR_HANDLED', 'HANDLED_GENERAL')
            next_term_calced = sympify(0.)
        self.homotopy_terms.append(next_term_calced)

    def calc_solution(self, n):
        ans = None
        if n >= 0:
            self._calc_terms(n)
            print("Homotopy terms analytically: ", self.homotopy_terms)
            # TODO here return the answer only by adding n terms and not all stored
            ans = sum(self.homotopy_terms)
        else:
            ham_error('HAM_ERROR_HANDLED', 'HANDLED_REQUEST_FOR_NEGATIVE_APPROXIMATIONS_NUMBER')
        return ans


class DeformEqODE:
    def __init__(self, odeq, bcs, x0, lin_op, c0):
        self.eq = odeq
        self.bcs = bcs
        self.x0 = x0
        self.lin_op = lin_op
        self.c0 = c0
        self.homotopy_terms = [x0]
        self.deform_eqs = []
        self.ode_order = self._ode_order()
        self.sympy_solved = True
        self.deqs_to_solve_export = []

    @staticmethod
    def _n_op_replace_homotopy_fun(n_op_expr, ode_fder_list):
        lexicon = {}
        for i in range(len(ode_fder_list)):
            curr_fder_symbol = ode_fder_list[i]
            lexicon[curr_fder_symbol] = f(x, q).diff(x, i)
        ret_expr = n_op_expr.subs(lexicon)
        return ret_expr

    def _calc_replacement(self, der):
        d_metrics = der.variable_count
        x_times = 0
        q_times = 0
        for var_count in d_metrics:
            if var_count[0] == x:
                x_times = var_count[1]
            elif var_count[0] == q:
                q_times = var_count[1]
        replacer = factorial(q_times)*self.homotopy_terms[q_times]
        replacer = replacer.diff(x, x_times)
        if q in replacer.atoms(Symbol):
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        return replacer

    def _eval_at_zero_q(self, rhs):
        der_list = rhs.atoms(Derivative)
        repl_lexicon = {}
        for der in der_list:
            repl_der = self._calc_replacement(der)
            repl_lexicon[der] = repl_der
        rhs = rhs.xreplace(repl_lexicon)
        rhs = simplify(rhs.xreplace({f(x, q): self.homotopy_terms[0]}))
        return rhs

    def _l_op_replace_fun(self):
        lexicon = {}
        for i in range(len(ODE_F_EXPRESSIONS)):
            curr_fder_symbol = ODE_F_EXPRESSIONS[i]
            if curr_fder_symbol in self.lin_op.free_symbols:
                lexicon[curr_fder_symbol] = f(x).diff(x, i)
        return self.lin_op.subs(lexicon)

    @staticmethod
    def _check_if_solution_explicit(deq, sol):
        return True
        # TODO find a way to check it. possibly differentiate

    def _parse_solution(self, deq):
        sol = None
        state = True
        try:
            ans = dsolve(deq, f(x), real=True)
            if isinstance(ans, Eq):
                sol = ans.rhs
                sol = try_clear_expr(sol)
            elif isinstance(ans, list):
                # TODO treat later the multiple solution case. It should check
                # TODO which solution is explicit and then choose somehow among them
                ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
            else:
                ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        except NotImplementedError:
            state = False
            sol = None
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        if state:
            state = self._check_if_solution_explicit(deq, sol)
        return sol, state

    def _produce_integr_constants_list(self):
        a_list = []
        for i in range(self.ode_order):
            a_list.append(INTEGRATION_CONSTANTS[i])
        return a_list

    def _calc_integration_constants(self, solution):
        ret_expr = solution
        lin_eqs = []
        constants = self._produce_integr_constants_list()
        for curr_bcs in self.bcs:  # type: hamTask.BoundaryCondition
            a_new_lin_eq = solution.diff(x, curr_bcs.fx_order)
            a_new_lin_eq = a_new_lin_eq.subs({x: curr_bcs.x})
            lin_eqs.append(a_new_lin_eq)
        assert len(constants) == len(lin_eqs)
        try:
            ans = linsolve(lin_eqs, constants)
        except ValueError:
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
            ans = nonlinsolve(lin_eqs, constants)
        if isinstance(ans, FiniteSet):
            vals_lex = {}
            ans = (list(ans))[0]
            for i in range(len(ans)):
                vals_lex[constants[i]] = ans[i]
            ret_expr = solution.subs(vals_lex)
        elif isinstance(ans, EmptySet):
            self.sympy_solved = False
        return ret_expr

    def _populate(self):
        m = len(self.homotopy_terms)
        c0 = self.c0
        if m > 1:
            last_term = self.homotopy_terms[m - 1]
        else:
            last_term = S(0)
        rhs = self._n_op_replace_homotopy_fun(self.eq, ODE_F_EXPRESSIONS)
        rhs = rhs.diff(q, m - 1)
        rhs = self._eval_at_zero_q(rhs)
        # ct = Symbol('ct', real=True)
        ct = c0
        rhs *= ct
        lhs = self._l_op_replace_fun()
        self.deqs_to_solve_export.append(lhs-rhs)
        if self.sympy_solved:
            solution, self.sympy_solved = self._parse_solution(lhs-rhs)
            if self.sympy_solved:
                solution = simplify(solution + last_term)
                solution = simplify(self._calc_integration_constants(solution))
                self.homotopy_terms.append(solution)
                # print('solution is ', solution)

    def _calc_terms(self, n):
        prev_num = len(self.homotopy_terms)
        if n > prev_num - 1:
            d_times = n - prev_num + 1
            for i in range(d_times):
                self._populate()

    def _ode_order(self):
        ans = 0
        symbols_in_ode = self.eq.free_symbols
        gen = (var for var in symbols_in_ode if var in ODE_F_EXPRESSIONS)  # generator expression
        for curr_symbol in gen:
            index = ODE_F_EXPRESSIONS.index(curr_symbol)
            if index > ans:
                ans = index
        return ans

    def _check_bcs_number(self):
        return len(self.bcs) == self.ode_order

    def _check_max_derivative(self):
        ans = True
        for curr_bcs in self.bcs:
            if curr_bcs.fx_order >= self.ode_order:
                ans = False
                break
        return ans

    @staticmethod
    def _non_coincident_bcs(bcs_a, bcs_b):
        return not (bcs_a.x == bcs_b.x and bcs_a.fx_order == bcs_b.fx_order)

    def _bcs_coincident(self, ref_bcs):
        ans = True
        gen = (bcs for bcs in self.bcs if bcs != ref_bcs)
        for another_bcs in gen:
            if not self._non_coincident_bcs(ref_bcs, another_bcs):
                ans = False
                break
        return ans

    def _check_non_coincidence(self):
        ans = True
        if len(self.bcs) > 1:
            for curr_bcs in self.bcs:
                if not self._bcs_coincident(curr_bcs):
                    ans = False
                    break
        return ans

    def check_bcs(self):
        # TODO remember to allow expressions as BCS in the future (e.g. Cauchy conditions)
        return self._check_bcs_number() and self._check_max_derivative() and self._check_non_coincidence()

    def valid_x0_for_bcs(self):
        ans = True
        for curr_bcs in self.bcs:
            if not curr_bcs.satisfied_by_expr(self.x0):
                ans = False
                break
        return ans

    def calc_solution(self, n):
        # TODO should return tuple (value, bool=sympy_solved)
        ans = None
        if n >= 0:
            self._calc_terms(n)
            ans = sum(self.homotopy_terms)
        else:
            ham_error('HAM_ERROR_HANDLED', 'HANDLED__GENERAL')
        return ans, True


class DeformEqPDE:
    def __init__(self, pdeq, bcs, x0, lin_op, c0):
        self.eq = pdeq
        self.bcs = bcs
        self.x0 = x0
        self.lin_op = lin_op
        self.c0 = c0
        self.homotopy_terms = [x0]
        self.deform_eqs = []
        self.pde_orders = self._pde_orders()
        self.sympy_solved = True
        self.deqs_to_solve_export = []

    def _pde_orders(self):
        x_ord = 0
        y_ord = 0
        t_ord = 0
        symbols_in_pde = self.eq.free_symbols
        gen = (var for var in symbols_in_pde if var in PDE_U_EXPRESSIONS)
        for curr_symbol in gen:
            ord_tup = decode_der_orders_pde(curr_symbol)
            if ord_tup[0] > x_ord:
                x_ord = ord_tup[0]
            if ord_tup[1] > y_ord:
                y_ord = ord_tup[1]
            if ord_tup[2] > t_ord:
                t_ord = ord_tup[2]

    def _check_if_sufficient_bcs(self):
        ans = False
        lin_op_orders = orders_of_pde(self.lin_op)
        the_sum = 0
        if lin_op_orders[0] > 0:
            the_sum += lin_op_orders[0]
        if lin_op_orders[1] > 0:
            the_sum += lin_op_orders[1]
        if lin_op_orders[2] > 0:
            the_sum += lin_op_orders[2]
        if the_sum > 0:
            ans = the_sum == len(self.bcs)
        else:
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        return ans

    def _check_for_right_order_bcs(self):
        ans = True
        pdeq_orders = orders_of_pde(self.eq)
        print('pde orders: ', pdeq_orders)
        for curr_bcs in self.bcs:
            bcs_orders = orders_of_pde(curr_bcs.value)
            print('bcs orders: ', bcs_orders)
            if pdeq_orders[0] and bcs_orders[0] >= pdeq_orders[0]:
                ans = False
                break
            if pdeq_orders[1] and bcs_orders[1] >= pdeq_orders[1]:
                ans = False
                break
            if pdeq_orders[2] and bcs_orders[2] >= pdeq_orders[2]:
                ans = False
                break
        return ans

    def check_bcs(self):
        return self._check_if_sufficient_bcs() and self._check_for_right_order_bcs()

    def valid_x0_for_bcs(self):
        ans = True
        for curr_bcs in self.bcs: # type: hamTask.BCPDE
            ans = curr_bcs.satisfied_by(self.x0)
            if not ans:
                break
        return ans

    @staticmethod
    def _n_op_replace_homotopy_fun(n_op_expr):
        lexicon = {}
        the_symbols = n_op_expr.atoms(Symbol)
        gen = (var for var in the_symbols if var in PDE_U_EXPRESSIONS)
        for u_fun in gen:
            x_ord, y_ord, t_ord = decode_der_orders_pde(u_fun)
            u_expr = u(x, y, t, q).diff(x, x_ord, y, y_ord, t, t_ord)
            lexicon[u_fun] = u_expr
        ret_expr = n_op_expr.subs(lexicon)
        return ret_expr

    def _calc_replacement(self, der):
        d_metrics = der.variable_count
        x_times = 0
        y_times = 0
        t_times = 0
        q_times = 0
        for var_count in d_metrics:
            if var_count[0] == x:
                x_times = var_count[1]
            elif var_count[0] == y:
                y_times = var_count[1]
            elif var_count[0] == t:
                t_times = var_count[1]
            elif var_count[0] == q:
                q_times = var_count[1]
            else:
                ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        replacer = factorial(q_times)*self.homotopy_terms[q_times]
        replacer = replacer.diff(x, x_times, y, y_times, t, t_times)
        if q in replacer.atoms(Symbol):
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        return replacer

    def _eval_at_zero_q(self, rhs):
        der_list = rhs.atoms(Derivative)
        repl_lexicon = {}
        for der in der_list:
            repl_der = self._calc_replacement(der)
            repl_lexicon[der] = repl_der
        rhs = rhs.xreplace(repl_lexicon)
        rhs = simplify(rhs.xreplace({u(x, y, t, q): self.homotopy_terms[0]}))
        return rhs

    def _l_op_replace_fun(self):
        lexicon = {}
        the_symbols = self.lin_op.atoms(Symbol)
        gen = (var for var in the_symbols if var in PDE_U_EXPRESSIONS)
        for u_fun in gen:
            x_ord, y_ord, t_ord = decode_der_orders_pde(u_fun)
            lexicon[u_fun] = u(x, y, t).diff(x, x_ord, y, y_ord, t, t_ord)
        return self.lin_op.subs(lexicon)

    def _check_if_solution_explicit(self, deq, sol):
        # TODO find a way to check if solution is explicit
        return True

    def _parse_solution(self, deq):
        sol = None
        status = True
        try:
            ans = pdsolve(deq, u(x, y, t), real=True)
            if isinstance(ans, Eq):
                sol = ans.rhs
                sol = try_clear_expr(sol)
                # TODO handle later the multiple solution case
            elif isinstance(ans, list):
                ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__NOT_IMPLEMENTED_YET')
            else:
                ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        except NotImplementedError:
            status = False
            sol = None
            # TODO temp treatment. FALL BACK to EQ_GEN_MODE
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        if status:
            status = self._check_if_solution_explicit(deq, sol)
        return sol, status

    def _populate(self):
        m = len(self.homotopy_terms)
        c0 = self.c0
        if m > 1:
            last_term = self.homotopy_terms[m - 1]
        else:
            last_term = S(0)
        rhs = self._n_op_replace_homotopy_fun(self.eq)
        pass

    def _calc_terms(self, n):
        prev_num = len(self.homotopy_terms)
        if n > prev_num - 1:
            d_times = n - prev_num + 1
            for i in range(d_times):
                self._populate()

    def calc_solution(self, n):
        # TODO should return tuple (value, sympy_solved)
        ans = None
        if n >= 0:
            self._calc_terms(n)
            print("Homotopy terms analytically: ", self.homotopy_terms)
            # TODO return only n terms and not all sum
            ans = sum(self.homotopy_terms)
        else:
            ham_error('HAM_ERROR_HANDLED', "HANDLED_REQUEST_FOR_NEGATIVE_APPROXIMATIONS_NUMBER")
        return ans, True

