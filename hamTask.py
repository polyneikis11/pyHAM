import os
import csv
from enum import Enum
from sympy import S, diff, SympifyError, sympify, simplify, LambertW, exp, solveset
from sympy import FiniteSet
from hamErrors import ham_error
from hamBase import HamEnt, ham_any
from algebraicEqSolver import solve_kernel_algebraic
from odeSolver import solve_kernel_ODE, OdeErrorCalculator
from pdeSolver import solve_kernel_PDE, PdeErrorCalculator
from StefanSolver import solve_kernel_Stefan
from hamSymbols import *
from geometry import HamPoint, str_to_float

# TODO next: print linear operators in __str__ of ODE task, BC class (base + derived) and prints
# GLOBAL CONSTANTS

MAX_APPROX_ALGEBRAIC = 10


# *******************************************************************************
# GLOBAL FUNCTIONS

def is_csv_file(filename):
    ans = False
    if filename and os.path.exists(filename):
        res = os.path.splitext(filename)
        ans = res[1] == '.csv'
    return ans




def check_eq_algebr_valid(s_expr):
    return s_expr.diff(x) != S(0)


def check_eq_ode_valid(s_expr):
    ans = False
    for curr_f_der in range(len(ODE_F_EXPRESSIONS)):
        if 0 != curr_f_der:
            der = s_expr.diff(ODE_F_EXPRESSIONS[curr_f_der])
            if der != S(0):
                ans = True
                break
    return ans


def check_eq_pde_valid(s_expr):
    ans = False
    symbols_in_expr = s_expr.atoms(Symbol)
    for curr_symbol in symbols_in_expr:
        if curr_symbol != PDE_U_EXPRESSIONS[0] and curr_symbol in PDE_U_EXPRESSIONS:
            ans = True
            break
    return ans


def check_ode_initial_guess_valid(s_expr):
    ans = True
    for curr_f_der in ODE_F_EXPRESSIONS:
        der = s_expr.diff(curr_f_der)
        if der != S(0):
            ans = False
            break
    return ans


def try_to_sympify_and_ret_val(str_expr):
    try:
        ret_val = ham_sympify(str_expr)
    except SympifyError:
        ret_val = None
    return ret_val


def check_lin_op_valid(expr):
    ans = True
    # TODO test if linear operator is a) linear, b) L(0)=0. If not,
    # TODO c, abort specific subproblem and refer it by appropriate error enumeration.
    return ans


def parse_possible_expr_str_iterable(obj, store_here):
    if type(obj) == list:
        for a_val in obj:
            ret_val = try_to_sympify_and_ret_val(a_val)
            if ret_val:
                store_here.append(ret_val)
    else:
        ret_val = try_to_sympify_and_ret_val(obj)
        if ret_val is not None:
            store_here.append(ret_val)

# ***************************************************************
# ENUMERATIONS


class TaskType(Enum):
    ALGEBRAIC_EQ_TASK = 0
    ODE_TASK = 1
    PDE_TASK = 2
    STEFAN_TASK = 3

# ----------------------------------------------------------------`


class SubSolutionStatus(Enum):
    SOLVED = 0
    NOT_SET_TO_SOLVE = 1
    INVALID_EQUATION = 2
    BAD_INITIAL_GUESS = 3
    SERIES_NOT_CONVERGED = 4
    SYMPY_FAILED_FALL_BACK = 5  # just produce the list of BVP's to be solved
    NOT_SUFFICIENT_BCS = 6
    NOT_INIT = 7


# ***************************************************************
# CLASSES
class BoundaryCondition(HamEnt):
    def __init__(self, bvp_type):
        super().__init__('HAM_BC')
        self.x = None
        self.y = None
        self.z = None
        self.t = None
        self.fx_order = None
        self.fy_order = None
        self.fz_order = None
        self.ft_order = None
        self.value = None
        self.bvp_type = bvp_type

    def _check_ODE_bc(self):
        assert free_of_symbols_except(self.x, GLOBAL_PARAMETERS)
        self.y = None
        self.z = None
        self.t = None
        assert type(self.fx_order) is int and self.fx_order >= 0 and self.fx_order < 10
        self.fy_order = None
        self.fz_order = None
        self.ft_order = None

    def _check_PDE_bc(self):
        pass
        #TODO

    def check_type(self):
        if TaskType.ALGEBRAIC_EQ_TASK == self.bvp_type:
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        elif TaskType.ODE_TASK == self.bvp_type:
            self._check_ODE_bc()
        elif TaskType.PDE_TASK == self.bvp_type or TaskType.STEFAN_TASK == self.bvp_type:
            self._check_PDE_bc()
        return

    def set_f_der_val(self, f_symbol, f_val):
        index = ODE_F_EXPRESSIONS.index(f_symbol)
        if index > -1:
            self.fx_order = index
            self.value = f_val

    def satisfied_by_expr(self, expr):
        diff_x_times = 0
        diff_y_times = 0
        diff_z_times = 0
        diff_t_times = 0
        if self.fx_order:
            diff_x_times = self.fx_order
        if self.fy_order:
            diff_y_times = self.fy_order
        if self.fz_order:
            diff_z_times = self.fz_order
        if self.ft_order:
            diff_t_times = self.ft_order
        use_expr = expr.diff(x, diff_x_times, y, diff_y_times, z, diff_z_times, t, diff_t_times)
        temp_lex = {x: self.x, y: self.y, z: self.z, t: self.t}
        use_expr = simplify(use_expr.subs(temp_lex))
        return simplify(self.value) == simplify(use_expr)

    def _symbol_str(self):
        ret_val = None
        if self.bvp_type == TaskType.ODE_TASK:
            ret_val = str(ODE_F_EXPRESSIONS[self.fx_order])
        elif self.bvp_type == TaskType.PDE_TASK:
            symb_list = ['u', 'x', str(self.fx_order), 'y', str(self.fy_order), 't', str(self.ft_order)]
            ret_val = ''.join(symb_list)
        else:
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        return ret_val

    def __str__(self):
        ret_str = ''
        if self.bvp_type == TaskType.ODE_TASK:
            ret_str = '( x: ' + str(self.x) + ', '
            ret_str += self._symbol_str() + ': ' + str(self.value) + ')\n'
        elif self.bvp_type == TaskType.PDE_TASK:
            ret_str += '( '
            if self.x is not None:
                ret_str += 'x: ' + str(self.x) + ', '
            if self.y is not None:
                ret_str += 'y: ' + str(self.y) + ', '
            if self.t is not None:
                ret_str += 't: ' + str(self.t) + ', '
            ret_str += self._symbol_str() + ': ' + str(self.value) + ')\n'
        else:
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        return ret_str


# -------------------------------------------------------------------------
class BCPDE:
    def __init__(self):
        self.boundary = None
        self.value = None
        self.explicit = False
        self.explicit_symbol = None

    def __str__(self):
        ret_str = '----- PDE Boundary Condition -----\n'
        if self.explicit_symbol is not None:
            ret_str += '----- Boundary: {} = {}'.format(self.explicit_symbol, self.boundary) + '\n'
        elif self.boundary is not None:
            ret_str += '----- Boundary: {}'.format(self.boundary) + '\n'
        else:
            ret_str += '----- Boundary: Not defined\n'
        if self.value is not None:
            ret_str += '----- Value: {}'.format(self.value) + '\n'
        else:
            ret_str += '----- Value: Not defined\n'
        ret_str += '===================================='
        return ret_str

    def replace_params_vals(self, params_lex):
        self.boundary = self.boundary.subs(params_lex)
        self.value = self.value.subs(params_lex)

    @staticmethod
    def _check_boundary(boundary_in):
        ans = False
        symbols_list = boundary_in.atoms(Symbol)
        if x in symbols_list or y in symbols_list or t in symbols_list:
            if not contains_ufun_symbols(boundary_in):
                ans = True
        return ans

    @staticmethod
    def _check_value(value_in):
        return contains_ufun_symbols(value_in)

    def _try_explicit_boundary(self):
        ans = None
        the_symbols = self.boundary.atoms(Symbol)
        solved = False
        if x in the_symbols:
            try:
                ans = solveset(self.boundary, x, domain=S.Reals)
                if isinstance(ans, FiniteSet):
                    solved = True
                    self.explicit_symbol = x
            except:
                pass
        elif y in the_symbols and not solved:
            try:
                ans = solveset(self.boundary, y, domain=S.Reals)
                if isinstance(ans, FiniteSet):
                    solved = True
                    self.explicit_symbol = y
            except:
                pass
        elif t in the_symbols and not solved:
            try:
                ans = solveset(self.boundary, t, domain=S.Reals)
                if isinstance(ans, FiniteSet):
                    solved = True
                    self.explicit_symbol = t
            except:
                pass
        else:
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        if solved:
            ans = list(ans)
            if len(ans) > 1:
                ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
            self.boundary = simplify(ans[0])
            self.explicit = True

    def set_boundary(self, boundary_in):
        if self._check_boundary(boundary_in):
            self.boundary = boundary_in
            self._try_explicit_boundary()

    def set_value(self, value_in):
        if self._check_value(value_in):
            self.value = value_in

    def is_well_defined(self):
        return self.boundary is not None and self.value is not None

    def max_orders(self):
        x_ord = -1
        y_ord = -1
        t_ord = -1
        symb_list = self.value.atoms(Symbol)
        gen = (var for var in symb_list if var in PDE_U_EXPRESSIONS)
        for the_symb in gen:
            temp_x_ord, temp_y_ord, temp_t_ord = decode_der_orders_pde(the_symb)
            if temp_x_ord > x_ord:
                x_ord = temp_x_ord
            if temp_y_ord > y_ord:
                y_ord = temp_y_ord
            if temp_t_ord > t_ord:
                t_ord = temp_t_ord
        return x_ord, y_ord, t_ord

    def satisfied_by(self, expr):
        ans = False
        if self.explicit:
            subs_lexicon = {}
            the_symbols = self.value.atoms(Symbol)
            gen = (var for var in the_symbols if var in PDE_U_EXPRESSIONS)
            for u_fun in gen:
                x_ord, y_ord, t_ord = decode_der_orders_pde(u_fun)
                expr_to_use = expr.diff(x, x_ord)
                expr_to_use = expr_to_use.diff(y, y_ord)
                expr_to_use = expr_to_use.diff(t, t_ord)
                subs_lexicon[u_fun] = simplify(expr_to_use)
            gen = (var for var in the_symbols if var in MOVING_BOUNDARY_EXP)
            for r_fun in gen:
                indx = MOVING_BOUNDARY_EXP.index(r_fun)
                subs_lexicon[r_fun] = expr.diff(t, indx)
            subs_lexicon[self.explicit_symbol] = self.boundary
            outcome = simplify(self.value.subs(subs_lexicon))
            print('outcome is ', outcome)
            ans = outcome == S(0) or float(outcome) < 1e-10
        else:
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        return ans
# --------------------------------------------------------------------


class HamTaskSubSolution:
    def __init__(self, task_type):
        self.task_parent = None
        self.c0 = None
        self.initial_guess = None
        self.approximation = None
        self.answer = None
        self.task_type = task_type
        self.lin_op = None
        self.bcs = None
        self.force = False

    def _additional_str_ode(self):
        ret_str = ''
        ret_str += self.task_parent.boundary_conditions_str()
        ret_str += 'L operator: ' + str(self.lin_op) + '\n'
        return ret_str

    def __str__(self):
        a_val = ''
        a_val += str(self.task_type) + '\n'
        a_val += 'Eq: ' + str(self.task_parent.eq) + '\n'
        if self.task_type == TaskType.ODE_TASK:
            a_val += self._additional_str_ode()
        elif self.task_type == TaskType.PDE_TASK:
            a_val += self._additional_str_ode()
        if self.c0:
            a_val += 'C0: ' + str(self.c0) + '\n'
        else:
            a_val += 'C0: N/A\n'
        if self.initial_guess:
            a_val += 'initial_guess: ' + str(self.initial_guess) + '\n'
        else:
            a_val += 'initial_guess: N/A\n'
        if self.approximation:
            a_val += 'approximation order: ' + str(self.approximation) + '\n'
        else:
            a_val += 'approximation order: N/A\n'
        if self.answer:
            a_val += 'answer: ' + str(self.answer) + '\n'
        else:
            a_val += 'answer: N/A\n'
        a_val += '-'*10 + '\n'
        return a_val

    def solve_kernel(self):
        if TaskType.ALGEBRAIC_EQ_TASK == self.task_type:
            solve_kernel_algebraic(self)
        elif TaskType.ODE_TASK == self.task_type:
            solve_kernel_ODE(self)
        elif TaskType.PDE_TASK == self.task_type:
            solve_kernel_PDE(self)
        elif TaskType.STEFAN_TASK == self.task_type:
            solve_kernel_Stefan(self)
        else:
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')

# -------------------------------------------------------------------------


class SolutionAnswer:
    def __init__(self):
        self.status = SubSolutionStatus.NOT_INIT
        self.value = None
        self.error = None

    def __str__(self):
        return 'value: {}, error: {}'.format(self.value, self.error)

# -------------------------------------------------------------------------


class HamTask(HamEnt):
    def __init__(self, task_type):
        if task_type == 'HAM_TASK':
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__VIRTUAL_METHOD_CALLED')
        else:
            super().__init__(task_type)
            self.eq = None
            self.parsing_status = False
            self.initial_guesses = []
            self.c0_vals = []
            self.approximations = []
            self.to_solve = True
            self.solutions = []

    def set_parsing_status(self, status):
        self.parsing_status = status

    @staticmethod
    def auto_calc_initial_guess():
        ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__VIRTUAL_METHOD_CALLED')

    @staticmethod
    def solve():
        ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__VIRTUAL_METHOD_CALLED')

    @staticmethod
    def build(self, inp_dict):
        ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__VIRTUAL_METHOD_CALLED')

    def __str__(self):
        ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__VIRTUAL_METHOD_CALLED')
        return "a_blank_line"

# -------------------------------------------------------------------------


class HamTaskAlgebraic(HamTask):
    def __init__(self):
        super().__init__('HAM_ALG_TASK')
        self.num_iter = 0

    def set_eq(self, expr):
        # TODO add an exception here to catch SympifyError in order to handle it.
        inc_expr = ham_sympify(expr)
        if not check_eq_algebr_valid(inc_expr):
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__BAD_ALGEBRAIC_EQUATION')
        else:
            self.eq = inc_expr

    def set_initial_guess(self, initial_guess_expr):
        # cast to float.
        # TODO if eq present. check that the derivative is not zero
        self.initial_guess = float(initial_guess_expr)

    def set_approximation_rank(self, approx_rank):
        # check if it is a positive integer and that is below a threshold
        if type(approx_rank) == int:
            if approx_rank < 1:
                ham_error('HAM_ERROR_HANDLED', 'HANDLED__INVALID_APPROX_RANK_GIVEN_NON_POSITIVE_INTEGER')
            elif approx_rank > MAX_APPROX_ALGEBRAIC:
                ham_error('HAM_ERROR_HANDLED', 'HANDLED__INVALID_APPROX_RANK_GIVEN_BIGGER_THAN_MAX')
            else:
                self.approximation = approx_rank
        else:
            ham_error('HAM_ERROR_HANDLED', 'HANDLED__INVALID_APPROX_RANK_GIVEN_NON_POSITIVE_INTEGER')

    @staticmethod
    def set_num_iter(self, num_iter):
        # not implemented yet error
        ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__NOT_IMPLEMENTED_YET')

    def auto_calc_initial_guess(self):
        # TODO for the time just return zero. later we will think of something
        self.initial_guess = sympify("0.")

    def solve(self):
        if self.parsing_status and self.to_solve:
            for initial_guess in self.initial_guesses:
                for c0_curr_val in self.c0_vals:
                    for approx_rank in self.approximations:
                        subproblem = HamTaskSubSolution(TaskType.ALGEBRAIC_EQ_TASK)
                        subproblem.c0 = c0_curr_val
                        subproblem.approximation = approx_rank
                        subproblem.initial_guess = initial_guess
                        subproblem.task_parent = self
                        subproblem.solve_kernel()
                        self.solutions.append(subproblem)
        else:
            print("Task {:d} is not marked for solution".format(self.id))

    def build_optionals(self, inp_dict):
        t_solve = inp_dict.get('solve')
        if t_solve:
            if 'n' == t_solve:
                self.to_solve = False
        t_init_guess = inp_dict.get('initial_guess')
        if t_init_guess:
            parse_possible_expr_str_iterable(t_init_guess, self.initial_guesses)
        else:
            self.auto_calc_initial_guess()
        t_c0 = inp_dict.get('c0')
        if t_c0:
            parse_possible_expr_str_iterable(t_c0, self.c0_vals)
        else:
            self.c0_vals.append(sympify("-1"))
        t_approx = inp_dict.get('approximation')
        if t_approx:
            parse_possible_expr_str_iterable(t_approx, self.approximations)
        else:
            self.approximations.append(2)

    def build(self, inp_dict):
        t_eq = inp_dict.get('equation')
        if t_eq:
            try:
                self.eq = ham_sympify(t_eq)
                self.parsing_status = True
                self.build_optionals(inp_dict)
            except SympifyError:
                self.parsing_status = False
        else:
            self.parsing_status = False

    def __str__(self):
        ret_str = ''
        ret_str += 'OBJECT: HAM task algebraic\n'
        ret_str += 'ID: ' + str(self.id) + '\n'
        if self.parsing_status:
            ret_str += 'Equation: ' + str(self.eq) + '\n'
            ret_str += 'Initial Guess: {'
            for cnt in range(len(self.initial_guesses)):
                if cnt != (len(self.initial_guesses) - 1):
                    ret_str += str(self.initial_guesses[cnt]) + ', '
                else:
                    ret_str += str(self.initial_guesses[cnt]) + '}\n'
            ret_str += 'C0: {'
            for cnt in range(len(self.c0_vals)):
                if cnt != (len(self.c0_vals) - 1):
                    ret_str += str(self.c0_vals[cnt]) + ', '
                else:
                    ret_str += str(self.c0_vals[cnt]) + '}\n'
            ret_str += 'Approximations: {'
            for cnt in range(len(self.approximations)):
                if cnt != (len(self.approximations) - 1):
                    ret_str += str(self.approximations[cnt]) + ', '
                else:
                    ret_str += str(self.approximations[cnt]) + '}\n'
            ret_str += '------ SOLUTIONS ({}) ------\n'.format(len(self.solutions))
            for a_solution in self.solutions:
                ret_str += a_solution.__str__()
            ret_str += '-----------------------\n'
        else:
            ret_str += '*** Failed to parse data. Empty object ***'
            # TODO possibly here refer the lines in the input file
        ret_str += '='*20
        return ret_str

# -------------------------------------------------------------------------


class HamTaskODE(HamTask):
    def __init__(self):
        super().__init__('HAM_ODE_TASK')
        self.bcs = []  # boundary conditions
        self.lin_operator = []  # linear operator
        self.error_calculator = OdeErrorCalculator()
        self.force = False

    def set_eq(self, expr):
        inc_expr = ham_sympify(expr)
        if not check_eq_ode_valid(inc_expr):
            ham_error('HAM_ERROR_CRITICAL', "CRITICAL__ASSERT_NOT_REACHED")
        else:
            self.eq = inc_expr

    def auto_calc_initial_guess(self):
        ham_error("HAM_ERROR_CRITICAL", "CRITICAL__NOT_IMPLEMENTED_YET")

    def auto_calc_lin_op(self):
        ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__NOT_IMPLEMENTED_YET')

    def solve(self):
        # print('initial guesses: ', self.initial_guesses)
        # print("co vals: ", self.c0_vals)
        # print("approx: ", self.approximations)
        # print("lin op: ", self.lin_operator)
        if self.parsing_status and self.to_solve:
            for initial_guess in self.initial_guesses:
                for c0_curr_val in self.c0_vals:
                    for approx_rank in self.approximations:
                        for curr_lin_op in self.lin_operator:
                            # print('ready to solve ODE subproblem')
                            subproblem = HamTaskSubSolution(TaskType.ODE_TASK)
                            subproblem.lin_op = curr_lin_op
                            subproblem.c0 = c0_curr_val
                            subproblem.approximation = approx_rank
                            subproblem.initial_guess = initial_guess
                            subproblem.bcs = self.bcs
                            subproblem.task_parent = self
                            subproblem.solve_kernel()
                            self.solutions.append(subproblem)
        else:
            print("Task {:d} is not marked for solution".format(self.id))

    def build_optionals(self, inp_dict):
        t_solve = inp_dict.get('solve')
        if t_solve:
            if 'n' == t_solve:
                self.to_solve = False
        t_lin_op = inp_dict.get('linear_operator')
        if t_lin_op:
            parse_possible_expr_str_iterable(t_lin_op, self.lin_operator)
        else:
            self.auto_calc_lin_op()
        t_initial_guess = inp_dict.get('initial_guess')
        if t_initial_guess:
            parse_possible_expr_str_iterable(t_initial_guess, self.initial_guesses)
        else:
            self.auto_calc_initial_guess()
        t_c0 = inp_dict.get('c0')
        if t_c0:
            parse_possible_expr_str_iterable(t_c0, self.c0_vals)
        else:
            self.c0_vals.append(sympify("-1."))
        t_approx = inp_dict.get('approximation')
        if t_approx:
            parse_possible_expr_str_iterable(t_approx, self.approximations)
        else:
            self.approximations.append(2)
        t_endtime = inp_dict.get('ending_time')
        if t_endtime:
            self.error_calculator.domain_end = float(t_endtime)
        t_eval_pnts = inp_dict.get('eval_points_num')
        if t_eval_pnts:
            self.error_calculator.points = int(t_eval_pnts)
        t_interval_step = inp_dict.get('eval_interval_step')
        if t_interval_step:
            self.error_calculator.step = float(t_interval_step)
        t_force = inp_dict.get('force')
        if t_force and str(t_force) == 'y':
            self.force = True

    def parse_bcs(self, t_bcs_obj):
        for bcs_dict in t_bcs_obj:
            if 2 == len(bcs_dict):
                a_new_bcs = BoundaryCondition(TaskType.ODE_TASK)
                for key, val in bcs_dict.items():
                    if 'x' == key:
                        a_new_bcs.x = ham_sympify(val)
                    else:
                        key_s = ham_sympify(key)
                        val_s = ham_sympify(val)
                        a_new_bcs.set_f_der_val(key_s, val_s)
                        a_new_bcs.check_type()
                self.bcs.append(a_new_bcs)
            else:
                ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')

    def build(self, inp_dict):
        t_eq = inp_dict.get('equation')
        t_bcs = inp_dict.get('boundary_conditions')
        if t_eq and t_bcs:
            try:
                self.eq = ham_sympify(t_eq)
                self.parse_bcs(t_bcs)
                self.parsing_status = True
                self.build_optionals(inp_dict)
            except SympifyError:
                self.parsing_status = False
        else:
            self.parsing_status = False

    def boundary_conditions_str(self):
        ret_str = '-------- Boundary Conditions ----------\n'
        for curr_bcs in self.bcs:
            ret_str += curr_bcs.__str__()
        ret_str += '--------------------------------------\n'
        return ret_str

    def linear_operator_str(self):
        ret_str = '-'*5 + 'Linear Operators' + '-'*5 + '\n'
        for curr_lin_op in self.lin_operator:
            ret_str += 'L: ' + str(curr_lin_op) + '\n'
        return ret_str

    def initial_guesses_str(self):
        ret_str = '-'*5 + 'Initial Guesses' + '-'*5 + '\n'
        for curr_initial_guess in self.initial_guesses:
            ret_str += '\t ' + str(curr_initial_guess) + '\n'
        return ret_str

    def c0_vals_str(self):
        ret_str = 'C0 vals: '
        for curr_C0_val in self.c0_vals:
            ret_str += str(curr_C0_val) + ', '
        ret_str += ' ]\n'
        return ret_str

    def approximations_str(self):
        ret_str = 'Approximations: '
        for curr_approx in self.approximations:
            ret_str += str(curr_approx) + ', '
        ret_str += ' ]\n'
        return ret_str

    def __str__(self):
        ret_str = ''
        ret_str += 'OBJECT: HAM Task ODE\n'
        ret_str += 'ID: ' + str(self.id) + '\n'
        if self.parsing_status:
            ret_str += 'Equation: ' + str(self.eq) + '\n'
            ret_str += self.boundary_conditions_str()
            ret_str += self.linear_operator_str()
            ret_str += self.initial_guesses_str()
            ret_str += self.c0_vals_str()
            ret_str += self.approximations_str()
            ret_str += '---------- SOLUTIONS({})-------------\n'.format(len(self.solutions))
            for a_sol in self.solutions:
                ret_str += a_sol.__str__()
            ret_str += '-------------------------------------\n'
        else:
            ret_str += '*** Failed to parse data. Empty object ***\n'
            ret_str += '='*20
        return ret_str


# -------------------------------------------------------------------------


class HamTaskPDE(HamTask):
    def __init__(self):
        super().__init__('HAM_PDE_TASK')
        self.bcs = []
        self.lin_operator = []
        self.error_calculator = PdeErrorCalculator()
        self.force = False
        self.solution_instruction = []

    def set_eq(self, expr):
        inc_expr = ham_sympify(expr)
        if not check_eq_pde_valid(inc_expr):
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        else:
            self.eq = inc_expr

    def auto_calc_initial_guess(self):
        ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__NOT_IMPLEMENTED_YET')

    def auto_calc_lin_op(self):
        ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__NOT_IMPLEMENTED_YET')

    def solve(self):
        print('initial guesses: ', self.initial_guesses)
        print('co vals: ', self.c0_vals)
        print('approx: ', self.approximations)
        print('lin_op: ', self.lin_operator)
        if self.parsing_status and self.to_solve:
            for initial_guess in self.initial_guesses:
                for c0_curr_val in self.c0_vals:
                    for approx_rank in self.approximations:
                        for curr_lin_op in self.lin_operator:
                            print('ready to solve PDE subproblem')
                            subproblem = HamTaskSubSolution(TaskType.PDE_TASK)
                            subproblem.lin_op = curr_lin_op
                            subproblem.c0 = c0_curr_val
                            subproblem.approximation = approx_rank
                            subproblem.initial_guess = initial_guess
                            subproblem.bcs = self.bcs
                            subproblem.task_parent = self
                            subproblem.solve_kernel()
                            self.solutions.append(subproblem)
        else:
            print('Task {} is not marked for solution'.format(self.id))

    @staticmethod
    def _create_point_from_str(str_line):
        ret_pnt = None
        if len(str_line) == 3:
            ret_pnt = HamPoint()
            ret_pnt.x = str_to_float(str_line[0])
            ret_pnt.y = str_to_float(str_line[1])
            ret_pnt.t = str_to_float(str_line[2])
            if not ham_any([ret_pnt.x, ret_pnt.y, ret_pnt.t]):
                del ret_pnt
                ret_pnt = None
        return ret_pnt

    def _parse_error_eval_file(self, filename):
        csv_file = open(filename)
        reader = csv.reader(csv_file)
        pnts_list = list(reader)
        for pnt in pnts_list:
            pnt_got = self._create_point_from_str(pnt)
            if isinstance(pnt_got, HamPoint):
                self.error_calculator.eval_points.append(pnt_got)

    def build_optionals(self, inp_dict):
        t_solve = inp_dict.get('solve')
        if t_solve is not None:
            if 'n' == t_solve:
                self.to_solve = False
        t_lin_op = inp_dict.get('linear_operator')
        if t_lin_op is not None:
            parse_possible_expr_str_iterable(t_lin_op, self.lin_operator)
        else:
            self.auto_calc_lin_op()
        t_initial_guess = inp_dict.get('initial_guess')
        if t_initial_guess is not None:
            parse_possible_expr_str_iterable(t_initial_guess, self.initial_guesses)
        else:
            self.auto_calc_initial_guess()
        t_c0 = inp_dict.get('c0')
        if t_c0 is not None:
            parse_possible_expr_str_iterable(t_c0, self.c0_vals)
        else:
            self.c0_vals.append(sympify("-1."))
        t_approx = inp_dict.get("approximation")
        if t_approx is not None:
            parse_possible_expr_str_iterable(t_approx, self.approximations)
        else:
            self.approximations.append(2)
        t_sol_instr = inp_dict.get('solution_instruction')
        if t_sol_instr:
            parse_possible_expr_str_iterable(t_sol_instr, self.solution_instruction)
        t_force = inp_dict.get('force')
        if t_force and str(t_force) == 'y':
            self.force = True
        t_error_eval_file = inp_dict.get('error_eval_file')
        if t_error_eval_file:
            if is_csv_file(t_error_eval_file):
                self._parse_error_eval_file(t_error_eval_file)
            pass

    def _parse_a_bcs(self, bcs_dict):
        ret_bcs = BCPDE()
        if "boundary" in bcs_dict.keys():
            boundary_expr = ham_sympify(bcs_dict["boundary"])
            ret_bcs.set_boundary(boundary_expr)
        if "value" in bcs_dict.keys():
            value_expr = ham_sympify(bcs_dict["value"])
            ret_bcs.set_value(value_expr)
        print(ret_bcs)
        return ret_bcs

    def parse_bcs(self, t_bcs_obj):
        for bcs_dict in t_bcs_obj:
            if 2 == len(bcs_dict):
                a_new_bcs = self._parse_a_bcs(bcs_dict)
                if a_new_bcs.is_well_defined():
                    self.bcs.append(a_new_bcs)
            else:
                ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')

    def build(self, inp_dict):
        t_eq = inp_dict.get('equation')
        t_bcs = inp_dict.get('boundary_conditions')
        if t_eq is not None and t_bcs is not None:
            try:
                self.eq = ham_sympify(t_eq)
                self.parse_bcs(t_bcs)
                self.parsing_status = True
                self.build_optionals(inp_dict)
            except SympifyError:
                self.parsing_status = False
        else:
            self.parsing_status = False

    def boundary_conditions_str(self):
        ret_str = '-------------- Boundary Conditions ---------\n'
        for curr_bcs in self.bcs:
            ret_str += curr_bcs.__str__()
        ret_str += '-------------------------------------------\n'
        return ret_str

    def linear_operator_str(self):
        ret_str = '-'*5 + 'Linear Operator' + '-'*5 + '\n'
        for curr_lin_op in self.lin_operator:
            ret_str += 'L: ' + str(curr_lin_op) + '\n'
        return ret_str

    def initial_guesses_str(self):
        ret_str = '-'*5 + 'Linear Operator' + '-'*5 + '\n'
        for curr_initial_guess in self.initial_guesses:
            ret_str += '\t' + str(curr_initial_guess) + '\n'
        return ret_str

    def c0_vals_str(self):
        ret_str = 'C0 vals: '
        for curr_c0_val in self.c0_vals:
            ret_str += str(curr_c0_val) + ', '
        ret_str += ' ]\n'
        return ret_str

    def approximations_str(self):
        ret_str = 'Approximations: '
        for curr_approx in self.approximations:
            ret_str += str(curr_approx) + ', '
        ret_str += ' ]\n'
        return ret_str

    def __str__(self):
        ret_str = ''
        ret_str += 'OBJECT: HAM Task PDE\n'
        ret_str += 'ID: ' + str(self.id) + '\n'
        if self.parsing_status:
            ret_str += 'Equation: ' + str(self.eq) +'\n'
            ret_str += self.boundary_conditions_str()
            ret_str += self.linear_operator_str()
            ret_str += self.initial_guesses_str()
            ret_str += self.c0_vals_str()
            ret_str += self.approximations_str()
            ret_str += '------------ SOLUTIONS ({}) -----------\n'.format(len(self.solutions))
            for a_sol in self.solutions:
                ret_str += a_sol.__str__()
            ret_str += '---------------------------------------\n'
        else:
            ret_str += '*** Failed to parse data. Empty object ***\n'
            ret_str += '='*20
        return ret_str


# -------------------------------------------------------------------------

class HamTaskStefan(HamTask):
    def __init__(self):
        ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__NOT_IMPLEMENTED_YET')
        super().__init__('HAM_STEFAN_TASK')

