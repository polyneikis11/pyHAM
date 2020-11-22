from enum import Enum
from hamTask import HamTaskAlgebraic, HamTaskSubSolution, HamTaskODE, SubSolutionStatus
from hamCommons import *
from hamBcs import HamBcs, construct_inner_ode_bcs_from_HamBcs
from hamResult import HamResult, HamResultResolution


class HamSolver:
    def __init__(self):
        self.equation = init_expr()
        self.bcs = []
        self.ig = init_expr()
        self.lin_op = init_expr()
        self.C0 = S(-1.)
        self.ftol = 1e-6
        self.approx_order = 2
        self.ham_type = HamSolveTypes.UNKNOWN_TYPE

    def set_equation(self, eq_in):
        self.equation = h_sympify(eq_in)

    def get_equation(self):
        return str(self.equation)

    def set_initial_guess(self, ig_in):
        self.ig = h_sympify(ig_in)

    def get_initial_guess(self):
        return str(self.ig)

    def set_condition(self, ham_bc_in):
        if isinstance(ham_bc_in, HamBcs):
            if ham_bc_in.validate():
                self.bcs.append(ham_bc_in)
            else:
                raise ValueError("HamBcs object invalid.")
        else:
            raise TypeError("Invalid type: HamBcs needed.")

    def get_conditions(self):
        return self.bcs

    def clear_conditions(self):
        self.bcs.clear()

    def set_ftol(self, f_in):
        if isinstance(f_in, float):
            if f_in > 0.:
                self.ftol = f_in
            else:
                raise ValueError("Float tolerance should be possible")
        else:
            raise TypeError("Tolerance should be a float.")

    def get_ftol(self):
        return self.ftol

    def set_C0(self, C0_in):
        self.C0 = h_sympify(C0_in)

    def get_C0(self):
        return str(self.C0)

    def set_lin_op(self, lin_op_in):
        self.lin_op = h_sympify(lin_op_in)

    def get_lin_op(self):
        return str(self.lin_op)

    def set_approx_order(self, ord_in):
        ord_msg = "Approximation order should be a positive integer"
        if isinstance(ord_in, int):
            if ord_in >= 0:
                self.approx_order = ord_in
            else:
                raise ValueError(ord_msg)
        else:
            raise TypeError(ord_msg)

    def get_approx_order(self):
        return self.approx_order

    def _detect_ham_type(self):
        if is_ode_expr(self.equation):
            self.ham_type = HamSolveTypes.ODE_BVP
        elif is_algebraic_expr(self.equation):
            self.ham_type = HamSolveTypes.ALGEBRAIC_EQUATION
        else:
            self.ham_type = HamSolveTypes.UNKNOWN_TYPE

    def get_type(self):
        self._detect_ham_type()
        return str(self.ham_type)

    def _valid_ig_for_ode(self):
        return is_expr_valid(self.ig) and not is_ode_expr(self.ig)

    def _valid_ig_for_algebraic(self):
        return self._valid_ig_for_ode() and not is_algebraic_expr(self.ig)
        pass

    def _valid_linear_operator(self):
        return is_expr_valid(self.lin_op) and is_ode_expr(self.lin_op)
        pass

    def _valid_sufficient_bcs(self):
        if HamSolveTypes.ODE_BVP == self.ham_type:
            max_order = 0
            f_symbols = self.equation.atoms(Symbol)
            for current_symbol in f_symbols:
                if current_symbol in ODE_F_EXPRESSIONS:
                    current_index = ODE_F_EXPRESSIONS.index(current_symbol)
                    if current_index > max_order:
                        max_order = current_index
            ans = (len(self.bcs) == max_order)
        else:
            raise TypeError("Equation should be of ODE type.")
        return ans

    def _validate_inner(self):
        ret_err_code = HamErrorCode.HAM_OK
        self._detect_ham_type()
        if HamSolveTypes.UNKNOWN_TYPE == self.ham_type:
            ret_err_code = HamErrorCode.HAM_UNKNOWN_EQUATION_TYPE
        elif HamSolveTypes.ODE_BVP == self.ham_type:
            if not self._valid_ig_for_ode():
                ret_err_code = HamErrorCode.HAM_ODE_INVALID_INITIAL_GUESS
            elif not self._valid_linear_operator():
                ret_err_code = HamErrorCode.HAM_ODE_WRONG_LINEAR_OPERATOR
            elif not self._valid_sufficient_bcs():
                ret_err_code = HamErrorCode.HAM_ODE_WRONG_NUMBER_OF_BCS
        elif HamSolveTypes.ALGEBRAIC_EQUATION == self.ham_type:
            if not self._valid_ig_for_algebraic():
                ret_err_code = HamErrorCode.HAM_ALGEBRAIC_INVALID_INITIAL_GUESS
        return ret_err_code

    def validate(self, display_error=False, raise_exc=False):
        err_code = self._validate_inner()
        ans = HamErrorCode.HAM_OK == err_code
        if not ans:
            if display_error:
                print(ham_error_code_to_msg(err_code))
            if raise_exc:
                raise TypeError(ham_error_code_to_msg(err_code))
        return ans

    def _solve_algebraic(self):
        task = HamTaskAlgebraic()
        task.set_eq(self.equation)
        task.initial_guesses.append(self.ig)
        task.c0_vals.append(self.C0)
        task.approximations.append(self.approx_order)
        task.set_parsing_status(True)
        task.to_solve = True
        task.solve()
        ret_result = HamResult()
        if 1 == len(task.solutions):
            solution_obj = task.solutions[0]  # type:HamTaskSubSolution
            ret_result.value = solution_obj.answer.value
            ret_result.error = solution_obj.answer.error
            ret_result.tolerance = self.ftol
            ret_result.ham_type = self.ham_type
            ret_result.eq = self.equation
            if abs(ret_result.error) > ret_result.tolerance:
                ret_result.status = HamResultResolution.SOLVED_WITH_ERROR_OVER_TOLERANCE
            else:
                ret_result.status = HamResultResolution.SOLVED
        else:
            ret_result.status = HamResultResolution.NOT_RUN
            msg = 'PyHAM: Error in generating solution for this case.\n'
            msg += 'Please report bug to ' + pyham_bug_report_email
            print(msg)
        return ret_result

    def _solve_ode_bvp(self):
        task = HamTaskODE()
        task.set_parsing_status(True)
        task.to_solve = True
        task.set_eq(self.equation)
        task.initial_guesses.append(self.ig)
        task.approximations.append(self.approx_order)
        task.c0_vals.append(self.C0)
        for current_bcs_obj in self.bcs:
            a_new_bcs = construct_inner_ode_bcs_from_HamBcs(current_bcs_obj)
            task.bcs.append(a_new_bcs)
        task.lin_operator.append(self.lin_op)
        task.solve()
        ret_result = HamResult()
        if 1 == len(task.solutions):
            solution_obj = task.solutions[0]
            ret_result.value = solution_obj.answer.value
            ret_result.error = solution_obj.answer.error
            ret_result.tolerance = self.ftol
            ret_result.ham_type = self.ham_type
            ret_result.eq = self.equation
            for current_bcs_obj in self.bcs:
                copy_bcs = HamBcs()
                copy_bcs.value = current_bcs_obj.value
                copy_bcs.der_order = current_bcs_obj.der_order
                copy_bcs.point = current_bcs_obj.point
                ret_result.bcs.append(copy_bcs)
            if solution_obj.answer.status == SubSolutionStatus.SOLVED:
                if abs(ret_result.error) < ret_result.tolerance:
                    ret_result.status = HamResultResolution.SOLVED
                else:
                    ret_result.status = HamResultResolution.SOLVED_WITH_ERROR_OVER_TOLERANCE
            else:
                ret_result.status = HamResultResolution.NOT_SOLVED
        else:
            ret_result.status = HamResultResolution.NOT_RUN
            msg = 'PyHAM: Error in generating solution for this case.\n'
            msg += 'Please report bug to ' + pyham_bug_report_email
            print(msg)
        return ret_result

    def solve(self):
        ret_result = None
        if self.validate():
            if HamSolveTypes.ALGEBRAIC_EQUATION == self.ham_type:
                ret_result = self._solve_algebraic()
            elif HamSolveTypes.ODE_BVP == self.ham_type:
                ret_result = self._solve_ode_bvp()
        return ret_result

