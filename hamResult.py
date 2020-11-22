from enum import Enum
import numpy as np
from sympy import lambdify
from hamCommons import *


class HamResultResolution(Enum):
    SOLVED = 0
    SOLVED_WITH_ERROR_OVER_TOLERANCE = 1
    NOT_RUN = 2
    NOT_SOLVED = 3


class HamResult:
    def __init__(self):
        self.value = None
        self.error = None
        self.tolerance = None
        self.status = HamResultResolution.NOT_RUN
        self.ham_type = HamSolveTypes.UNKNOWN_TYPE
        self.eq = None
        self.bcs = []

    def valid_result(self):
        return self.status == HamResultResolution.SOLVED or \
               self.status == HamResultResolution.SOLVED_WITH_ERROR_OVER_TOLERANCE

# TODO
    # def eval_error(self, discrete=True, l_bound=None, r_bound=None, numpoints=-1):
    #     if self.valid_result():
    #         ret_error = 0.
    #         if HamSolveTypes.ALGEBRAIC_EQUATION == self.ham_type:
    #             ret_error = self.error
    #         elif HamSolveTypes.ODE_BVP == self.ham_type:
    #             pass
    #         else:
    #             print('PyHAM:HamResult:eval_error: Something went wrong.')
    #             print('Please report bug to', pyham_bug_report_email)
    #         return ret_error
    #     else:
    #         print('HamResult not valid. Cannot evaluate error')
    #         return None

    @staticmethod
    def _resol_to_message(self, resol):
        if resol == HamResultResolution.SOLVED:
            ret_msg = 'SOLVED'
        elif resol == HamResultResolution.SOLVED_WITH_ERROR_OVER_TOLERANCE:
            ret_msg = 'Solved but the error is over the defined tolerance.'
        elif resol == HamResultResolution.NOT_CONVERGED:
            ret_msg = 'Series did not converge.'
        elif resol == HamResultResolution.NOT_RUN:
            ret_msg = 'Process did not run.'
        else:
            ret_msg = 'Invalid error code. Please report bug.'
        return ret_msg

    def resolution(self):
        return self._resol_to_message(self, self.status)

    def get_plot_data(self, l_bound=None, r_bound=None, numpoints=100, display_errors=True):
        ret_obj = None
        if HamSolveTypes.ALGEBRAIC_EQUATION == self.ham_type:
            if not is_obj_valid(l_bound) or not is_obj_valid(r_bound):
                msg = 'PyHAM:HamResult:get_plot_data: For algebraic equations'
                msg += ' l_bound and r_bound must be provided'
                if display_errors:
                    print(msg)
            else:
                if float(l_bound) > float(r_bound) or is_same_point_1d_flt(l_bound, r_bound):
                    msg = 'PyHAM:HamResult:get_plot_data:'
                    msg += ' l_bound should be less than r_bound'
                    if display_errors:
                        print(msg)
                else:
                    x_array = np.linspace(float(l_bound), float(r_bound), numpoints)
                    f_lambda = lambdify(x, self.eq, 'numpy')
                    f_array = f_lambda(x_array)
                    ret_obj = x_array, f_array

        elif HamSolveTypes.ODE_BVP == self.ham_type:
            print("Not implemented yet for ODE bvp.")
        else:
            print('PyHAM:HamResult:get_plot_data: Something went wrong.')
            print('Please report bug to', pyham_bug_report_email)
        # return a dict title: tuple (of NumPy arrays)
        return ret_obj

    def get_error(self):
        return self.error

    def get_value(self):
        return self.value

    def get_tolerance(self):
        return self.tolerance
