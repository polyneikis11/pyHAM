from hamCommons import h_sympify, is_expr_valid, is_obj_valid
from hamTask import BoundaryCondition, TaskType


class HamBcs:
    def __init__(self, point_in=None, value_in=None, der_order=-1):
        self.point = None
        self.value = None
        self.der_order = -1
        if is_obj_valid(point_in):
            self.point = h_sympify(point_in)
        if is_obj_valid(value_in):
            self.value = h_sympify(value_in)
        if der_order >= 0:
            self.der_order = der_order

    def set_point(self, point_str):
        self.point = h_sympify(point_str)

    def get_point(self):
        if is_expr_valid(self.point):
            ret_str = str(self.point)
        else:
            raise ValueError("Point in HamBcs is invalid.")
        return ret_str

    def set_value(self, value_str):
        self.value = h_sympify(value_str)

    def get_value(self):
        if is_expr_valid(self.value):
            ret_str = str(self.value)
        else:
            raise ValueError("Value in HamBcs is invalid.")
        return ret_str

    def get_der_order(self):
        if self.der_order >=0:
            return self.der_order
        else:
            raise ValueError("Order of derivative in HamBcs is invalid.")

    def _report_error(self):
        print('HamBcs validation report:')
        if not is_expr_valid(self.point):
            print('Invalid point of HamBcs.')
        if not is_expr_valid(self.value):
            print('Invalid value of HamBcs.')
        if self.der_order < 0:
            print('Invalid order of derivative of HamBcs.')
        print('-------------------------')

    def validate(self, display_error=False):
        ans = self.der_order >= 0 and is_expr_valid(self.point) and is_expr_valid(self.value)
        if display_error:
            self._report_error()
        return ans


def construct_inner_ode_bcs_from_HamBcs(ham_bcs_obj):
    if isinstance(ham_bcs_obj, HamBcs):
        a_new_bcs = BoundaryCondition(TaskType.ODE_TASK)
        a_new_bcs.x = ham_bcs_obj.point
        a_new_bcs.value = ham_bcs_obj.value
        a_new_bcs.fx_order = ham_bcs_obj.der_order
        return a_new_bcs
    else:
        raise TypeError("Object is not of HamBcs type.")




