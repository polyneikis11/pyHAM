from hamSolver import HamSolver
from matplotlib import pyplot as plt
from hamBcs import HamBcs
from hamResult import HamResult


def pyham_test_algebraic_1():
    my_solver = HamSolver()
    my_solver.set_equation("exp(-x**2)*cos(4*x)")
    my_solver.set_initial_guess("0.4")
    my_solver.set_ftol(1e-8)
    my_solver.set_approx_order(2)
    my_solver.validate()
    # print('equation to solve:', my_solver.get_equation())
    # print('initial guess:', my_solver.get_initial_guess())
    # print('linear op:', my_solver.get_lin_op())
    # print('C0:', my_solver.get_C0())
    # print('ftol:', my_solver.get_ftol())
    # print('Order of Approximation:', my_solver.get_approx_order())
    # print('validation:', my_solver.validate())
    # print('type of problem:', my_solver.get_type())
    ans = my_solver.solve()
    print('value:', ans.get_value())
    print('error:', ans.get_error())
    print('tol:', ans.get_tolerance())
    print('RESOLUTION:', ans.resolution())


def pyham_test_ode_1():
    my_solver = HamSolver()
    my_solver.set_equation("fx1 - 1/x**2 + fx0/x + fx0**2")
    my_solver.set_initial_guess("-1/x**2")
    my_solver.set_ftol(1e-8)
    my_solver.set_approx_order(2)
    my_solver.set_C0(-1.)
    my_solver.set_lin_op("fx1")
    a_bc = HamBcs("1.", "-1.", 0)
    my_solver.set_condition(a_bc)
    print('validation:', my_solver.validate(display_error=True))
    print('type is', my_solver.ham_type)
    res = my_solver.solve()
    print('resolution:', res.resolution())
    print('value:', res.value)
    print('error:', res.error)
    print('tol:', res.tolerance)


def run_external_pyham_tests():
    # pyham_test_algebraic_1()
    pyham_test_ode_1()
    pass
