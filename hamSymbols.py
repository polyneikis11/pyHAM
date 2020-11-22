import re
from sympy import Symbol, Function, Derivative, S, sympify
from math import factorial


def assure_symbols_containers_sizes():
    xterms_len = len(XTERMS_LIST)
    assert xterms_len == len(FDER_EVAL_LIST)
    assert len(ALG_LEXICON) == 3*(xterms_len - 1) - 1


# INITIATING SYMBOLS
ALL_SYMBOLS = []
HAM_LOCALS = {}

q = Symbol('q', real=True)  # homotopy parameter , q in [0, 1]
C0 = Symbol('C0', real=True)  # convergence-control parameter
h = Function('h')
f = Function('f', real=True)
x = Symbol('x', real=True, positive=True)
y = Symbol('y', real=True)
z = Symbol('z', real=True)
t = Symbol('t', real=True)
u = Function('u', real=True)
R = Function('R', real=True)

BASIC_SYMBOLS = [
    q,
    C0,
    x,
    y,
    z,
    t
]

# ====================================================================

# GLOBAL PARAMETERS
p1 = Symbol('p1', real=True)
p2 = Symbol('p2', real=True)
p3 = Symbol('p3', real=True)
p4 = Symbol('p4', real=True)
p5 = Symbol('p5', real=True)
p6 = Symbol('p6', real=True)
p7 = Symbol('p7', real=True)
p8 = Symbol('p8', real=True)
p9 = Symbol('p9', real=True)
p10 = Symbol('p10', real=True)
p11 = Symbol('p11', real=True)
p12 = Symbol('p12', real=True)
p13 = Symbol('p13', real=True)
p14 = Symbol('p14', real=True)
p15 = Symbol('p15', real=True)
p16 = Symbol('p16', real=True)
p17 = Symbol('p17', real=True)
p18 = Symbol('p18', real=True)
p19 = Symbol('p19', real=True)
p20 = Symbol('p20', real=True)
p21 = Symbol('p21', real=True)
p22 = Symbol('p22', real=True)
p23 = Symbol('p23', real=True)
p24 = Symbol('p24', real=True)
p25 = Symbol('p25', real=True)
p26 = Symbol('p26', real=True)
p27 = Symbol('p27', real=True)
p28 = Symbol('p28', real=True)
p29 = Symbol('p29', real=True)
p30 = Symbol('p30', real=True)
p31 = Symbol('p31', real=True)
p32 = Symbol('p32', real=True)
p33 = Symbol('p33', real=True)
p34 = Symbol('p34', real=True)
p35 = Symbol('p35', real=True)
p36 = Symbol('p36', real=True)
p37 = Symbol('p37', real=True)
p38 = Symbol('p38', real=True)
p39 = Symbol('p39', real=True)
p40 = Symbol('p40', real=True)

GLOBAL_PARAMETERS = [
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    p7,
    p8,
    p9,
    p10,
    p11,
    p12,
    p13,
    p14,
    p15,
    p16,
    p17,
    p18,
    p19,
    p20,
    p21,
    p22,
    p23,
    p24,
    p25,
    p26,
    p27,
    p28,
    p29,
    p30,
    p31,
    p32,
    p33,
    p34,
    p35,
    p36,
    p37,
    p38,
    p39,
    p40
]

C1 = Symbol('C1')
C2 = Symbol('C2')
C3 = Symbol('C3')
C4 = Symbol('C4')
C5 = Symbol('C5')
C6 = Symbol('C6')
C7 = Symbol('C7')
C8 = Symbol('C8')
C9 = Symbol('C9')
C10 = Symbol('C10')
C11 = Symbol('C11')
C12 = Symbol('C12')
C13 = Symbol('C13')
C14 = Symbol('C14')
C15 = Symbol('C15')
C16 = Symbol('C16')
C17 = Symbol('C17')
C18 = Symbol('C18')
C19 = Symbol('C19')
C20 = Symbol('C20')

INTEGRATION_CONSTANTS = [
    C1,
    C2,
    C3,
    C4,
    C5,
    C6,
    C7,
    C8,
    C9,
    C10,
    C11,
    C12,
    C13,
    C14,
    C15,
    C16,
    C17,
    C18,
    C19,
    C20
]
# ======================================================================

# ALGEBRAIC EQUATIONS SYMBOLS
x0t = Symbol('x0t', real=True)
x1t = Symbol('x1t', real=True)
x2t = Symbol('x2t', real=True)
x3t = Symbol('x3t', real=True)
x4t = Symbol('x4t', real=True)
x5t = Symbol('x5t', real=True)
x6t = Symbol('x6t', real=True)
x7t = Symbol('x7t', real=True)
x8t = Symbol('x8t', real=True)
x9t = Symbol('x9t', real=True)
x10t = Symbol('x10t', real=True)
x11t = Symbol('x11t', real=True)
x12t = Symbol('x12t', real=True)
x13t = Symbol('x13t', real=True)
x14t = Symbol('x14t', real=True)
x15t = Symbol('x15t', real=True)
x16t = Symbol('x16t', real=True)
x17t = Symbol('x17t', real=True)
x18t = Symbol('x18t', real=True)
x19t = Symbol('x19t', real=True)
x20t = Symbol('x20t', real=True)


fh0 = Symbol('fh0', real=True)
fh1 = Symbol('fh1', real=True)
fh2 = Symbol('fh2', real=True)
fh3 = Symbol('fh3', real=True)
fh4 = Symbol('fh4', real=True)
fh5 = Symbol('fh5', real=True)
fh6 = Symbol('fh6', real=True)
fh7 = Symbol('fh7', real=True)
fh8 = Symbol('fh8', real=True)
fh9 = Symbol('fh9', real=True)
fh10 = Symbol('fh10', real=True)
fh11 = Symbol('fh11', real=True)
fh12 = Symbol('fh12', real=True)
fh13 = Symbol('fh13', real=True)
fh14 = Symbol('fh14', real=True)
fh15 = Symbol('fh15', real=True)
fh16 = Symbol('fh16', real=True)
fh17 = Symbol('fh17', real=True)
fh18 = Symbol('fh18', real=True)
fh19 = Symbol('fh19', real=True)
fh20 = Symbol('fh20', real=True)

XTERMS_LIST = [
    x0t,
    x1t,
    x2t,
    x3t,
    x4t,
    x5t,
    x6t,
    x7t,
    x8t,
    x9t,
    x10t,
    x11t,
    x12t,
    x13t,
    x14t,
    x15t,
    x16t,
    x17t,
    x18t,
    x19t,
    x20t
]

FDER_EVAL_LIST = [
    fh0,
    fh1,
    fh2,
    fh3,
    fh4,
    fh5,
    fh6,
    fh7,
    fh8,
    fh9,
    fh10,
    fh11,
    fh12,
    fh13,
    fh14,
    fh15,
    fh16,
    fh17,
    fh18,
    fh19,
    fh20
]

ALG_LEXICON = {
    Derivative(h(q), q): x1t,
    Derivative(h(q), (q, 2)): factorial(2)*x2t,
    Derivative(h(q), (q, 3)): factorial(3)*x3t,
    Derivative(h(q), (q, 4)): factorial(4)*x4t,
    Derivative(h(q), (q, 5)): factorial(5)*x5t,
    Derivative(h(q), (q, 6)): factorial(6)*x6t,
    Derivative(h(q), (q, 7)): factorial(7)*x7t,
    Derivative(h(q), (q, 8)): factorial(8)*x8t,
    Derivative(h(q), (q, 9)): factorial(9)*x9t,
    Derivative(h(q), (q, 10)): factorial(10)*x10t,
    Derivative(h(q), (q, 11)): factorial(11)*x11t,
    Derivative(h(q), (q, 12)): factorial(12)*x12t,
    Derivative(h(q), (q, 13)): factorial(13)*x13t,
    Derivative(h(q), (q, 14)): factorial(14)*x14t,
    Derivative(h(q), (q, 15)): factorial(15)*x15t,
    Derivative(h(q), (q, 16)): factorial(16)*x16t,
    Derivative(h(q), (q, 17)): factorial(17)*x17t,
    Derivative(h(q), (q, 18)): factorial(18)*x18t,
    Derivative(h(q), (q, 19)): factorial(19)*x19t,
    Derivative(h(q), (q, 20)): factorial(20)*x20t,


    Derivative(f(h(q)), h(q)): fh1,
    Derivative(fh1, h(q)): fh2,
    Derivative(f(h(q)), (h(q), 2)): fh2,
    Derivative(fh2, h(q)): fh3,
    Derivative(f(h(q)), (h(q), 3)): fh3,
    Derivative(fh3, h(q)): fh4,
    Derivative(f(h(q)), (h(q), 4)): fh4,
    Derivative(fh4, h(q)): fh5,
    Derivative(f(h(q)), (h(q), 5)): fh5,
    Derivative(fh5, h(q)): fh6,
    Derivative(f(h(q)), (h(q), 6)): fh6,
    Derivative(fh6, h(q)): fh7,
    Derivative(f(h(q)), (h(q), 7)): fh7,
    Derivative(fh7, h(q)): fh8,
    Derivative(f(h(q)), (h(q), 8)): fh8,
    Derivative(fh8, h(q)): fh9,
    Derivative(f(h(q)), (h(q), 9)): fh9,
    Derivative(fh9, h(q)): fh10,
    Derivative(f(h(q)), (h(q), 10)): fh10,
    Derivative(fh10, h(q)): fh11,
    Derivative(f(h(q)), (h(q), 11)): fh11,
    Derivative(fh11, h(q)): fh12,
    Derivative(f(h(q)), (h(q), 12)): fh12,
    Derivative(fh12, h(q)): fh13,
    Derivative(f(h(q)), (h(q), 13)): fh13,
    Derivative(fh13, h(q)): fh14,
    Derivative(f(h(q)), (h(q), 14)): fh14,
    Derivative(fh14, h(q)): fh15,
    Derivative(f(h(q)), (h(q), 15)): fh15,
    Derivative(fh15, h(q)): fh16,
    Derivative(f(h(q)), (h(q), 16)): fh16,
    Derivative(fh16, h(q)): fh17,
    Derivative(f(h(q)), (h(q), 17)): fh17,
    Derivative(fh17, h(q)): fh18,
    Derivative(f(h(q)), (h(q), 18)): fh18,
    Derivative(fh18, h(q)): fh19,
    Derivative(f(h(q)), (h(q), 19)): fh19,
    Derivative(fh19, h(q)): fh20,
    Derivative(f(h(q)), (h(q), 20)): fh20
}

ALG_LEXICON_SECOND = {
    f(h(q)): fh0
}

ALG_LEXICON_THIRD = {
    q: 0.
}

# ====================================================================
# ORDINARY DIFFERENTIAL EQUATIONS SYMBOLS

fx0 = Symbol('fx0')
fx1 = Symbol('fx1')
fx2 = Symbol('fx2')
fx3 = Symbol('fx3')
fx4 = Symbol('fx4')
fx5 = Symbol('fx5')
fx6 = Symbol('fx6')
fx7 = Symbol('fx7')
fx8 = Symbol('fx8')
fx9 = Symbol('fx9')

ODE_F_EXPRESSIONS = [
    fx0,
    fx1,
    fx2,
    fx3,
    fx4,
    fx5,
    fx6,
    fx7,
    fx8,
    fx9
]

# ====================================================================
# PARTIAL DIFFERENTIAL EQUATIONS SYMBOLS
PDE_U_EXPRESSIONS = []

u0 = Symbol('u0')
# ---------------
ux = Symbol('ux')
uy = Symbol('uy')
ut = Symbol('ut')
# ---------------
ux2 = Symbol('ux2')
uy2 = Symbol('uy2')
ut2 = Symbol('ut2')
uxy = Symbol('uxy')
uxt = Symbol('uxt')
uyt = Symbol('uyt')
# ---------------
ux3 = Symbol('ux3')
uy3 = Symbol('uy3')
ut3 = Symbol('ut3')
ux2y = Symbol('ux2y')
ux2t = Symbol('ux2t')
uxy2 = Symbol('uxy2')
uxt2 = Symbol('uxt2')
uy2t = Symbol('uy2t')
uyt2 = Symbol('uyt2')
# ---------------
ux4 = Symbol('ux4')
uy4 = Symbol('uy4')
ut4 = Symbol('ut4')
ux3y = Symbol('ux3y')
ux3t = Symbol('ux3t')
uxy3 = Symbol('uxy3')
uxt3 = Symbol('uxt3')
ux2y2 = Symbol('ux2y2')
ux2t2 = Symbol('ux2t2')
uy3t = Symbol('uy3t')
# =====================================================================
# STEFAN PROBLEMS SPECIAL SYMBOLS
r0t = Symbol('r0t')
r1t = Symbol('r1t')
r2t = Symbol('r2t')
r3t = Symbol('r3t')
r4t = Symbol('r4t')
r5t = Symbol('r5t')
r6t = Symbol('r6t')
r7t = Symbol('r7t')
r8t = Symbol('r8t')
r9t = Symbol('r9t')

MOVING_BOUNDARY_EXP = [
    r0t,
    r1t,
    r2t,
    r3t,
    r4t,
    r5t,
    r6t,
    r7t,
    r8t,
    r9t,
]

# =====================================================================


def build_all_symbols():
    _create_pde_symbols(4, PDE_U_EXPRESSIONS)
    for curr_symbol in BASIC_SYMBOLS:
        ALL_SYMBOLS.append(curr_symbol)
    for curr_symbol in GLOBAL_PARAMETERS:
        ALL_SYMBOLS.append(curr_symbol)
    for curr_symbol in XTERMS_LIST:
        ALL_SYMBOLS.append(curr_symbol)
    for curr_symbol in ODE_F_EXPRESSIONS:
        ALL_SYMBOLS.append(curr_symbol)
    for curr_symbol in PDE_U_EXPRESSIONS:
        ALL_SYMBOLS.append(curr_symbol)
    for curr_symbol in MOVING_BOUNDARY_EXP:
        ALL_SYMBOLS.append(curr_symbol)
    for curr_symbol in ALL_SYMBOLS:
        HAM_LOCALS[str(curr_symbol)] = curr_symbol


def free_of_symbols_except(s_expr, except_symbols_list):
    ans = True
    symbols_in_expr = s_expr.atoms(Symbol)
    for curr_symbol in ALL_SYMBOLS:
        if curr_symbol not in except_symbols_list:
            if curr_symbol in symbols_in_expr:
                ans = False
                break
    return ans


def contains_ufun_symbols(s_expr):
    ans = False
    symbols_in_expr = s_expr.atoms(Symbol)
    for curr_symbol in symbols_in_expr:
        if curr_symbol in PDE_U_EXPRESSIONS or curr_symbol in MOVING_BOUNDARY_EXP:
            ans = True
            break
    return ans


def decode_der_orders_pde(the_symb):
    ret_orders = None
    if the_symb in PDE_U_EXPRESSIONS:
        num_list = str(the_symb)
        splitted = re.split(r'[xyt]', num_list)
        x_ord = splitted[1]
        y_ord = splitted[2]
        z_ord = splitted[3]
        ret_orders = int(x_ord), int(y_ord), int(z_ord)
    return ret_orders


def orders_of_pde(pde_expr):
    the_symbols = pde_expr.atoms(Symbol)
    x_ord = -1
    y_ord = -1
    t_ord = -1
    gen = (var for var in the_symbols if var in PDE_U_EXPRESSIONS)
    for u_fun in gen:
        temp_x_ord, temp_y_ord, temp_t_ord = decode_der_orders_pde(u_fun)
        if temp_x_ord > x_ord:
            x_ord = temp_x_ord
        if temp_y_ord > y_ord:
            y_ord = temp_y_ord
        if temp_t_ord > t_ord:
            t_ord = temp_t_ord
    return x_ord, y_ord, t_ord


def _create_pde_symbols(order, container):
    for x_num in range(order + 1):
        rest_of_x = order - x_num
        for y_num in range(rest_of_x + 1):
            rest_of_xy = rest_of_x - y_num
            for t_num in range(rest_of_xy + 1):
                symbols_list = []
                symbols_list.append('u')
                symbols_list.append('x')
                symbols_list.append(str(x_num))
                symbols_list.append('y')
                symbols_list.append(str(y_num))
                symbols_list.append('t')
                symbols_list.append(str(t_num))
                symbol_str = ''.join(symbols_list)
                a_new_symbol = Symbol(symbol_str)
                container.append(a_new_symbol)


def ham_sympify(expr):
    return sympify(expr, locals=HAM_LOCALS)
