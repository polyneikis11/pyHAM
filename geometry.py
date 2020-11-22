from hamSymbols import *
from hamBase import HamEnt
from hamErrors import ham_error


def str_to_float(str_in):
    try:
        res = float(str_in)
    except ValueError:
        res = None
    return res


class HamPoint:
    def __init__(self):
        self.x = None
        self.y = None
        self.t = None

    def __str__(self):
        x_to_print = 'n'
        y_to_print = 'n'
        t_to_print = 't'
        if self.x is not None:
            x_to_print = self.x
        if self.y is not None:
            y_to_print = self.y
        if self.t is not None:
            t_to_print = self.t
        return '({}, {}, {})'.format(x_to_print, y_to_print, t_to_print)

    def p_list(self, symb_in):
        symb_out = []
        for symb in symb_in:
            if x == symb:
                symb_out.append(self.x)
            elif y == symb:
                symb_out.append(self.y)
            elif t == symb:
                symb_out.append(self.t)
            else:
                ham_error('HAM_ERROR_HANDLED', 'HANDLED__SOFT_ASSERT')
        return symb_out

# -------------------------------------------------------------------------------------------


class HamManifold(HamEnt):
    def __init__(self):
        super().__init__('HAM_MANIFOLD')
        self.dims = None

