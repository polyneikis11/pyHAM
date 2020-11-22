from hamErrors import ham_error
from hamSymbols import assure_symbols_containers_sizes as assure_fun
from hamSymbols import build_all_symbols


HAM_TYPES = [
    'BASE_ENT',
    'HAM_UI',
    'HAM_GUI',
    'HAM_NOGUI',
    'HAM_TASK',
    'HAM_ALG_TASK',
    'HAM_ODE_TASK',
    'HAM_PDE_TASK',
    'HAM_STEFAN_TASK',
    'HAM_PARSER',
    'HAM_REPORTER',
    'HAM_ERROR',
    'ODE',
    'PDE',
    'BVP',
    'MBVP',
    'HAM_BC',
    'HAM_MANIFOLD',
    'HAM_PLOTTER'
]

HAM_TREES = {}

HAM_TREES_COUNTERS = {}


def ham_not_none(obj):
    return obj is not None


def ham_any(obj_iter):
    ans = False
    for obj in obj_iter:
        if ham_not_none(obj):
            ans = True
            break
    return ans


def build_trees():
    assure_fun()
    for ham_type in HAM_TYPES:
        a_new_htree = {}
        HAM_TREES[ham_type] = a_new_htree
        HAM_TREES_COUNTERS[ham_type] = 0


def init_base():
    build_trees()
    build_all_symbols()


def is_valid_type(ham_type):
    return ham_type in HAM_TYPES


class HamEnt:
    def __init__(self, ham_type):
        self.id = 0
        self.h_type = ham_type
        if not is_valid_type(ham_type):
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__BAD_HAMTYPE')
        self.reg()

    def reg(self):
        tree_to_enter = HAM_TREES[self.h_type]
        prev_cnt = HAM_TREES_COUNTERS[self.h_type]
        self.id = prev_cnt + 1
        tree_to_enter[self.id] = self
        HAM_TREES_COUNTERS[self.h_type] = self.id

    def unreg(self):
        tree_to_leave = HAM_TREES[self.h_type]
        del tree_to_leave[self.id]

    def __del__(self):
        if self.id > 0:
            self.unreg()


