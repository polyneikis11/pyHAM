from functools import total_ordering
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation as animation
from sympy import lambdify, simplify
from enum import Enum
from hamBase import HamEnt, ham_not_none
from hamSymbols import *
from geometry import HamPoint
from hamErrors import ham_error


class PlotType(Enum):
    UNAVAILABLE = 0,
    PLOT = 1,
    CONTOUR_PLOT = 2,
    ANIMATED_PLOT = 3,
    ANIMATED_CONTOUR_PLOT = 4


def decide_plot_type(points):
    ret_type = PlotType.UNAVAILABLE
    x_absent = False
    y_absent = False
    t_absent = False
    for pnt in points:  # type: HamPoint
        if pnt.x is None:
            x_absent = True
            break
        if pnt.y is None:
            y_absent = True
        else:
            if y_absent:
                x_absent = True
                break  # mixed point case
        if pnt.t is None:
            t_absent = True
        else:
            if t_absent:
                x_absent = True
                break
    if not x_absent:
        if y_absent:
            if t_absent:
                ret_type = PlotType.PLOT
            else:
                ret_type = PlotType.ANIMATED_PLOT
        else:
            if t_absent:
                ret_type = PlotType.CONTOUR_PLOT
            else:
                ret_type = PlotType.ANIMATED_CONTOUR_PLOT
    return ret_type


def vars_by_plot_type(plot_type):
    ret_vars = None
    if plot_type == PlotType.PLOT:
        ret_vars = [x]
    elif plot_type == PlotType.CONTOUR_PLOT:
        ret_vars = [x, y]
    elif plot_type == PlotType.ANIMATED_PLOT:
        ret_vars = [x, t]
    elif plot_type == PlotType.ANIMATED_CONTOUR_PLOT:
        ret_vars = [x, y, t]
    else:
        ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
    return ret_vars


@total_ordering
class TimeData:
    def __init__(self):
        self.time = None
        self.data = []

    def __eq__(self, other):
        return self.time == other.time

    def __lt__(self, other):
        return self.time < other.time


# -----------------------------------------------------
class HamPlotter(HamEnt):
    def __init__(self):
        super().__init__('HAM_PLOTTER')
        self.style = None
        self.fun = None
        self.points = None
        self.parametric = None
        self.status = True
        self.plt_type = PlotType.UNAVAILABLE
        self.data = []
        self.ani_data = []
        self.computed_ani_data = []
        self.current_plot = None
        self.iter_in_ani_plot = 0
        self.max_iters_in_ani_plot = -1

    def _reset(self):
        self.status = True
        self.plt_type = PlotType.UNAVAILABLE
        self.data.clear()
        self.ani_data.clear()
        self.computed_ani_data.clear()
        self.current_plot = None
        self.iter_in_ani_plot = 0
        self.max_iters_in_ani_plot = -1

    def _set_fun(self, fun):
        self.fun = fun

    def _set_points(self, points):
        self.points = points

    def _set_parametric(self, parametric):
        self.parametric = parametric

    def _preliminary_checks(self):
        self.plt_type = decide_plot_type(self.points)
        if self.plt_type != PlotType.UNAVAILABLE:
            vars_in = vars_by_plot_type(self.plt_type)
            vars_in_expr = self.fun.atoms(Symbol)
            bad_gen = (var for var in vars_in_expr if var not in vars_in)
            for var in bad_gen:
                self.status = False
                ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
                break
        else:
            # TODO : bad request for plot. Mention back
            self.status = False
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')

    def _build_static_data(self):
        for point in self.points:  # type: HamPoint
            a_list = [point.x]
            if ham_not_none(point.y):
                a_list.append(point.y)
            self.data.append(a_list)

    def _compute_ani_data(self, ani_data_in):
        if PlotType.ANIMATED_PLOT == self.plt_type:
            the_vars = [x]
        else:
            the_vars = [x, y]
        curr_time = ani_data_in[0].time
        fun_to_use = simplify(self.fun.subs({t: curr_time}))
        u_l = lambdify(the_vars, fun_to_use, modules=['math', 'sympy'])
        x_vals = []
        u_vals = []
        if PlotType.ANIMATED_PLOT == self.plt_type:
            for t_data in ani_data_in:  # type: TimeData
                x_vals.append(t_data.data[0])
                u_vals.append(u_l(t_data.data[0]))
            x_vals_to_ret = np.array(x_vals)
            u_vals_to_ret = np.array(u_vals)
            ret_val = (x_vals_to_ret, u_vals_to_ret)
        else:
            # TODO fill when ready for contour plots
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        return ret_val

    def _build_ani_data(self):
        all_t_data = []
        for point in self.points:  # type: HamPoint
            t_data = TimeData()
            t_data.time = point.t
            t_data.data = [point.x]
            if ham_not_none(point.y):
                t_data.data.append(point.y)
            all_t_data.append(t_data)
        all_t_data.sort()
        curr_time = None
        curr_list = None
        for curr_t_data in all_t_data:
            if curr_t_data.time == curr_time:
                curr_list.append(curr_t_data)
            else:
                curr_time = curr_t_data.time
                curr_list = [curr_t_data]
                self.ani_data.append(curr_list)
        for curr_ani_data in self.ani_data:
            self.computed_ani_data.append(self._compute_ani_data(curr_ani_data))

    def _is_ani_plot(self):
        return PlotType.ANIMATED_PLOT == self.plt_type or PlotType.ANIMATED_CONTOUR_PLOT == self.plt_type

    def _build_data(self):
        if self.status:
            if self._is_ani_plot():
                self._build_ani_data()
            else:
                self._build_static_data()

    def _plot_plot(self):
        x_vals = []
        u_vals = []
        u_l = lambdify(x, self.fun, modules=['math', 'sympy'])
        for data_point in self.data:
            x_vals.append(data_point[0])
            u_vals.append(u_l(data_point[0]))
        x_vals_np = np.array(x_vals)
        u_vals_np = np.array(u_vals)
        plt.plot(x_vals_np, u_vals_np)
        plt.show()

    def _plot_cplot(self):
        ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__NOT_IMPLEMENTED_YET')
        pass

    def _update_for_ani(self, i):
        if self.iter_in_ani_plot < self.max_iters_in_ani_plot:
            self.current_plot.set_data(self.computed_ani_data[self.iter_in_ani_plot][0], self.computed_ani_data[self.iter_in_ani_plot][1])
            self.iter_in_ani_plot += 1
        return self.current_plot,

    def _get_lims_of_ani_data(self):
        xmin = np.amin(self.computed_ani_data[0][0])
        xmax = np.amax(self.computed_ani_data[0][0])
        umin = np.amin(self.computed_ani_data[0][1])
        umax = np.amax(self.computed_ani_data[0][1])
        for i in range(len(self.computed_ani_data)):
            temp_x_array = self.computed_ani_data[i][0]
            temp_u_array = self.computed_ani_data[i][1]
            temp_x_min = np.amin(temp_x_array)
            temp_x_max = np.amax(temp_x_array)
            temp_u_min = np.amin(temp_u_array)
            temp_u_max = np.amax(temp_u_array)
            if temp_x_min < xmin:
                xmin = temp_x_min
            if temp_x_max > xmax:
                xmax = temp_x_max
            if temp_u_min < umin:
                umin = temp_u_min
            if temp_u_max > umax:
                umax = temp_u_max
        return xmin, xmax, umin, umax

    def _plot_ani_plot(self):
        fig, ax = plt.subplots()
        line, = ax.plot(self.computed_ani_data[0][0], self.computed_ani_data[0][1])
        xmin, xmax, umin, umax = self._get_lims_of_ani_data()
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(umin, umax)
        self.current_plot = line
        self.max_iters_in_ani_plot = len(self.computed_ani_data)
        ani = animation.FuncAnimation(fig, self._update_for_ani, interval=1000, blit=True, frames=self.max_iters_in_ani_plot, repeat=False)
        # plt.show()
        ani.save('test.gif', writer='imagemagick', fps=1)

    def _plot_ani_cplot(self):
        ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__NOT_IMPLEMENTED_YET')
        pass

    def _plot(self):
        if self.status:
            if PlotType.PLOT == self.plt_type:
                self._plot_plot()
            elif PlotType.CONTOUR_PLOT == self.plt_type:
                self._plot_cplot()
            elif PlotType.ANIMATED_PLOT == self.plt_type:
                self._plot_ani_plot()
            elif PlotType.ANIMATED_CONTOUR_PLOT == self.plt_type:
                self._plot_ani_cplot()
            else:
                ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')

    def do_the_plot(self, u_fun, points, parametric=None, style=None):
        if ham_not_none(parametric):
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__NOT_IMPLEMENTED_YET')
        if ham_not_none(style):
            self.style = style
        self._reset()
        self._set_fun(u_fun)
        self._set_points(points)
        self._set_parametric(parametric)
        self._preliminary_checks()
        self._build_data()
        self._plot()

