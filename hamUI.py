import os
import datetime
from hamBase import HamEnt
from hamErrors import ham_error
from hamParser import HamParser


def is_valid_filename(file_name):
    return file_name and os.path.exists(file_name)


def is_valid_output_filename(file_name):
    return file_name and file_name != ""


def ask_for_input_filename():
    ret_str = None
    complain = False
    while not is_valid_filename(ret_str):
        if complain:
            ham_error('HAM_ERROR_HANDLED', 'HANDLED__INVALID_INPUT_FILENAME_TERMINAL')
        ret_str = input("pyHAM: provide a valid path for the input file:\n")
        complain = True
    return ret_str


def generate_current_filename():
    ret_str = str(datetime.datetime.now())
    ret_str = ret_str.replace(' ', '_')
    ret_str = ret_str.replace(':', '_')
    ret_str = ret_str.replace('.', '_')
    ret_str += '.txt'
    return ret_str


def ask_for_output_filename():
    output_filename = input("pyHAM: provide path for output file or leave blank:\n")
    if not is_valid_output_filename(output_filename):
        ham_error('HAM_ERROR_HANDLED', 'HANDLED__INVALID_OUTPUT_FILENAME_TERMINAL')
        output_filename = generate_current_filename()
        print("Results will be written to file: " + os.getcwd() + output_filename)
    return output_filename


class HamUI(HamEnt):
    def __init__(self, ui_type, mngr):
        if ui_type == 'HAM_UI':
            ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__VIRTUAL_METHOD_CALLED')
        super().__init__(ui_type)
        self.mngr = mngr

    def exec(self):
        ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__VIRTUAL_METHOD_CALLED')


class HamGUI(HamUI):
    def __init__(self, mngr):
        ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__NOT_IMPLEMENTED_YET')
        super().__init__('HAM_GUI', mngr)

    def exec(self):
        # TODO
        pass


class HamNoGUI(HamUI):
    def __init__(self, mngr, input_filename=None, output_filename=None):
        super().__init__('HAM_NOGUI', mngr)
        self.input_filename = input_filename
        self.output_filename = output_filename
        self.parser = HamParser()

    def exec(self):
        if not is_valid_filename(self.input_filename):
            self.input_filename = ask_for_input_filename()

        if not is_valid_filename(self.output_filename):
            self.output_filename = ask_for_output_filename()

        self.parse()
        self.solve()
        self.report()

    def parse(self):
        ret_tasks = self.parser.parse_json_file(self.input_filename)
        for a_task in ret_tasks:
            self.mngr.tasks.append(a_task)
        pass

    def solve(self):
        self.mngr.solve_tasks()

    def report(self):
        pass

