from hamUI import HamGUI, HamNoGUI, is_valid_filename, is_valid_output_filename
from hamErrors import ham_error


class HamManager:
    def __init__(self, args):
        self.args = args
        self.ui = None
        self.tasks = []

        self.configure()

    def configure(self):
        build_gui = True
        input_f_name = None
        output_f_name = None
        for sys_arg in self.args[1:]:
            if sys_arg == '--nogui':
                build_gui = False
            elif sys_arg.startswith('--input_file='):
                tokens = sys_arg.split('=')
                if len(tokens) > 1:
                    input_f_name = tokens[1]
                    if not is_valid_filename(input_f_name):
                        ham_error('HAM_ERROR_HANDLED', 'HANDLED_INVALID_INPUT_FILENAME_SYSARGS')
            elif sys_arg.startswith('--output_file='):
                tokens = sys_arg.split('=')
                if len(tokens) > 1:
                    output_f_name = tokens[1]
                    if not is_valid_output_filename(output_f_name):
                        ham_error('HAM_ERROR_HANDLED', 'HANDLED_INVALID_OUTPUT_FILENAME_SYSARGS')

        if build_gui:
            self.ui = HamGUI(self)
        else:
            self.ui = HamNoGUI(self, input_filename=input_f_name, output_filename=output_f_name)

    def execute(self):
        self.ui.exec()

    def solve_tasks(self):
        for a_task in self.tasks:
            a_task.solve()
