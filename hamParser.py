import json
from hamBase import HamEnt
from hamErrors import ham_error
from hamTask import HamTaskAlgebraic, HamTaskODE, HamTaskPDE, HamTaskStefan


class HamParser(HamEnt):
    def __init__(self):
        super().__init__('HAM_PARSER')

    def build_task_type(self, inp_dict):
        ret_task = None
        t_type = inp_dict.get('type')
        if t_type:
            if 'algebraic' == t_type:
                ret_task = HamTaskAlgebraic()
            elif 'ODE' == t_type:
                ret_task = HamTaskODE()
            elif 'PDE' == t_type:
                ret_task = HamTaskPDE()
            elif 'Stefan' == t_type:
                ret_task = HamTaskStefan
            else:
                ham_error('HAM_ERROR_CRITICAL', 'CRITICAL__ASSERT_NOT_REACHED')
        return ret_task

    def build_task(self, inp_dict):
        ret_task = self.build_task_type(inp_dict)
        if ret_task:
            ret_task.build(inp_dict)
        return ret_task

    def parse_json_file(self, inp_file):
        f = open(inp_file)
        ret_tasks = []
        input_str = f.read()
        input_json_data = json.loads(input_str)
        for inp_dict in input_json_data:
            a_new_task = self.build_task(inp_dict)
            if a_new_task:
                ret_tasks.append(a_new_task)
        f.close()
        return ret_tasks


