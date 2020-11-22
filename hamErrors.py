MAX_APPROX_ALGEBRAIC = 10

HAM_ERROR_RANKS = [
    'HAM_ERROR_CRITICAL',
    'HAM_ERROR_HANDLED',
    'HAM_WARNING'
]

HAM_CRITICAL_ERRORS_TYPOLOGY = {
    'CRITICAL__ASSERT_NOT_REACHED': 'Assert not reached',
    'CRITICAL__BAD_ERROR_RANK': 'Invalid error rank given',
    'CRITICAL__BAD_ERROR_TYPE': 'Invalid error type given',
    'CRITICAL__BAD_HAMTYPE': 'Invalid ham entity type to create',
    'CRITICAL__VIRTUAL_METHOD_CALLED': 'Method of base virtual class called',
    'CRITICAL__NOT_IMPLEMENTED_YET': 'Not implemented yet. Abort',
    'CRITICAL__BAD_ALGEBRAIC_EQUATION': 'Invalid algebraic equation to solve.'
}

HAM_HANDLED_ERRORS_TYPOLOGY = {
    'HANDLED__NOT_IMPLEMENTED_YET': 'Not implemented yet. Fall through...',
    'HANDLED__INVALID_INPUT_FILENAME_SYSARGS': 'Input filename provided in script parameters is invalid',
    'HANDLED__INVALID_OUTPUT_FILENAME_SYSARGS': 'Output filename provided in script parameters is invalid',
    'HANDLED__INVALID_INPUT_FILENAME_TERMINAL': 'Input filename provided is invalid.',
    'HANDLED__INVALID_OUTPUT_FILENAME_TERMINAL': 'Output filename provided is invalid',
    'HANDLED_INVALID_APPROX_RANK_GIVEN_NON_POSITIVE_INTEGER': '''
    Approximation rank given is invalid. It should be a positive integer.
    Default value will be used
    ''',
    'HANDLED_INVALID_APPROX_RANK_GIVEN_BIGGER_THAN_MAX': '''
    Approximation rank given is invalid. It should be less than MAX available = 
    ''' + str(MAX_APPROX_ALGEBRAIC),
    'HANDLED_REQUEST_FOR_NEGATIVE_APPROXIMATIONS_NUMBER': 'Approximation number to solution negative',
    'HANDLED_GENERAL': 'General Error, to be explained',
    'HANDLED__SOFT_ASSERT': 'This should not happen here, but we continue'
}

HAM_WARNINGS_TYPOLOGY = {

}


class HamException(Exception):
    def __init__(self, msg):
        super().__init__(msg)


class HamError:
    def __init__(self):
        self.err_rank = None
        self.err_type = None
        self.err_msg = None
        self.ui = None

    def propagate(self):
        raise HamException("Should not be reached")


class HamErrorCritical(HamError):
    def __init__(self, error_type, error_msg, ui=None):
        super().__init__()
        self.err_rank = 'HAM_ERROR_CRITICAL'
        self.err_type = error_type
        self.err_msg = error_msg
        self.ui = ui

    def propagate(self):
        exception_str = self.err_rank + ', Type: ' + self.err_type
        raise HamException(exception_str)
        # TODO: if ui , ui.report_error(self)


class HamErrorHandled(HamError):
    def __init__(self, error_type, error_msg, ui=None):
        super().__init__()
        self.err_rank = 'HAM_ERROR_HANDLED'
        self.err_type = error_type
        self.err_msg = error_msg
        self.ui = ui

    def propagate(self):
        pass
        # TODO: if ui, ui.report_error(self)


class HamWarning(HamError):
    def __init__(self, error_type, error_msg, ui=None):
        super().__init__()
        self.err_rank = 'HAM_WARNING'
        self.err_type = error_type
        self.err_msg = error_msg
        self.ui = ui

    def propagate(self):
        pass
        # TODO: if ui, ui.report_error(self)


def valid_err_type(err_type, err_rank):
    ans = False
    if err_rank == 'HAM_ERROR_CRITICAL':
        ans = err_type in HAM_CRITICAL_ERRORS_TYPOLOGY
    elif err_rank == 'HAM_ERROR_HANDLED':
        ans = err_type in HAM_HANDLED_ERRORS_TYPOLOGY
    elif err_rank == 'HAM_WARNING':
        ans = err_type in HAM_WARNINGS_TYPOLOGY
    return ans


def ham_error(rank, err_type, err_msg='', ui=None):
    """
    external routine to report an error in PyHAM
    :param rank: significance of an error. Error-Handled-Warning
    :param err_type: see Typologies of errors per rank in this file
    :param err_msg: additional message to be displayed beyond the default by type. (default empty string)
    :param ui: specify if there is a HamUI object to handle the error report (default None)
    :return:
    """
    error_obj = None
    msg_to = ''
    if valid_err_type(err_type, rank):
        if rank == 'HAM_ERROR_CRITICAL':
            msg_to += HAM_CRITICAL_ERRORS_TYPOLOGY[err_type] + err_msg
            error_obj = HamErrorCritical(err_type, msg_to, ui)
        elif rank == 'HAM_ERROR_HANDLED':
            msg_to += HAM_HANDLED_ERRORS_TYPOLOGY[err_type] + err_msg
            error_obj = HamErrorHandled(err_type, msg_to, ui)
        elif rank == 'HAM_WARNING':
            msg_to += HAM_WARNINGS_TYPOLOGY[err_type] + err_msg
            error_obj = HamWarning(err_type, msg_to, ui)

        if error_obj:
            error_obj.propagate()
        else:
            in_err = 'CRITICAL__BAD_ERROR_RANK'
            in_err_msg = HAM_CRITICAL_ERRORS_TYPOLOGY[in_err]
            error_obj = HamErrorCritical(in_err, in_err_msg, ui)
            error_obj.propagate()
    else:
        in_err = 'CRITICAL__BAD_ERROR_TYPE'
        in_err_msg = HAM_CRITICAL_ERRORS_TYPOLOGY[in_err]
        error_obj = HamErrorCritical(in_err, in_err_msg, ui)
        error_obj.propagate()

