import os
import logging

class XenomakeError(Exception):
    def __init__(self,msg=None):
        self.msg = msg

    def __str__(self):
        msg = 'ERROR: ' + str(self.__clas__.__name__) + '\n'

        if hasattr(self,'msg') and self.msg is not None:
            msg += self.msg


class FileWrongExtensionError(XenomakeError):
    def __init__(self,filename,expected_extension):
        self.filename = filename
        self.expected_extension = expected_extension

    def __str__(self):
        msg = super().__str__()
        msg += f'File {self.filename} has wrong extension.\n'
        msg += f'The extension should be {self.expected_extension}.\n'
        return msg


class ConfigVariableError(XenomakeError):
    def __init__(self,variable_name,variable_value):
        self.variable_name = variable_name
        self.variable_value = variable_value

class EmptyConfigVariableError(XenomakeError):
    def __init__(self, variable_name):
        self.variable_name = variable_name

    def __str__(self):
        msg = super().__str__()
        msg += f'cannot remove, or set {self.variable_name} to emtpy list, or None\n'
        msg += 'this ERROR could happen in two cases: \n'
        msg += f'1) you tried to remove a {self.variable_name}, '
        msg += f'and as a result the sample would not have'
        msg += f' any {self.variable_name} available.\n'
        msg += f'2) you tried to remove the `default` value of'
        msg += f' {self.variable_name} from the configuration.\n'
        return msg
    
class NoConfigError(XenomakeError):
    def __init__(self):
        pass
    def __str__(self):
        msg = super().__str__()
        msg += f'no projects or samples provided.\n'
        return msg

class ConfigVariableError(XenomakeError):
    def __init__(self, variable_name, variable_value):
        self.variable_name = variable_name
        self.variable_value = variable_value

class ConfigVariableNotFoundError(ConfigVariableError):
    def __str__(self):
        msg = super().__str__()
        msg += f'{self.variable_name}: {self.variable_value} not found.\n'
        return msg

class UnrecognizedConfigVariable(ConfigVariableError):
    def __str__(self):
        msg = super().__str__()
        msg += f'{self.variable_name} not recognized'
        return msg