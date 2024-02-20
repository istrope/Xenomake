import os

class XenomakeError(Exception):
    def __init__(self,msg=None):
        self.msg = msg

    def __str__(self):
        if self.msg is None:
            return self.msg
        
        msg = 'ERROR: ' + str(self.__class__.__name__) + '\n'
        return msg + self.msg
        


class FileWrongExtensionError(XenomakeError):
    def __init__(self,filename,expected_extension):
        self.filename = filename
        self.expected_extension = expected_extension

    def __str__(self):
        msg = super().__str__()
        msg += f'File {self.filename} has wrong extension.\n'
        msg += f'The extension should be {self.expected_extension}.\n'
        return msg

class FileNotFoundError(XenomakeError):
    def __init__(self,path):
        self.file_path = path
    
    def __str__(self):
        msg = super().__str__()
        msg += f'File {self.file_path} not found in system'

class ConfigVariableError(XenomakeError):
    def __init__(self,variable_name,variable_value):
        self.variable_name = variable_name
        self.variable_value = variable_value

class EmptyConfigVariableError(XenomakeError):
    def __init__(self, variable_name):
        self.variable_name = variable_name

    def __str__(self):
        msg = super().__str__()
        msg += f'cannot set {self.variable_name} to emtpy list, or None\n'
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