import os
import sys
class XenomakeError(Exception):
    def __init__(self, msg=None):
        sys.tracebacklimit = 0
        self.msg = msg or "An error occurred in Xenomake."  # Default message if none provided

    def __str__(self):
        return f"ERROR: {self.__class__.__name__}\n{self.msg}"  # Simplified error message handling


class FileWrongExtensionError(XenomakeError):
    def __init__(self, filename, expected_extension):
        super().__init__(f"File {filename} has the wrong extension. Expected: {expected_extension}")
        self.filename = filename
        self.expected_extension = expected_extension


class FileNotFoundError(XenomakeError):
    def __init__(self, path):
        super().__init__(f"File {path} not found in the system.")
        self.file_path = path


class EmptyConfigVariableError(XenomakeError):
    def __init__(self, variable_name):
        super().__init__(f"Cannot set {variable_name} to an empty list or None.")
        self.variable_name = variable_name


class NoConfigError(XenomakeError):
    def __init__(self):
        super().__init__("No projects or samples provided.")


class ConfigVariableError(XenomakeError):
    def __init__(self, variable_name, variable_value):
        super().__init__(f"Config variable '{variable_name}' has an invalid value: {variable_value}.")
        self.variable_name = variable_name
        self.variable_value = variable_value


class ConfigVariableNotFoundError(ConfigVariableError):
    def __init__(self, variable_name, variable_value):
        super().__init__(variable_name, variable_value)
        self.msg = f"{variable_name}: {variable_value} not found in the configuration."

class UnrecognizedConfigVariable(ConfigVariableError):
    def __init__(self, variable_name, variable_value):
        super().__init__(variable_name, variable_value)
        self.msg = f"{variable_name} has an unrecognized value: {variable_value}."

    def __str__(self):
        return self.msg
