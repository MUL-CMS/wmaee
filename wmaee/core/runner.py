import os
import logging
import json
import sys
from os.path import join
from os import getcwd
from wmaee.core.types import *
from wmaee.core.common import working_directory, get_configuration_directory, Shell

__END_MARK__ = '__VASP_FINISHED__'
APPLICATION_CONFIG = 'application.json'
DEFAULT_APPLICATION = 'vasp_std'
DEFAULT_HOSTNAME = 'local'
DEFAULT_PARTITION = 'default'
__LOADED_CONFIGURATION = None

logger = logging.getLogger('wmaee.core.runner')


def get_vasp_configuration(application=None, hostname=None, partition=None):
    if application is None:
        # search in environment variables
        application = os.environ.get('WMAEE_APPLICATION', default=None)
        if application is None:
            application = DEFAULT_APPLICATION
    if hostname is None:
        hostname = os.environ.get('WMAEE_HOSTNAME', default=None)
        if hostname is None:
            hostname = DEFAULT_HOSTNAME
    if partition is None:
        partition = os.environ.get('WMAEE_PARTITION', default=None)
        if partition is None:
            partition = DEFAULT_PARTITION

    configuration_file = join(get_configuration_directory(), APPLICATION_CONFIG)
    global __LOADED_CONFIGURATION
    if __LOADED_CONFIGURATION is None:
        with open(configuration_file, 'r') as config_file:
            configuration = json.load(config_file)
            __LOADED_CONFIGURATION = configuration
    else:
        configuration = __LOADED_CONFIGURATION

    try:
        current_configuration = configuration[application][hostname][partition]
    except KeyError:
        raise RuntimeError('Failed to configure')

    preamble, binary, command = current_configuration['preamble'], current_configuration['binary'], \
                                current_configuration['command']
    return preamble, command, binary


def _run_vasp_internal(directory: Optional[Directory] = None, cpus: Optional[int] = 2,
                       show_output: Optional[bool] = True, return_stdout: Optional[bool] = False,
                       application: Optional[Union[None, str]] = None, hostname: Optional[Union[None, str]] = None,
                       partition: Optional[Union[None, str]] = None):
    """
    Call VASP programm and obtain the output, by executing the contents of __VASP_PREAMBLE and __VASP_COMMAND in
    bash subprocess and piping the output
    :param directory: (str or working_directory) the directory where to execute vasp
    :param cpus: (int) the number of core used to run the VASP subprocess (default: 2)
    :param show_output: (bool) wether to print the output on sys.stdout, passed on to _read_output (default: True)
    :param return_stdout: (bool) wether to return the output as a list of strings. Passed on to _send_command (default: False)
    :return: (bool) or (bool, list of str) exitcode== 0 and output depending on the setting of propagate_stdout
    """
    logger = logging.getLogger('VASPRunner')

    if directory is None:
        directory = working_directory(getcwd())

    if isinstance(directory, working_directory):
        pass
    elif isinstance(directory, str):
        directory = working_directory(directory)

    with directory:
        preamble, command, binary = get_vasp_configuration(application=application, hostname=hostname,
                                                           partition=partition)
        command = command.format(cores=cpus, binary=binary)
        change_command = 'cd {}'.format(getcwd())

        with Shell(stdout=sys.stdout if show_output else None, stderr=sys.stderr) as shell:
            shell.run(change_command)
            shell.run(preamble)
            output = shell.run(command, return_out=return_stdout, return_exit=True)

        if return_stdout:
            _, exitcode, output = output
        else:
            _, exitcode = output

        if exitcode != 0:
            logger.warning('VASP did not execute successfully! Exit-Code: "{}"'.format(exitcode))

        return exitcode if not return_stdout else (exitcode, output)
        # Write input files


def vasp(directory: Optional[Directory] = None, cpus: Optional[int] = 2,
         show_output: Optional[bool] = True, return_stdout: Optional[bool] = False,
         application: Optional[Union[None, str]] = None, hostname: Optional[Union[None, str]] = None,
         partition: Optional[Union[None, str]] = None):
    """
    Run VASP in a given directory. The maximum number of cores is limited to six. If no directory is specified VASP
    will be executed on os.getcwd() directory
    :param directory: (str or working_directory) the directory where to execute vasp
    :param cpus: (int) the number of core used to run the VASP subprocess (default: 2)
    :param show_output: (bool) wether to print the output on sys.stdout, passed on to _read_output (default: True)
    :param return_stdout: (bool) wether to return the output as a list of strings. Passed on to _send_command (default: False)
    :param application: (str) the string of the application name (default: None)
    :param hostname: (str) the string of the current hostname (default: None)
    :param partition: (str) the partition name as a string (default: None)
    :return: (bool) or (bool, list of str) exitcode== 0 and output depending on the setting of propagate_stdout
    """
    if cpus < 0:
        cpus = 1
    else:
        cpus = int(cpus)
    return _run_vasp_internal(directory=directory, cpus=cpus, show_output=show_output, return_stdout=return_stdout,
                              application=application, hostname=hostname, partition=partition)
