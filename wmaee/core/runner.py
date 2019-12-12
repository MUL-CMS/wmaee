
import os
import logging
import shlex
import json
from subprocess import PIPE, Popen
from wmaee.core.common import working_directory, ThreadWithReturnValue, get_configuration_directory
from wmaee.core.io import StringStream
from os.path import join
from os import getcwd
from io import TextIOWrapper


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
    global  __LOADED_CONFIGURATION
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

    preamble, binary, command = current_configuration['preamble'], current_configuration['binary'], current_configuration['command']
    return preamble, command, binary

def _shell_alive(shell_handle):
    """
    Thest wether a shell is still open
    :param shell_handle: (subprocess.Popen) the process handle
    :return: (bool) a boolean flag indicating if the shell is still opened
    """
    return shell_handle and not isinstance(shell_handle.poll(), int)


def _send_command(shell, cmd, return_stdout=False, propagate_stdout=True):
    """
    Execute a command on a system shell and return the output
    :param shell: (subprocess.Popen, subprocess.PIPE, subprocess.PIPE, StringStream) the process handles and output pipes
    :param cmd: (str) the command to execute
    :param return_stdout: (bool) wether to return the commands output (default: False)
    :param propagate_stdout: (bool) wether to send the commands output to the shell_output pipe (default: True)
    :return: (bool) or (bool, list of str) exitcode== 0 and output depending on the setting of propagate_stdout
    """
    shell_handle, shell_stdin, shell_stdout, shell_output = shell
    shell_input = shell_stdin
    if not _shell_alive(shell_handle):
        raise RuntimeError('Shell is not open')

    finish = '__local_command_finish_mark__'
    # self._shell_stdin.write(cmd + '\n')
    echo_cmd = 'echo {} $?\n'.format(finish)
    command = cmd.strip('\n')
    shell_input.write('; '.join([cmd, echo_cmd]))
    # Very important - flush() otherwise nothing will happen and the reader thread will get stuck in an infinite
    shell_input.flush()

    shout = []
    exit_status = 0
    for line in shell_stdout:
        if finish in line:
            try:
                mark, code = line.lstrip().rstrip().split(' ')
                exit_status = int(code)
                assert mark == finish
            except:
                shout.append(line)
                if propagate_stdout:
                    shell_output.write(line)
                    shell_stdout.flush()
            else:
                break
        else:
            shout.append(line)
            if propagate_stdout:
                shell_output.write(line)
                shell_stdout.flush()
    return exit_status == 0 if not return_stdout else (exit_status == 0, shout)


def _read_output(shell, log_file, show_output):
    """
    Reads the output of until __END_MARK__ is reached, terminates by returning the exit code
    :param shell: (subprocess.PIPE or StringStream) the shell output handle
    :param log_file: (file) the filehandle to the log file
    :param show_output: (bool) wether to print the output to sys.stdout
    :return: (int) the exitcode
    """
    _, _, _, f = shell
    line = f.readline()
    log_file.write(line)
    if show_output:
        print(line, end='')
    while __END_MARK__ not in line.strip():
        line = f.readline()
        if line:
            log_file.write(line)
            if show_output:
                print(line, end='')
        # If everything worked old_line should contain the exit status
    line, exitcode = line.strip().split(' ')
    log_file.flush()
    try:
        exitcode = int(exitcode)
    except ValueError:
        log_file.close()
        return None
    else:
        log_file.close()
        return exitcode


def _run_vasp_internal(directory=None, cpus=2, show_output=True, return_stdout=False, application=None, hostname=None, partition=None):
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
    shell_handle = Popen(shlex.split('/bin/bash'), stdin=PIPE, stdout=PIPE)
    shell_stdin = TextIOWrapper(shell_handle.stdin, encoding='utf-8')
    shell_stdout = TextIOWrapper(shell_handle.stdout, encoding='utf-8')
    shell_output = StringStream()
    shell = (shell_handle, shell_stdin, shell_stdout, shell_output)

    if directory is None:
        directory = working_directory(getcwd())

    if isinstance(directory, working_directory):
        pass
    elif isinstance(directory, str):
        directory = working_directory(directory)

    with directory:
        log_file = '{}.vasp.log'.format(directory.name)
        log_file_fd = open(log_file, 'w')

        thr_read_log = ThreadWithReturnValue(target=_read_output, args=(shell, log_file_fd, show_output))
        thr_read_log.start()
        preamble, command, binary = get_vasp_configuration(application=application, hostname=hostname, partition=partition)
        command = command.format(cores=cpus, binary=binary)
        echo_cmd = 'echo "{} $?"'.format(__END_MARK__)
        change_command = 'cd {}'.format(getcwd())
        _send_command(shell, change_command, return_stdout=False, propagate_stdout=False)
        for preamble_ in preamble:
            _send_command(shell, preamble_, return_stdout=False, propagate_stdout=False)

        output = _send_command(shell, command, return_stdout=return_stdout, propagate_stdout=True)
        _send_command(shell, echo_cmd, return_stdout=False, propagate_stdout=True)
        thr_read_log.join()

        if return_stdout:
            exitcode, output = output
        else:
            exitcode = output

        if not exitcode:
            logger.warning('VASP did not execute successfully! Exit-Code: "{}"'.format(exitcode))
            return

        return exitcode if not return_stdout else (exitcode, output)
        # Write input files

def vasp(directory=None, cpus=2, show_output=True, return_stdout=False, application=None, hostname=None, partition=None):
    """
    Run VASP in a given directory. The maximum number of cores is limited to six. If no directory is specified VASP
    will be executed on os.getcwd() directory
    :param directory: (str or working_directory) the directory where to execute vasp
    :param cpus: (int) the number of core used to run the VASP subprocess (default: 2)
    :param show_output: (bool) wether to print the output on sys.stdout, passed on to _read_output (default: True)
    :param return_stdout: (bool) wether to return the output as a list of strings. Passed on to _send_command (default: False)
    :return: (bool) or (bool, list of str) exitcode== 0 and output depending on the setting of propagate_stdout
    """
    if cpus < 0:
        cpus = 1
    elif cpus > 6:
        cpus = 6
        logging.getLogger('VASPRunner').warning('Since you have to share our nodew with you colleagues we limited the maximum number of core to six!\n'
                                                'Weil Ihr euch die Systemrssourcen teilen müsst, haben wir die Anzahl der maximalen Kerne pro Job auf sechs limitiert!')
    else:
        cpus = cpus
    return _run_vasp_internal(directory=directory, cpus=cpus, show_output=show_output, return_stdout=return_stdout,
                              application=application, hostname=hostname, partition=partition)

