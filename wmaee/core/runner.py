
import os
import logging
import shlex
import json
from subprocess import PIPE, Popen
from wmaee.core.common import working_directory, ThreadWithReturnValue, get_configuration_directory
from wmaee.core.io import StringStream
from os.path import join
from os import getcwd
from time import sleep
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


def _read_output(shell, log_file, show_output, time):
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
                if __END_MARK__ not in line:
                    print(line, end='')
            # sleep here otherwise the thread will use 100% cpu power
        sleep(time)
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


class Shell:
    __instance = None

    def __init__(self, restart=True, show_output=True, time=0.01, propagate_output=True):
        self.shell_handle = Popen(shlex.split('/bin/bash'), stdin=PIPE, stdout=PIPE)
        self.shell_stdin = TextIOWrapper(self.shell_handle.stdin, encoding='utf-8')
        self.shell_stdout = TextIOWrapper(self.shell_handle.stdout, encoding='utf-8')
        self.shell_output = StringStream()
        self.shell = (self.shell_handle, self.shell_stdin, self.shell_stdout, self.shell_output)
        self.restart = restart
        self.time = time
        self.log_file_fd = open('shell.log', 'w')
        self.thr_read_log = ThreadWithReturnValue(target=_read_output,
                                                  args=(self.shell, self.log_file_fd, show_output, self.time))
        self.thr_read_log.start()
        self.show_output = show_output
        self.propagate_output = propagate_output

    @classmethod
    def get(cls):
        if cls.__instance is None:
            cls.__instance = Shell()
        return cls.__instance

    def run(self, cmd, return_stdout=False):
        if isinstance(cmd, str):
            cmd = [cmd]
        output = []
        if not _shell_alive(self.shell_handle):
            self.__init__(restart=self.restart, show_output=self.show_output, time=self.time)
        for i, command in enumerate(cmd):
            o = _send_command(self.shell, command, return_stdout=return_stdout, propagate_stdout=self.propagate_output)
            output.append(o)
        return output

    def __call__(self, *args):
        return self.run(*args)

    @property
    def alive(self):
        return _shell_alive(self.shell_handle)

    def close(self):
        if self.alive:
            self.shell_output.write('%s 0' % __END_MARK__)
            self.thr_read_log.join()
            self.shell_stdin.write('exit\n')

    def __enter__(self, *args):
        if self.alive:
            return self
        else:
            if self.restart:
                self.__init__(restart=self.restart, show_output=self.show_output, time=self.time)
            else:
                raise RuntimeError('Shell is not alive')

    def __exit__(self, *args):
        if self.alive:
            self.close()


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


    if directory is None:
        directory = working_directory(getcwd())

    if isinstance(directory, working_directory):
        pass
    elif isinstance(directory, str):
        directory = working_directory(directory)

    with directory:
        preamble, command, binary = get_vasp_configuration(application=application, hostname=hostname, partition=partition)
        command = command.format(cores=cpus, binary=binary)
        echo_cmd = 'echo "{} $?"'.format(__END_MARK__)
        change_command = 'cd {}'.format(getcwd())

        with Shell(show_output=show_output) as shell:
            shell.show_output = False
            shell.propagate_output = False
            shell.run(change_command, return_stdout=False)
            shell.propagate_output = True
            shell.show_output = show_output
            for preamble_ in preamble:
                shell.run(preamble_, return_stdout=False)
            output = shell.run(command, return_stdout=return_stdout)

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
    else:
        cpus = int(cpus)
    return _run_vasp_internal(directory=directory, cpus=cpus, show_output=show_output, return_stdout=return_stdout,
                              application=application, hostname=hostname, partition=partition)

