
import os
import logging
import shlex
import json
import sys
from subprocess import PIPE, Popen
from wmaee.core.common import working_directory, get_configuration_directory, LoggerMixin
from wmaee.core.event import Event
from os.path import join
from os import getcwd
from time import sleep, time as current_time
from io import TextIOWrapper, StringIO
from threading import Thread, Event as ThreadingEvent
try:
    from queue import Queue, Empty
except ImportError:
    from Queue import Queue, Empty


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


class Shell(LoggerMixin):

    __instance = None

    def __init__(self, restart=True, stdout=sys.stdout, stderr=sys.stderr, stdin=None, out_log=None, err_log=None, time=0.01, timing=True):
        super(Shell, self).__init__()
        self.shell_handle = Popen(shlex.split('/bin/bash'), stdin=PIPE, stdout=PIPE, stderr=PIPE)
        self.shell_stdin = TextIOWrapper(self.shell_handle.stdin, encoding='utf-8')
        self.shell_stdout = TextIOWrapper(self.shell_handle.stdout, encoding='utf-8')
        self.shell_stderr = TextIOWrapper(self.shell_handle.stderr, encoding='utf-8')
        self.shell = (self.shell_handle, self.shell_stdin, self.shell_stdout, self.shell_stderr)
        self._restart = restart
        self._time = time
        self._timing = timing
        self._command_queue = Queue()
        self.command_finished = Event()
        self.command_started = Event()
        self._out_log_fd = open('%s.err.log' % out_log, 'w') if isinstance(out_log, str) else out_log
        self._err_log_fd = open('%s.err.log' % err_log, 'w') if isinstance(err_log, str) else err_log
        self._history = []
        if stdin is None:
            self._default_stdin = StringIO()
            stdin = self._default_stdin
        self._stdin = stdin
        self._stderr = stderr
        self._stdout = stdout
        self._close = ThreadingEvent()
        streams = (self.shell_stdout, self.shell_stderr)
        self._queues = [Queue() for _ in streams]
        self._threads = []
        self._finish_mark = '__local_command_finish_mark__'
        for queue, stream in zip(self._queues, streams):
            thread = Thread(target=self._enqueue_output, args=(stream, queue))
            thread.daemon = True
            self._threads.append(thread)
            thread.start()

        # TODO: Implement command listening and output buffering in distributor thread
        self._output_hooks = []
        self._input_hooks = []
        self.output_streams = []
        if self._stdout is not None:
            self.output_streams.append(self._stdout)
        if self._out_log_fd is not None:
            self.output_streams.append(self._out_log_fd)
        self.error_streams = []
        if self._stderr is not None:
            self.error_streams.append(self._stderr)
        if self._err_log_fd is not None:
            self.error_streams.append(self._err_log_fd)
        self.input_streams = [self.shell_stdin]

        self._thr_dist = Thread(target=self._distribute)
        self._thr_dist.start()

    def _enqueue_output(self, out, queue):
        self.logger.debug('Forwarder thread started')
        for line in iter(out.readline, ''):
            # keep forwarding if command is active, as long as commands are in the queue
            if self._command_queue.qsize() == 0:
                if self._close.is_set():
                    break
            if line:
                queue.put(line)
        self.logger.debug('Forwarder thread finished')

    def _distribute(self):
        outq, errq = self._queues
        while True:
            if self._command_queue.qsize() == 0:
                # command are in the queue block until we have worked on them
                if self._close.is_set():
                    break
            for q, s in zip((outq, errq), (self.output_streams, self.error_streams)):
                try:
                    line = q.get_nowait()
                except Empty:
                    sleep(self._time)
                else:
                    if line:
                        if self._finish_mark in line:
                            # it is the echo finish mark command -> do not propagate it
                            # we allow the loop to exit
                            # get the timing
                            if self._timing:
                                cmd_id, cmd, start_time = self._command_queue.get()
                            else:
                                cmd_id, cmd = self._command_queue.get()
                            # try to fetch the exit code
                            try:
                                mark, exit_code = line.split(' ')
                                exit_code = int(exit_code)
                            except:
                                exit_code = float('nan')
                            if self._timing:
                                elapsed = current_time() - start_time
                                result = (cmd_id, cmd, elapsed, exit_code)
                            else:
                                result = (cmd_id, cmd, exit_code)
                            self._history.append(result)
                            self.command_finished.fire(*result)
                        else:
                            for stream in s:
                                stream.write(line)
                    else:
                        # try the next queue
                        sleep(self._time)

    @classmethod
    def get(cls):
        if cls.__instance is None:
            cls.__instance = Shell()
        return cls.__instance

    def restart(self):
        if not self.alive:
            self.close()
            self.__init__(restart=self._restart, stdout=self._stdout, stderr=self._stderr, stdin=self._stdin, out_log=self._out_log_fd, err_log=self._err_log_fd)

    def _send_command(self, cmd, raw=False):
        if self._restart:
            self.restart()
        if cmd.endswith('\n'):
            cmd = cmd.rsplit()
        if not raw:
            # self._shell_stdin.write(cmd + '\n')
            # we have to set the threading event and tell the distributor thread
            # that we have started a command and prevent it from joining
            cmd_id = len(self._history) + self._command_queue.qsize()
            if self._timing:
                current_timing = current_time()
                self._command_queue.put((cmd_id, cmd, current_timing))
            else:
                self._command_queue.put((cmd_id, cmd))
            echo_cmd = 'echo {} $?\n'.format(self._finish_mark)
            cmd = '; '.join([cmd, echo_cmd])

        self.shell_stdin.write(cmd)
        self.command_started.fire()
        self.shell_stdin.flush()

    def run(self, cmd, output):
        if isinstance(cmd, str):
            cmd = [cmd]
        for i, command in enumerate(cmd):
            self._send_command(command)

    def __call__(self, *args, **kwargs):
        return self.run(*args, **kwargs)

    @property
    def alive(self):
        return _shell_alive(self.shell_handle)

    @property
    def history(self):
        return iter(reversed(self._history))

    def close(self):
        if self.alive:
            # we do not want the threads to wait for a response
            self._send_command('exit', raw=True)
            self._close.set()
            self._thr_dist.join()

    def __enter__(self, *args):
        if self.alive:
            return self
        else:
            if self._restart:
                self.restart()
            return self

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

        with Shell(stdout=sys.stdout if show_output else None) as shell:
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

