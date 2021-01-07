import os
import logging
import shlex
import sys
from uuid import uuid4
from os import chdir, getcwd, mkdir
from os.path import join, exists
from shutil import rmtree
from subprocess import PIPE, Popen
from wmaee.core.event import Event, EventHandler
from wmaee.utils import collection, unpack_single
from time import sleep, time as current_time
from io import TextIOWrapper, StringIO
from threading import Thread, Event as ThreadingEvent
from typing import Union, TextIO, Collection, Optional, Tuple, List, Dict, NoReturn, Iterator

try:
    from queue import Queue, Empty
except ImportError:
    from Queue import Queue, Empty

ShellStreams = Union[None, TextIO, Collection[TextIO]]
ShellStream = Union[None, TextIO]


def remove_white(string: str) -> str:
    """
    Removes all whitespaces in a given string
    :param string: (str) the string
    :return: (str) a copy without whitespaces
    """
    whitespace = [' ', '\t', '\n']
    mystr = str(string)
    for removal in whitespace:
        mystr = mystr.replace(removal, '')
    return mystr


class ThreadWithReturnValue(Thread):
    """
    Small wraper around threading.Thread which stores the return value of the executed function
    """

    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}):
        Thread.__init__(self, group, target, name, args, kwargs)
        self._return = None

    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args,
                                        **self._kwargs)

    def join(self, *args, **kwargs):
        Thread.join(self, *args, **kwargs)
        return self._return


class working_directory(object):
    """
    A convenience class which syntactic sugar, allowing the user to change the directories.
    Can also be nested.
    """

    def __init__(self, name=None, prefix=None, delete=False):
        """
        Constructs a working_directory object
        :param name: (str) name of the directory if None is given os.getcwd() will be used (default: None)
        :param prefix: (str) a prefix where to locate the directory (default: None)
        :param delete: (bool) wether to delete the directory after a with clause (default: False)
        """
        self._name = str(uuid4()) if not name else name
        self._delete = delete
        self._curr_dir = getcwd()
        self._active = False
        if prefix is not None:
            self._name = join(prefix, self._name)

    def __enter__(self):
        if not exists(self._name):
            mkdir(self._name)
        chdir(self._name)
        self._active = True

    def __exit__(self, exc_type, exc_val, exc_tb):
        chdir(self._curr_dir)
        if self._delete:
            rmtree(self._name)
        self._active = False

    @property
    def name(self):
        return self._name

    @property
    def active(self):
        return self._active


class LoggerMixin(object):
    """
    A mixin for logger
    """

    def __init__(self, **kwargs):
        super(LoggerMixin, self).__init__(**kwargs)

    @property
    def logger(self):
        return logging.getLogger(self.fullname())

    @classmethod
    def fullname(cls):
        """
        Returns the fully qualified name string of the cls
        :return: (str) the class identifier
        """
        name = '.'.join([
            cls.__module__,
            cls.__name__
        ])
        return name


def is_ipython():
    try:
        from IPython import get_ipython
    except ImportError:
        return False

    return get_ipython() is not None


if is_ipython():
    from tqdm import tqdm_notebook
    tqdm = tqdm_notebook
else:
    from tqdm import tqdm as tqdm_
    tqdm = tqdm_


def get_configuration_directory():
    """
    Build the path the to configuration directory for this module
    :return: (str) the absolute path of the configuration directory
    """
    if 'WMAEE_CONFIG_DIR' not in os.environ:
        if exists('.config'):
            return join(os.getcwd(), '.config')
        else:
            raise RuntimeError('No configuration directory found')
    else:
        return os.environ['WMAEE_CONFIG_DIR']


class OutputEmitter(LoggerMixin):

    def __init__(self, trigger=''):
        super(OutputEmitter, self).__init__()
        self._trigger = trigger
        self.trigger = Event()

    def write(self, line):
        if self._trigger in line:
            self.trigger.fire()

class Shell(LoggerMixin):
    __instance = None

    def __init__(self, restart: Optional[bool] = True, stdout: Optional[ShellStreams] = sys.stdout,
                 stderr: Optional[ShellStreams] = sys.stderr,
                 stdin: Optional[ShellStreams] = None, out_log: Optional[Union[None, TextIO, str]] = None,
                 err_log: Optional[Union[None, TextIO, str]] = None, time: Optional[float] = 0.0001,
                 timing: Optional[float] = True, shell_cmd: Optional[str] = '/bin/bash'):
        """
        Constructs a bash/sh/zsh shell instance
        :param restart: (bool) flag wether to automatically restart the shell subprocess if it died [default: True]
        :param stdout: (TextIO or None) where to redirect the shells stdout stream [default: sys.stdout]
        :param stderr: (TextIO or None) where to redirect the shells stderr stream [default: sys.stderr]
        :param stdin: (TextIO or None) stdin for the shell [default: None] TODO: maybe we should remove this option here
        :param out_log: (TextIO, str or None) a log file where to write the shell stdout [default: None]
        :param err_log: (TextIO, str or None) a log file where to write the shell stderr [default: None]
        :param time: (float) time interval to check the streams for output, prevent core blocking [default: 0.0001]
        :param timing: (bool) wether to time the commands or not [default: True]
        :param shell_cmd: (str) the command of the passed on the to the Popen to create the subprocess [default: /bin/bash]
        """
        super(Shell, self).__init__()
        self._shell_cmd: str = shell_cmd
        self._shell_handle: Popen = Popen(shlex.split(self._shell_cmd), stdin=PIPE, stdout=PIPE, stderr=PIPE)
        self._shell_stdin: TextIO = TextIOWrapper(self._shell_handle.stdin, encoding='utf-8')
        self._shell_stdout: TextIO = TextIOWrapper(self._shell_handle.stdout, encoding='utf-8')
        self._shell_stderr: TextIO = TextIOWrapper(self._shell_handle.stderr, encoding='utf-8')
        self._shell: Tuple[Popen, TextIO, TextIO, TextIO] = (
            self._shell_handle, self._shell_stdin, self._shell_stdout, self._shell_stderr)
        self._restart: bool = restart
        self._time: float = time
        self._timing: float = timing
        self._command_queue: Queue = Queue()
        self.command_finished: Event = Event()
        self.command_started: Event = Event()
        self._out_log_fd: TextIO = open('%s.err.log' % out_log, 'w') if isinstance(out_log, str) else out_log
        self._err_log_fd: TextIO = open('%s.err.log' % err_log, 'w') if isinstance(err_log, str) else err_log
        self._history: Collection[Tuple] = []
        if stdin is None:
            self._default_stdin = StringIO()
            stdin = self._default_stdin
        self._stdin: Union[TextIO, None] = stdin
        self._stderr: Union[TextIO, None] = stderr
        self._stdout: Union[TextIO, None] = stdout
        self._close: ThreadingEvent = ThreadingEvent()
        self._cmd_active: ThreadingEvent = ThreadingEvent()
        self._cmd_data: Union[None, Tuple] = None
        streams = (self._shell_stdout, self._shell_stderr)
        self._queues: Collection[Queue] = [Queue() for _ in streams]
        self._threads: List[Thread] = []
        self._finish_mark: str = '__local_command_finish_mark__'
        for queue, stream in zip(self._queues, streams):
            thread = Thread(target=self._enqueue_output, args=(stream, queue))
            thread.daemon = True
            self._threads.append(thread)
            thread.start()

        self._output_hook: Event = Event()
        self._error_hook: Event = Event()
        self._input_hook: Event = Event()
        self._started_handler: EventHandler = EventHandler('run_block_started_handler',
                                                           self._set_current_command_active)
        self._finished_handler: EventHandler = EventHandler('run_block_finished_handler',
                                                            self._set_current_command_finished)
        self.command_started.set_event_handler(self._started_handler)
        self.command_finished.set_event_handler(self._finished_handler)
        self._remaining_commands: List[Tuple] = []  # local buffer for batch execution
        self._output_hooks: Dict[int, List[EventHandler]] = {}
        self._error_hooks: Dict[int, List[EventHandler]] = {}
        self.output_streams: List[TextIO] = []
        self.error_streams: List[TextIO] = []
        if self._stdout is not None:
            self.output_streams.append(self._stdout)
        if self._out_log_fd is not None:
            self.output_streams.append(self._out_log_fd)
        if self._stderr is not None:
            self.error_streams.append(self._stderr)
        if self._err_log_fd is not None:
            self.error_streams.append(self._err_log_fd)
        self.input_streams: List[TextIO] = [self._shell_stdin]

        self._thr_dist: Thread = Thread(target=self._distribute)
        self._thr_dist.start()

    def _enqueue_output(self, out: TextIO, queue: Queue) -> NoReturn:
        """
        Function to wrap streams into queues to allow for non-blocking I/O
        :param out: (TextIO) the stream used for non-blocking I/O
        :param queue: (Queue) the queue where the stream will be forwarded to
        """
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
        """
        Propagate and forward the input streams
        """
        outq, errq = self._queues
        while True:
            if self._command_queue.qsize() == 0:
                # command are in the queue block until we have worked on them
                if self._close.is_set():
                    break
            for q, s, hk, hks in zip((outq, errq), (self.output_streams, self.error_streams),
                                     (self._output_hook, self._error_hook), (self._output_hooks, self._error_hooks)):
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
                            # if a command is active we write it to the hooks
                            # if self._cmd_data for None
                            if self._cmd_data is not None:
                                if self._timing:
                                    cmd_id, _, _ = self._cmd_data
                                else:
                                    cmd_id, _ = self._cmd_data
                            else:
                                cmd_id = None

                            if cmd_id is not None:
                                if cmd_id in hks:
                                    for hndlr in hks[cmd_id]:
                                        hk.fire_handler(hndlr.name, line)
                    else:
                        # try the next queue
                        sleep(self._time)

    @classmethod
    def get(cls):
        """
        Getter for the global static Shell instance
        :return: (Shell) the global shell instance
        """
        if cls.__instance is None:
            cls.__instance = Shell()
        return cls.__instance

    def restart(self) -> NoReturn:
        """
        Restarts the shell, calls the constructor again, if the current shell is not alive
        """
        if not self.alive:
            self.close()
            # now restart the shell, and initialize everything
            self.__init__(restart=self._restart,
                          stdout=self._stdout,
                          stderr=self._stderr,
                          stdin=self._stdin,
                          out_log=self._out_log_fd,
                          err_log=self._err_log_fd,
                          time=self._time,
                          timing=self._timing,
                          shell_cmd=self._shell_cmd)

    def _send_command(self, cmd: str, raw: Optional[bool] = False) -> Union[None, int]:
        """
        Sends a command to the input stream, and generates the meta-data for the execution
        :param cmd: (str) the command to execute
        :param raw: (bool) if False the command will be executed silently and does not show up in the history
        """
        if self._restart:
            self.restart()
        if not raw:
            # self._shell_stdin.write(cmd + '\n')
            # we have to set the threading event and tell the distributor thread
            # that we have started a command and prevent it from joining
            cmd_id = self._get_next_cmd_id()
            if self._timing:
                current_timing = current_time()
                cmd_data = (cmd_id, cmd, current_timing)
            else:
                cmd_data = (cmd_id, cmd)
            self._command_queue.put(cmd_data)
            echo_cmd = 'echo {} $?\n'.format(self._finish_mark)
            cmd = '; '.join([cmd, echo_cmd])

        self._shell_stdin.write(cmd)
        if not raw:
            # we want to execute this line after we have written the command to the stdin
            self.command_started.fire(*cmd_data)
        self._shell_stdin.flush()
        return cmd_id if not raw else None

    def _set_current_command_active(self, *args):
        if not self._cmd_active.is_set():
            self._cmd_active.set()
            self._cmd_data = args
        else:
            pass
            self._remaining_commands.append(args)
            # raise RuntimeError('It is not possible to run shell commands on parallel')

    def _set_current_command_finished(self, *args) -> NoReturn:
        """
        Internal handler function for the command_finished Event
        :param args: (tuple) represents the command data
        """
        cmd_id, cmd = args[:2]
        if not self._cmd_active.is_set():
            raise RuntimeError('No command is currently active')
        else:
            if self._timing:
                active_id, active_cmd, _ = self._cmd_data
            else:
                active_id, active_cmd = self._cmd_data
            if cmd_id == active_id and active_cmd == cmd:
                # everything is alright
                self._cmd_active.clear()
                # Now remove the event handler of this command
                if cmd_id in self._output_hooks:
                    for handlr in self._output_hooks[cmd_id]:
                        self._output_hook.remove_event_handler(handlr)
                    del self._output_hooks[cmd_id]
                if cmd_id in self._error_hooks:
                    for handlr in self._error_hooks[cmd_id]:
                        self._error_hook.remove_event_handler(handlr)
                    del self._error_hooks[cmd_id]
                if len(self._remaining_commands) > 0:
                    # set the next command active
                    next_cmd_data = self._remaining_commands.pop(0)
                    self._set_current_command_active(*next_cmd_data)
            else:
                self.logger.warning('The command %i:"%s" finished although no listener was registered' % (cmd_id, cmd))

    def _block_until_command_finished(self) -> NoReturn:
        """
        Simple blocker function, which waits for a threading.Event to be set
        """
        while self._cmd_active.is_set():
            sleep(self._time)

    def _shell_alive(self) -> bool:
        """
        Thest wether a shell is still open or not
        """
        return self._shell_handle and not isinstance(self._shell_handle.poll(), int)

    def _get_next_cmd_id(self) -> int:
        """
        Computes the next cmd id fot the nex command
        :return: (int) the command id
        """
        return len(self._history) + self._command_queue.qsize()

    def run(self, cmd: Union[str, Collection[str]], block: Optional[bool] = True,
            output: Optional[ShellStream] = None,
            error: Optional[ShellStream] = None,
            return_out: Optional[bool] = False, return_err: Optional[bool] = False,
            return_id: Optional[bool] = True, return_exit: Optional[bool] = False) -> Union[tuple, Collection[tuple]]:
        """
        Execute a command in the current shell
        :param cmd: (str, list of str) a command or a list of commands
        :param block: (bool) wether to block execution until the command(s) is/are finished [default: True]
        :param output: (None, TextIO or list of TextIO) a stream of a list of stream where stdout should be forked to [default: None]
        :param error: (None, TextIO or list of TextIO) a stream of a list of stream where stderr should be forked to [default: None]
        :param return_out: (bool) a flag if stdout should be returned can only be used with block=True [default: False]
        :param return_err: (bool) a flag if stderr should be returned can only be used with block=True [default: False]
        :param return_id: (bool) a flag if command ids should be returned [default: True]
        :param return_exit: (bool) a flag if exit codes be returned can only be used with block=True [default: False]
        :return: (tuple, list of tuple) the fields requested as given by return_out, return_err, return_id and
        return_exit
        """
        propagator = lambda stream: lambda line: stream.write(line)
        cmd = collection(cmd)
        err_buffers = [StringIO() for _ in cmd] if return_err else []
        out_buffers = [StringIO() for _ in cmd] if return_out else []
        cmd_ids = []
        exit_codes = []

        def exit_code_extractor(cmd_id):
            def _handler(other_id, *args):
                if cmd_id == other_id:
                    exit_codes.append(args[-1])

            return _handler

        if any((return_out, return_err, return_exit)) and not block:
            block = True
            self.logger.warning(
                'The "return_exit", "return_err" and "return_out" keyword argument cannot be used in combination with block=False. '
                'The settings will be overridden!')

        for i, command in enumerate(cmd):
            next_cmd_id = self._get_next_cmd_id()
            if output is not None:
                output_handlers = [EventHandler('output_hook_%i_%i' % (next_cmd_id, j), propagator(o)) for j, o in
                                   enumerate(collection(output))]
                self._output_hooks[next_cmd_id] = output_handlers
                for handler in output_handlers:
                    self._output_hook.add_event_handler(handler)
            if error is not None:
                error_handlers = [EventHandler('error_hook_%i_%i' % (next_cmd_id, j), propagator(o)) for j, o in
                                  enumerate(collection(error))]
                self._error_hooks[next_cmd_id] = error_handlers
                for handler in error_handlers:
                    self._error_hook.add_event_handler(handler)
            if return_out:
                # also this guys will be removed at finished call, thus we keep everything clean
                return_out_handler = EventHandler('return_output_hook_%i' % next_cmd_id, propagator(out_buffers[i]))
                if next_cmd_id not in self._output_hooks:
                    self._output_hooks[next_cmd_id] = [return_out_handler]
                else:
                    self._output_hooks[next_cmd_id].append(return_out_handler)
                self._output_hook.add_event_handler(return_out_handler)
            if return_err:
                return_err_handler = EventHandler('return_error_hook_%i' % next_cmd_id, propagator(err_buffers[i]))
                if next_cmd_id not in self._error_hooks:
                    self._error_hooks[next_cmd_id] = [return_err_handler]
                else:
                    self._error_hooks[next_cmd_id].append(return_err_handler)
                self._error_hook.add_event_handler(return_err_handler)
            if return_exit:
                exit_handler_name = 'run_exit_handler_%i' % next_cmd_id
                self.command_finished.add_event_handler(
                    EventHandler(exit_handler_name, exit_code_extractor(next_cmd_id)))
                exit_code_length = len(exit_codes)
            self._send_command(command)
            if block:
                self._block_until_command_finished()
            if return_exit:
                self.command_finished.remove_event_handler(exit_handler_name)
                assert len(exit_codes) - exit_code_length == 1
            # now exit codes must be longer by exactly one command
            cmd_ids.append(next_cmd_id)

        if return_out:
            for out_buf in out_buffers:
                out_buf.seek(0)
            out_buffers = [out_buf.getvalue() for out_buf in out_buffers]
        if return_err:
            for err_buf in err_buffers:
                err_buf.seek(0)
            err_buffers = [err_buf.getvalue() for err_buf in err_buffers]
        if any((return_err, return_out, return_id)):
            r = [data for data, b in
                 zip((cmd_ids, out_buffers, err_buffers, exit_codes), (return_id, return_out, return_err, return_exit))
                 if b]
            r = list(zip(*r))
            return unpack_single(r)

    def __call__(self, cmd: Union[str, Collection[str]], block: Optional[bool] = True,
                 output: Optional[ShellStream] = None,
                 error: Optional[ShellStream] = None,
                 return_out: Optional[bool] = False, return_err: Optional[bool] = False,
                 return_id: Optional[bool] = True, return_exit: Optional[bool] = False) -> Union[
        tuple, Collection[tuple]]:
        """
        Syntactic sugar for self.run
        :param cmd: (str, list of str) a command or a list of commands
        :param block: (bool) wether to block execution until the command(s) is/are finished [default: True]
        :param output: (None, TextIO or list of TextIO) a stream of a list of stream where stdout should be forked to [default: None]
        :param error: (None, TextIO or list of TextIO) a stream of a list of stream where stderr should be forked to [default: None]
        :param return_out: (bool) a flag if stdout should be returned can only be used with block=True [default: False]
        :param return_err: (bool) a flag if stderr should be returned can only be used with block=True [default: False]
        :param return_id: (bool) a flag if command ids should be returned [default: True]
        :param return_exit: (bool) a flag if exit codes be returned can only be used with block=True [default: False]
        :return: (tuple, list of tuple) the fields requested as given by return_out, return_err, return_id and
        return_exit
        """
        return self.run(cmd, block=block, output=output, error=error, return_out=return_out, return_err=return_err,
                        return_id=return_id, return_exit=return_id)

    @property
    def alive(self) -> bool:
        """
        Queries the shell for its health state
        :return: (bool) alive flag
        """
        return self._shell_alive()

    @property
    def history(self) -> Iterator[tuple]:
        """
        Returns the history in reversed order
        :return: (iterator of tuple) the content depend  if timing is set
        """
        return iter(reversed(self._history))

    @staticmethod
    def _close_if_possible(obj):
        if hasattr(obj, "close"):
            obj.close()

    def close(self) -> NoReturn:
        """
        Closes the shell if it is still alive
        """
        if self.alive:
            # we do not want the threads to wait for a response
            self._send_command('exit', raw=True)
            self._close.set()
            self._thr_dist.join()
            for thread in self._threads:
                thread.join(1)
            # terminate the subprocess, after we have stopped the forwarder threads
            self._shell_handle.terminate()
            # now we have to close log files and so on
            things_to_close = [
                self._out_log_fd,
                self._err_log_fd,
                self._shell_stdout,
                self._shell_stdin,
                self._shell_stderr
            ]
            for thing in things_to_close:
                try:
                    self._close_if_possible(thing)
                except RuntimeError:
                    raise

    def __enter__(self, *args):
        """
        fucntion for with semantics, restarts the shell if restart flag is set
        :return: (Shell) the self Shell instance
        """
        if self.alive:
            return self
        else:
            if self._restart:
                self.restart()
            return self

    def __exit__(self, *args) -> NoReturn:
        """
        function for with semantics. Closes the shell
        """
        if self.alive:
            self.close()

