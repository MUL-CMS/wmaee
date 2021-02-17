import re
import sys
import uuid
import shlex
import asyncio
import datetime
import collections.abc
from wmaee.core.aio import create_standard_streams
from typing import Union, Optional, Tuple, Callable, NoReturn, Iterable, ByteString, AnyStr, IO, Dict

# create shortcut types
LineProcessor = Callable[[bytes, Iterable[asyncio.Event]], NoReturn]
OutputStream = Union[IO, asyncio.StreamWriter]
InputStream = Union[IO, asyncio.StreamReader]
OptionalEventLoop = Optional[Union[asyncio.AbstractEventLoop, None]]

# namedtuples for more consize code. Always in order (stdin, stdout, stderr)
Callbacks = collections.namedtuple("Callback", ("stdin", "stdout", "stderr"), defaults=([], [], []))
Forks = collections.namedtuple("Forks", ("stdin", "stdout", "stderr"), defaults=([], [], []))

# shortcut for brevity
PIPE = asyncio.subprocess.PIPE

START_MARK = "__start__"
FINISH_MARK = "__finish__"


class Shell:

    def __init__(self, stdin=None, stdout=None, stderr=None, callbacks=None, forks=None, shell_cmd="/bin/bash"):
        """
        Create a shell object which holds the subprocess handle
        :param handle: (asyncio.subprocess.Process) the internal shell handle
        :param stdin: (input stream) the input stream which is forwarded to handle.stdin
        :param stdout: (output stream) the output stream where handle.stdout is forwarded to
        :param stderr: (output stream) the output stream where handle.stderr is forwarded to
        :param callbacks: (Callbacks or tuple of callable) callbacks for the streams in order (in, out err)
        """
        self._handle = None
        self._cmd = shell_cmd
        self._forwarders = None
        self._stdout = stdout
        self._stderr = stderr
        self._stdin = stdin
        self._forks = forks or Forks()
        self._callbacks = callbacks or Callbacks()
        self.commands = {}

    def set_handle(self, handle):
        self._handle = handle

    def set_stdin(self, stdin):
        self._stdin = stdin

    def set_stdout(self, stdout):
        self._stdout = stdout

    def set_stderr(self, stderr):
        self._stderr = stderr

    def set_forwarders(self, coros):
        self._forwarders = coros

    @property
    def handle(self):
        return self._handle

    @property
    def stdout(self):
        return self._stdout

    @property
    def stderr(self):
        return self._stderr

    @property
    def stdin(self):
        return self._stdin

    @property
    def callbacks(self):
        return self._callbacks

    @property
    def forks(self):
        return self._forks

    async def __aenter__(self):
        handle, stdin, stdout, stderr = await create_shell_handle(self._stdin, self._stdout, self._stderr,
                                                                  shell_cmd=self._cmd)
        self.set_handle(handle)
        self.set_stdin(stdin)
        self.set_stdout(stdout)
        self.set_stderr(stderr)
        self.set_forwarders(create_io_forwarding_tasks(self))

        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        self._handle.terminate()
        await self._handle.wait()
        for coro in self._forwarders.values(): coro.cancel()


async def write_message(writer: OutputStream, msg: bytes, flush: Optional[bool] = False) -> NoReturn:
    """
    Helper method which can write to IO type ans Stream typed writer objects
    :param writer: (IO or output stream) where {msg} is written to
    :param msg: (bytes) the data
    :param flush: (bool) wether to flush after each writer or not (default: False)
    """
    writer.write(msg)
    if hasattr(writer, "drain"):
        await writer.drain()
    elif flush and hasattr(writer, "flush"):
        writer.flush()


def make_callback_forker(get_forks: Callable[[], Iterable[OutputStream]],
                         flush: Optional[bool] = False) -> LineProcessor:
    """
    Factory function for a line processor which forks an incoming message to a list of output streams
    :param forks: (callable -> iterable output stream) output writer tuple where the content of the reader will be forked to
    :param flush: (bool) wether to flush the writer after each message (default: False)
    :return: (line processor) the processing function
    """

    async def _line_processor(msg):
        for fork in get_forks():
            await write_message(fork, msg, flush=flush)

    return _line_processor


# TODO: Figure out if we really need this one
def make_callback_callbacks(get_callbacks: Callable[[], Iterable[Callable[[bytes], NoReturn]]]) -> LineProcessor:
    """
    Factory function for a line processor which calls a callback function whenever a message is received
    :param callback: (iterable of callable) logging callbacks execute for each piped line
    :return: (line processor) the processing function
    """

    def _line_processor(msg):
        for cb in get_callbacks():
            cb(msg)

    return _line_processor


async def processor_pipe(reader: InputStream, writer: OutputStream,
                         callbacks: Optional[Union[Iterable[LineProcessor], None]] = None,
                         flush: Optional[bool] = False):
    """
    a pipe which intercepts and looks at the data, to determine when a shell command was started and to determine its
    end time respectively
    :param reader: (input stream) the input reader
    :param writer: (output stream) output end of the pipe
    :param processors: (iterable of line processors) line processor functions (default=None)
    :param flush: (bool) wether to flush the writer after each message (default: False)
    """
    # create the regex patterns which will mach, start and end of the executed command
    processors = callbacks or []  # loop over empty list

    # we do not know if a line processor is an ordinary function or a coro, thus we handle both cases

    async def execute_line_processor(func, *args, **kwargs):
        return await func(*args, **kwargs) if asyncio.iscoroutinefunction(func) else func(*args, **kwargs)

    while not reader.at_eof():
        msg = await reader.readline()
        # execute the line processors
        for lp in processors:
            await execute_line_processor(lp, msg)

        # finally forward the msg to the writer
        await write_message(writer, msg, flush=flush)


async def create_shell_handle(stdin: Optional[Union[InputStream, None]] = None,
                              stdout: Optional[Union[OutputStream, None]] = None,
                              stderr: Optional[Union[OutputStream, None]] = None, *,
                              shell_cmd: Optional[AnyStr] = "/bin/bash",
                              loop: OptionalEventLoop = None) -> Tuple[
    asyncio.subprocess.Process, InputStream, OutputStream, OutputStream]:
    """
    Creates an asyncio subprocess handle and connects it with the streams specified in the arguments
    :param stdin: (input stream) stdin stream
    :param stdout: (output stream) stdout stream
    :param stderr: (output stream) stderr stream
    :param shell_cmd: (str) the command used to open the handle (default: /bin/bash)
    :param loop: (asyncio.EventLoop) the asnycio.EventLoop (default: None)
    :return: (Shell) the shell object which wraps streams, forks and callbacks
    """
    loop = asyncio.get_event_loop() if loop is None else loop
    # determine which streams need protection

    # wrap the standard IOs into async
    sstdin, sstdout, sstderr = await create_standard_streams(sys.stdin, sys.stdout, sys.stderr)

    shell_handle = await asyncio.create_subprocess_exec(*shlex.split(shell_cmd), stdout=PIPE, stderr=PIPE, stdin=PIPE,
                                                        loop=loop)
    # warp everything needed into an Shell object
    return shell_handle, stdin or sstdin, stdout or sstdout, stderr or sstderr


# For brevity we define a ForkType here
ForksType = Union[Tuple[Iterable[OutputStream], Iterable[OutputStream], Iterable[OutputStream]], Forks]
CallbacksType = Union[Tuple[Iterable[Callable[[bytes], NoReturn]], Iterable[Callable[[bytes], NoReturn]], Iterable[
    Callable[[bytes], NoReturn]]], Callbacks]


def create_io_forwarding_tasks(shell: Shell) -> Dict[
    AnyStr, asyncio.Task]:
    """
    We follow the stategy here that pipes are only open as long as commands are running. Thus when a command is to be
    executed we create the pipe tasks
    :param shell: (Shell) the shell object encapsulating the pipe streams
    :return: (dict of str and asnycio.Task) the pipe coros
    """

    # retrieve forks and possibly override it with the args passed in
    # create a line processor factory
    make_callback_list = lambda get_forks, get_cb: [
        #        mak(shell, encoding=encoding),
        make_callback_forker(get_forks, flush=False),
        make_callback_callbacks(get_cb)
    ]
    fork_getter_factory = lambda num: lambda: shell.forks[num]
    cb_getter_factory = lambda num: lambda: shell.callbacks[num]
    forwarding_tasks = dict(
        stdin=asyncio.create_task(
            processor_pipe(
                shell.stdin,
                shell.handle.stdin,
                callbacks=make_callback_list(fork_getter_factory(0), cb_getter_factory(0))
            )
        ),
        stdout=asyncio.create_task(
            processor_pipe(
                shell.handle.stdout,
                shell.stdout,
                callbacks=make_callback_list(fork_getter_factory(1), cb_getter_factory(1))
            )
        ),
        stderr=asyncio.create_task(
            processor_pipe(
                shell.handle.stderr,
                shell.stderr,
                callbacks=make_callback_list(fork_getter_factory(2), cb_getter_factory(2))
            )
        )
    )
    return forwarding_tasks


async def main(loop: OptionalEventLoop = None, shell_cmd: Optional[AnyStr] = "/bin/bash") -> NoReturn:
    """
    Runs a
    :param loop: (asyncio.EventLoop) the event loop
    :param shell_cmd: (str) the command used to open the handle (default: /bin/bash)
    """
    loop = loop or asyncio.get_event_loop()
    async with Shell(shell_cmd=shell_cmd) as shell:
        shell.handle.stdin.write(b"ls -lsfag\n")
        await shell.handle.stdin.drain()
        await shell.handle.wait()
    # await send_command(shell, shell_cmd)


# Create a clause if this script is running as main
if __name__ == '__main__':
    asyncio.get_event_loop().run_until_complete(main())
