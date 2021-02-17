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

# helper function creating a regex which matches the start mark of a command for a given uuid
START_PATTERN = re.compile(rf"^{START_MARK}\s*(?P<identifier>[a-z0-9]{32})$")
# helper function creating a regex which matches the end mark of a command for a given uuid
FINISH_PATTERN = re.compile(
    rf"^{FINISH_MARK}\s*(?P<identifier>[a-z0-9]{32})\s*(?P<returncode>\d+)$")


def make_command_marks() -> Tuple[AnyStr, AnyStr, AnyStr]:
    """
    Helper function to create start and finish mark:
    :return: (str, str, str) the start and finish mark with a unique id
    """
    identifier = uuid.uuid4().hex
    return identifier, f"{START_MARK} {identifier}", f"{FINISH_MARK} {identifier}"


def wrap_command(command: AnyStr, encoding: Optional[AnyStr] = "utf-8") -> Tuple[bytes, bytes]:
    """
    Builds an echo start_mark and echo finish_mark around {command}
    :param command: (str) the command to execute
    :param encoding: (str) the charset encoding (default: utf-8)
    :return: (bytes, bytes) the encoded identifier and the command envelope
    """
    identifier, sm, fm = make_command_marks()
    return identifier.encode(encoding=encoding), f"""echo "{sm}" && {command} ; echo "{fm}" $?

    """.encode(encoding=encoding)


class Shell:

    def __init__(self, handle, stdin, stdout, stderr, callbacks=None, forks=None):
        """
        Create a shell object which holds the subprocess handle
        :param handle: (asyncio.subprocess.Process) the internal shell handle
        :param stdin: (input stream) the input stream which is forwarded to handle.stdin
        :param stdout: (output stream) the output stream where handle.stdout is forwarded to
        :param stderr: (output stream) the output stream where handle.stderr is forwarded to
        :param callbacks: (Callbacks or tuple of callable) callbacks for the streams in order (in, out err)
        :param forks: (has write method) streams or IO to which the standard streams are foked to (in, out err)
        """
        self._handle = handle
        self._stdout = stdout
        self._stderr = stderr
        self._stdin = stdin
        self._forks = forks or Forks()
        self._callbacks = callbacks or Callbacks()
        self.commands = {}

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
    def forks(self):
        return self._forks

    @property
    def callbacks(self):
        return self._callbacks


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


def make_line_processor_command_wrapper(shell: Shell,  encoding: Optional[AnyStr] = "utf-8") -> LineProcessor:
    """
    Factory function for a line processor which takes care of a function
    :param shell: (Shell) Shell objets to store command execution data, such as e. g. return-codes and timings
    :param identifier: (str) a unique command identifier
    :param encoding: (str) the encoding of the stream (default: utf-8)
    :return: (line processor) the processing function
    """

    def _line_processor(msg, events):
        skip_line, *_ = events
        # we assume byte mode -> decode
        line = msg.decode(encoding=encoding).strip()
        mstart = START_PATTERN.match(line)
        mfinish = FINISH_PATTERN.match(line)
        if mstart:
            # store the info in the shell when we pass the starting line
            cmd_id = mstart.groupdict()["identifier"]
            shell.commands[cmd_id]["started"] = datetime.datetime.now()
            skip_line.set()
        elif mfinish:
            cmd_id = mfinish.groupdict()["identifier"]
            shell.commands[cmd_id]["finished"] = datetime.datetime.now()
            shell.commands[cmd_id]["returncode"] = int(mfinish.groupdict()["returncode"])
            # we do not forward this list but rather exit form here
            skip_line.set()

    return _line_processor


def make_line_processor_forker(forks: Iterable[OutputStream], flush: Optional[bool] = False) -> LineProcessor:
    """
    Factory function for a line processor which forks an incoming message to a list of output streams
    :param forks: (iterable of IO or output stream) output writer tuple where the content of the reader will be forked to
    :param flush: (bool) wether to flush the writer after each message (default: False)
    :return: (line processor) the processing function
    """

    async def _line_processor(msg, _):
        for fork in forks:
            await write_message(fork, msg, flush=flush)

    return _line_processor

# TODO: Figure out if we really need this one
def make_line_processor_callback(callbacks: Iterable[Callable[[bytes], NoReturn]]) -> LineProcessor:
    """
    Factory function for a line processor which calls a callback function whenever a message is received
    :param callback: (iterable of callable) logging callbacks execute for each piped line
    :return: (line processor) the processing function
    """

    def _line_processor(msg, _):
        for cb in callbacks:
            cb(msg)

    return _line_processor


async def processor_pipe(reader: InputStream, writer: OutputStream,
                         processors: Optional[Union[Iterable[LineProcessor], None]] = None,
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
    processors = processors or []  # loop over empty list
    skip_line = asyncio.Event()  # use to simulate continue
    break_loop = asyncio.Event()  # used to simulate break
    events = (skip_line, break_loop)

    # we do not know if a line processor is an ordinary function or a coro, thus we handle both cases

    async def execute_line_processor(func, *args, **kwargs):
        return await func(*args, **kwargs) if asyncio.iscoroutinefunction(func) else func(*args, **kwargs)

    while not reader.at_eof():
        msg = await reader.readline()
        # execute the line processors
        for lp in processors:
            await execute_line_processor(lp, msg, events)

        # perform flow control and reset the signals
        if skip_line.is_set():
            skip_line.clear()
            continue
        if break_loop.is_set():
            break_loop.clear()
            break

        # finally forward the msg to the writer
        await write_message(writer, msg, flush=flush)


async def create_shell_handle(stdin: Optional[Union[InputStream, None]] = None,
                              stdout: Optional[Union[OutputStream, None]] = None,
                              stderr: Optional[Union[OutputStream, None]] = None, *,
                              shell_cmd: Optional[AnyStr] = "/bin/bash",
                              loop: OptionalEventLoop = None) -> Shell:
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
    # TODO: deprecated remove this feature
    protect = tuple(map(lambda x: x is None, (stdin, stdout, stderr)))

    # wrap the standard IOs into async
    sstdin, sstdout, sstderr = await create_standard_streams(sys.stdin, sys.stdout, sys.stderr)

    shell_handle = await asyncio.create_subprocess_exec(*shlex.split(shell_cmd), stdout=PIPE, stderr=PIPE, stdin=PIPE,
                                                        loop=loop)
    # warp everything needed into an Shell object
    return Shell(shell_handle, stdin or sstdin, stdout or sstdout, stderr or sstderr)


# For brevity we define a ForkType here
ForksType = Union[Tuple[Iterable[OutputStream], Iterable[OutputStream], Iterable[OutputStream]], Forks]
CallbacksType = Union[Tuple[Iterable[Callable[[bytes], NoReturn]], Iterable[Callable[[bytes], NoReturn]], Iterable[
    Callable[[bytes], NoReturn]]], Callbacks]


def create_io_forwarding_tasks(shell: Shell, forks: Optional[Union[ForksType, None]] = None,
                               callbacks=None, encoding: Optional[Union[CallbacksType, None]] = "utf-8") -> Dict[
    AnyStr, asyncio.Task]:
    """
    We follow the stategy here that pipes are only open as long as commands are running. Thus when a command is to be
    executed we create the pipe tasks
    :param shell: (Shell) the shell object encapsulating the pipe streams
    :param identifier: (bytes) the command identifier
    :param forks: (tuple of iterable of output stream) stream to fork stdin, stdout, stderr (default: None)
    :param callbacks: (tuple of iterable of callable) cbs sor stdin, stdout, stderr when msg is received (default: None)
    :param encoding: (str) the charset encoding (default: utf-8)
    :return: (dict of str and asnycio.Task) the pipe coros
    """

    # retrieve forks and possibly override it with the args passed in
    stdin_cb, stdout_cb, stderr_cb = shell.callbacks if callbacks is None else callbacks
    stdin_fk, stdout_fk, stderr_fk = shell.forks if forks is None else forks
    print(stdin_fk, stdout_fk, stderr_fk)
    # create a line processor factory
    make_line_processors = lambda forks, cb: [
#        make_line_processor_command_wrapper(shell, encoding=encoding),
        make_line_processor_forker(forks),
        make_line_processor_callback(cb)
    ]
    # create the pipe tasks
    forwarding_tasks = dict(
        stdin=asyncio.create_task(
            processor_pipe(
                shell.stdin,
                shell.handle.stdin,
                processors=make_line_processors(stdin_fk, stdin_cb),
            )
        ),
        stdout=asyncio.create_task(
            processor_pipe(
                shell.handle.stdout,
                shell.stdout,
                processors=make_line_processors(stdout_fk, stdout_cb),
            )
        ),
        stderr=asyncio.create_task(
            processor_pipe(
                shell.handle.stderr,
                shell.stderr,
                processors=make_line_processors(stderr_fk, stderr_cb),
            )
        )
    )
    return forwarding_tasks


async def send_command(shell: Shell, command: AnyStr, forks: Optional[Union[ForksType, None]] = None,
                       callbacks: Optional[Union[CallbacksType, None]] = None) -> NoReturn:
    """
    Creates a tasks and waits until one of the pipe tasks manages to intercept a finish mark
    :param shell: (Shell) the shell where the command will run
    :param command: (str) the command to execute
    :param forks: (tuple of iterable of output stream) stream to fork stdin, stdout, stderr (default: None)
    :param callbacks: (tuple of iterable of callable) cbs sor stdin, stdout, stderr when msg is received (default: None)
    """

    identifier, cmd = wrap_command(command)
    # tell the shell that the command is acutally being executed
    shell.commands[identifier.decode()] = dict(content=command)
    # send it to its stdin
    shell.handle.stdin.write(cmd)
    await shell.handle.stdin.drain()
    # create the I/O pipe coros and wait for one of the to complete (usually stdout). Then cancel the others


async def main(loop: OptionalEventLoop = None, shell_cmd: Optional[AnyStr] = "/bin/bash") -> NoReturn:
    """
    Runs a
    :param loop: (asyncio.EventLoop) the event loop
    :param shell_cmd: (str) the command used to open the handle (default: /bin/bash)
    """
    loop = loop or asyncio.get_event_loop()
    shell = await create_shell_handle(shell_cmd=shell_cmd, loop=loop)
    tasks = create_io_forwarding_tasks(shell, forks=None, callbacks=None)
    # now we wait for the handle to finish
    done, pending = await asyncio.wait(list(tasks.values()) + [shell.handle.wait()], return_when=asyncio.FIRST_COMPLETED)
    # await send_command(shell, shell_cmd)


# Create a clause if this script is running as main
if __name__ == '__main__':
    asyncio.get_event_loop().run_until_complete(main())
