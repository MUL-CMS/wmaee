

from io import StringIO
from time import sleep


class StringStream(StringIO):
    """
    A class representing a dummy stream, which can be used to write data to and read from it
    """

    def __init__(self, string=''):
        """
        Create a string stream with a initial value
        :param string: (str) the initial value (default: "")
        """
        super(StringStream, self).__init__(initial_value=string)
        self._pos = 0
        self._remaining = 0
        self._length = 0

    def read(self, size=-1):
        """
        Performs a read operation on the StringStream object. Blocks if not enough data is available
        :param size: (int) number of characters to be read from the stream (default: -1)
        :return: (str) data read from the StringStream object
        """
        while self._remaining < size:
            # otherwise we block the python interpreter from doing anything
            sleep(0.05)
        result = super(StringStream, self).read()
        # Increase position, from current position seek( ..., 1)
        result_length = len(result)
        # Increase position, and cosume
        self._pos += result_length
        self._remaining = self._length - self._pos
        self.seek(self._pos)
        return result

    def write(self, s):
        """
        Write data to the StringStream object
        :param s: (str) the data
        """
        write_length = len(s)
        super(StringStream, self).write(s)
        # After write file is at the end
        # Seek from back and make it available
        self._length += write_length
        self._remaining = self._length - self._pos

    def readline(self, size=-1, block=True):
        """
        Reads a line from the StringStream object. Block if not a full line is available
        :param size: (int) number of characters to read (default: -1)
        :param block: (bool) wether to block until  a line is available (default: True)
        :return: (str) the data read from the StringStream object
        """
        if self.tell() != self._pos:
            self.seek(self._pos)
        result = super(StringStream, self).readline()
        result_length = len(result)

        self._pos += result_length
        self._remaining = self._length - self._pos
        # Seek new position
        self.seek(self._pos)
        # if block:
        #    while not result:
        #        result = super(StringStream, self).readline()
        #        sleep(0.025)
        return result

